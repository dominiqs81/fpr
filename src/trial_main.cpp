/**
 * @file trial_main.cpp
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#include "fixproprep/mip.h"
#include "fixproprep/worker.h"
#include "fixproprep/heurmgr.h"
#include "fixproprep/dfs.h"
#include "fixproprep/strategies.h"
#include "fixproprep/table_propagators.h"
#include "fixproprep/linear_propagator.h"
#include "fixproprep/version.h"
#include <utils/app.h>
#include <utils/cpxmodel.h>
#include <utils/consolelog.h>
#include <utils/path.h>
#include <utils/timer.h>
#include <utils/fileconfig.h>
#include <utils/str_utils.h>
#include <iostream>
#include <fstream>

using namespace dominiqs;

static const unsigned int MIN_NUMBER_OF_INPUTS = 1;


bool checkParams(const Params& params)
{
	// Check run type
	if (runTypes.find(params.run) == runTypes.end()) {
		consoleError("Unknown run type: {}", params.run);
		fmt::print("Available run types:");
		for (const auto& itr: runTypes)  fmt::print(" {}", itr.first);
		fmt::print("\n");
		return false;
	}
	// Check strategy type
	if (strategies.find(params.strategy) == strategies.end()) {
		consoleError("Unknown strategy type: {}", params.strategy);
		fmt::print("Available strategies:");
		for (const auto& itr: strategies)  fmt::print(" {}", itr.first);
		fmt::print("\n");
		return false;
	}
	return true;
}


class MyApp : public App
{
protected:
	// helpers
	bool checkUsage()
	{
		if (args.input.size() < MIN_NUMBER_OF_INPUTS)
		{
			consoleError("usage: trymh22 instance_file [run={{dfs,dfsrep,dive,divepro}} strategy={{...}}]");
			consoleLog("");
			consoleLog("\tFor a list of available run types, pass an empty run via: trymh22 run=  instance_file");
			consoleLog("\tFor a list of available strategy, pass an empty strategy via: trymh22 strategy=  instance_file");
			consoleLog("\tAdditional parameters (e.g., a time limit) can be passed to the command line as key=value pairs,");
			consoleLog("\tsee the list parameters printed to the console at the very beginning of an execution ([config] section).");
			consoleLog("\tExample: trymh22 instance_file timeLimit=10 saveSol=1");
			return false;
		}
		return true;
	}
	void exec()
	{
		// read params
		Params params;
		params.readConfig();
		if (!checkParams(params))  return;

		// log config
		consoleInfo("[config]");
		consoleLog("gitHash = {}", FIXPROPREP_GIT_HASH);
		consoleLog("probFile = {}", args.input[0]);
		consoleLog("run = {}", params.run);
		consoleLog("strategy = {}", params.strategy);
		params.logToConsole();
		consoleLog("");

		// read problem
		MIPModelPtr model = MIPModelPtr(new CPXModel());
		model->readModel(args.input[0]);
		consoleLog("originalProblem: #rows={} #cols={} #nnz={}", model->nrows(), model->ncols(), model->nnz());
		if (!params.paramFile.empty())  model->readParams(params.paramFile);

		gStopWatch().start();

		// get original data
		MIPInstance origMip = extract(model);

		// use MIP presolver
		model->presolve();
		MIPModelPtr premodel = model->presolvedModel();
		bool hasPresolvedModel = true;
		if (!params.mipPresolve || !premodel) {
			// presolve disabled or made no reduction: just clone the original model
			premodel = model->clone();
			hasPresolvedModel = false;
		}
		double presolveTime = gStopWatch().getElapsed();
		consoleLog("presolvedProblem: #rows={} #cols={} #nnz={}", premodel->nrows(), premodel->ncols(), premodel->nnz());
		consoleInfo("Presolve done after {}s", presolveTime);
		assert( premodel );

		if (!params.paramFile.empty())  premodel->readParams(params.paramFile);
		if (params.savePresolved)  premodel->writeModel(getProbName(args.input[0]) + "_pre.lp");
		consoleLog("");

		MIPData data(premodel);
		consoleInfo("Data extraction done after {}s", gStopWatch().getElapsed());

		// compute initial dual bound
		data.dualBound = trivialDualBound(data);
		consoleInfo("Trivial dual bound: {}", data.dualBound);

		// initial constraint propagation
		const MIPInstance& mip = data.mip;
		PropagationEngine engine{data};
		engine.add(PropagatorPtr{new CliquesPropagator{data.cliquetable}});
		engine.add(PropagatorPtr{new ImplPropagator{data.impltable}});
		engine.add(PropagatorPtr{new LinearPropagator{data}});
		engine.init(mip.lb, mip.ub, mip.xtype);

		bool infeas = engine.propagate(true);
		assert( !infeas );
		consoleInfo("Initial constraint propagation done after {}s", gStopWatch().getElapsed());

		// Fix trivially roundable binaries once and for all
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.uplocks[j] == 0) {
				bool infeas = engine.changeLowerBound(j, engine.getDomain().ub(j));
				assert( !infeas );
			}
			else if (data.dnlocks[j] == 0) {
				bool infeas = engine.changeUpperBound(j, engine.getDomain().lb(j));
				assert( !infeas );
			}
		}
		// Propagate those fixings, if any
		infeas = engine.propagate(false);
		assert( !infeas );

		// create LP relaxation
		data.lp = premodel->clone();
		data.lp->switchToLP();

		Strategy strategy = strategies[params.strategy];
		int tries = 0;
		double lpTime = 0.0;
		uint64_t maxWork = 0;

		try {
			// construct LP solution if needed by our var/value strategies
			double startTime = gStopWatch().getElapsed();
			int n = data.mip.ncols;
			bool zeroedObj = false;
			bool solvedLP = false;
			if (strategy.depends & DEP_ZEROOBJ) {
				consoleLog("Zeroing out LP objective");
				// first: temporarily drop objective
				std::vector<double> zeros(n, 0.0);
				std::vector<int> allIdx(n, 0);
				std::iota(allIdx.begin(), allIdx.end(), 0);
				data.lp->objcoefs(n, allIdx.data(), zeros.data());
				zeroedObj = true;
			}

			// compute deterministic timelimit for LP solves
			double maxTicks = 1e-2 * data.lp->nnz();
			maxTicks = std::max(maxTicks, 1e3);
			maxTicks = std::min(maxTicks, 1e5);
			consoleInfo("Deterministic LP budget: {} ticks", maxTicks);
			data.lp->dblParamInternal(CPXPARAM_DetTimeLimit, maxTicks);

			if (strategy.depends & DEP_COREPOINT) {
				data.lp->intParamInternal(CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
				if (zeroedObj) {
					data.zeroobj_corepoint = solveLP(data.lp, params, 'B', true);
					if (!data.zeroobj_corepoint.empty())  solvedLP = true;
				}
				else {
					data.corepoint = solveLP(data.lp, params, 'B', true);
					if (!data.corepoint.empty())  solvedLP = true;
				}
			}
			else if (strategy.depends & DEP_VERTEX) {
				data.lp->intParamInternal(CPXPARAM_SolutionType, CPX_BASIC_SOLN);
				if (zeroedObj) {
					data.zeroobj_xlp = solveLP(data.lp, params, 'S', true);
					if (!data.zeroobj_xlp.empty())  solvedLP = true;
				}
				else {
					data.xlp = solveLP(data.lp, params, 'S', true);
					if (!data.xlp.empty())  solvedLP = true;
				}
			}

			if (((strategy.depends & DEP_COREPOINT) || (strategy.depends & DEP_VERTEX)) && !solvedLP) {
				throw std::runtime_error("Could not solve LP to optimality");
			}

			if (solvedLP) {
				// Increase the budget after LP-free strategies
				params.maxWorkRatio *= 10;

				if (zeroedObj) {
					assert( equal(data.lp->objval(), data.mip.objOffset) );
					consoleLog("Zero-obj LP relaxation = {}", data.lp->objval());
				}
				else {
					consoleLog("LP relaxation = {}", data.lp->objval());
				}
			}
			else {
				consoleLog("No LP solved");
			}
			lpTime = gStopWatch().getElapsed() - startTime;

			// Initialize worker data
			// Worker Data Manager
			bool forceStop = false;
			WorkerDataManager wManager{data, engine, forceStop};
			while (tries < params.maxTries) {
				// run heuristic
				runHeur(params.run, params.strategy, data, wManager, params, tries++);

				// stop at first feasible
				if (data.solpool.hasFeas())  break;

				// stop on Ctrl-C or time limit
				if (UserBreak)  break;
				if (gStopWatch().getElapsed() >= params.timeLimit)  break;
			}

			maxWork = wManager.maxWork;
			// Uncrush whatever we have to original space
			if (data.solpool.hasFeas()) {
				SolutionPtr incumbent = data.solpool.getIncumbent();
				assert( incumbent );
				std::vector<double> x;
				if (hasPresolvedModel) {
					consoleLog("Uncrushing incumbent");
					x = model->postsolveSolution(incumbent->x);
				}
				else x = incumbent->x;

				// double check it is still feasible
				assert( x.size() == origMip.ncols );
				assert( isSolFeasible(origMip, x, FEASTOL) );

				if (params.saveSol) {
					// print solution file
					std::string solFile = getProbName(args.input[0]) + ".sol";
					std::ofstream out(solFile);
					for (int j = 0; j < origMip.ncols; j++) {
						out << fmt::format("{} {:.17g}", origMip.cNames[j], x[j]) << std::endl;
					}
				}
			}
			else {
				consoleInfo("Could not find a feasible solution!");
			}
			// Print solution pool
			data.solpool.print();
		}
		catch(std::exception& e) {
			consoleError(e.what());
		}

		consoleInfo("[results]");
		consoleLog("presolveTime = {}", presolveTime);
		consoleLog("lpTime = {}", lpTime);
		consoleLog("found = {}", (int)data.solpool.hasFeas());
		consoleLog("tries = {}", tries);
		consoleLog("primalBound = {}", data.solpool.primalBound());
		if (data.solpool.hasSols())  consoleLog("minViol = {}", data.solpool.minViolation());
		else                         consoleLog("minViol = {}", 1000000.0);
		if (data.solpool.hasFeas()) {
			SolutionPtr incumbent = data.solpool.getIncumbent();
			consoleLog("solWork = {}", incumbent->work);
			consoleLog("solTime = {}", incumbent->timeFound);
		}
		else {
			consoleLog("solWork = {}", -1);
			consoleLog("solTime = {}", -1);
		}
		consoleLog("maxWork = {}", maxWork);
		consoleLog("time = {}", gStopWatch().getElapsed());
	}
};


int main (int argc, char const *argv[])
{
	MyApp theApp;
	theApp.parseArgsAndConfig(argc, argv);
	return theApp.run();
}
