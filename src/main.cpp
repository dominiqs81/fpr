/**
 * @file main.cpp
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#include "fixproprep/mip.h"
#include "fixproprep/worker.h"
#include "fixproprep/heurmgr.h"
#include "fixproprep/thread_pool.h"
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
#include <thread>
#include <iostream>
#include <fstream>

using namespace dominiqs;

static const unsigned int MIN_NUMBER_OF_INPUTS = 1;


void runAll(MIPData& data, WorkerDataManager& wManager, Params& params)
{
	int n = data.mip.ncols;
	SolutionPool& solpool{data.solpool};

	// Thread Pool
	ThreadPool thpool(params.threads);

	// First LP-free heuristics
	thpool.enqueue([&]() { runHeur("dfs", "badobjcl", data, wManager, params, 0); });
	thpool.enqueue([&]() { runHeur("dfs", "locks2", data, wManager, params, 0); });
	thpool.enqueue([&]() { runHeur("dive", "locks2", data, wManager, params, 0); });
	thpool.enqueue([&]() { runHeur("dfsrep", "locks", data, wManager, params, 0); });
	thpool.enqueue([&]() { runHeur("dfsrep", "badobjcl", data, wManager, params, 0); });
	thpool.enqueue([&]() { runHeur("diveprop", "random", data, wManager, params, 0); });
	thpool.wait();

	consoleInfo("LP-free heuristics done after {}s [sols={} feas={}]",
			gStopWatch().getElapsed(),
			solpool.getSols().size(),
			solpool.hasFeas());

	if (gStopWatch().getElapsed() >= params.timeLimit)  return;
	if (solpool.hasFeas())  return;

	// Increase the budget after LP-free strategies
	params.maxWorkRatio *= 10;

	// Now we go for zero-objective LPs
	// drop objective
	std::vector<double> zeros(n, 0.0);
	std::vector<int> allIdx(n, 0);
	std::iota(allIdx.begin(), allIdx.end(), 0);
	data.lp->objcoefs(n, allIdx.data(), zeros.data());

	// compute deterministic timelimit for LP solves
	double maxTicks = 1e-2 * data.lp->nnz();
	maxTicks = std::max(maxTicks, 1e3);
	maxTicks = std::min(maxTicks, 1e5);
	consoleInfo("Deterministic LP budget: {} ticks", maxTicks);

	// solve with barrier without crossover
	consoleInfo("Compute zero obj corepoint");
	data.lp->intParamInternal(CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
	data.lp->dblParamInternal(CPXPARAM_DetTimeLimit, maxTicks);
	data.zeroobj_corepoint = solveLP(data.lp, params, 'B', true);

	// restore obj
	data.lp->objcoefs(n, allIdx.data(), data.mip.obj.data());

	if (gStopWatch().getElapsed() >= params.timeLimit)  return;

	if (!data.zeroobj_corepoint.empty()) {
		// run those that depend on core zeroobj solutions
		thpool.enqueue([&]() { runHeur("dfs", "zerocore", data, wManager, params, 0); });
		thpool.enqueue([&]() { runHeur("dive", "zerocore", data, wManager, params, 0); });
		thpool.enqueue([&]() { runHeur("diveprop", "zerocore", data, wManager, params, 0); });
		thpool.wait();

		consoleInfo("Zero-obj corepoint (basic) heuristics done after {}s [sols={} feas={}]",
				gStopWatch().getElapsed(),
				solpool.getSols().size(),
				solpool.hasFeas());

		if (gStopWatch().getElapsed() >= params.timeLimit)  return;

		// Decide whether we run the radomized DFs based on cliques
		bool cliqueStruct = (data.nBinaries == data.mip.ncols) &&
			(data.cliquecover.nCovered() >= (0.5 * data.nBinaries)) &&
			((10 * data.cliquecover.nCliques()) <= data.cliquecover.nCovered());
		if (!solpool.hasFeas() && cliqueStruct) {
			// the number of tries depends on how big the model is...
			int maxTries = params.maxTries;
			if (data.lp->nnz() > 1'000'000)  maxTries /= 10;
			int trial = 0;
			while (trial < maxTries) {
				thpool.enqueue([&]() { runHeur("dfs", "cliques", data, wManager, params, trial + 0); });
				thpool.enqueue([&]() { runHeur("dfs", "cliques", data, wManager, params, trial + 1); });
				thpool.enqueue([&]() { runHeur("dfs", "cliques", data, wManager, params, trial + 2); });
				thpool.enqueue([&]() { runHeur("dfs", "cliques", data, wManager, params, trial + 3); });
				thpool.wait();
				trial += 4;
				if (solpool.hasFeas() || UserBreak)  break;
				if (gStopWatch().getElapsed() >= params.timeLimit)  break;
			}
		}

		consoleInfo("Zero-obj corepoint (advanced) done after {}s [sols={} feas={}]",
				gStopWatch().getElapsed(),
				solpool.getSols().size(),
				solpool.hasFeas());
	}
	else {
		consoleInfo("Could not solve zeroobj-corepoint LP: skipping runs depending on it");
	}

	if (gStopWatch().getElapsed() >= params.timeLimit)  return;
	if (solpool.hasFeas())  return;

	// now enable crossover to get a vertex
	// drop objective (again)
	data.lp->objcoefs(n, allIdx.data(), zeros.data());

	// solve with concurrent LP (should we force primal??)
	consoleInfo("Compute zero obj vertex");
	data.lp->intParamInternal(CPXPARAM_SolutionType, CPX_BASIC_SOLN);
	data.lp->intParamInternal(CPXPARAM_Advance, 0); //< disable warmstart
	data.zeroobj_xlp = solveLP(data.lp, params, 'S', true);
	
	// restore parameters to default
	data.lp->intParamInternal(CPXPARAM_SolutionType, CPX_AUTO_SOLN);
	data.lp->intParamInternal(CPXPARAM_Advance, 1);

	// restore obj
	data.lp->objcoefs(n, allIdx.data(), data.mip.obj.data());
	if (gStopWatch().getElapsed() >= params.timeLimit)  return;
	
	if (!data.zeroobj_xlp.empty()) {
		// finally run those that depend on zeroobj solutions
		thpool.enqueue([&]() { runHeur("dfs", "zerolp", data, wManager, params, 0); });
		thpool.enqueue([&]() { runHeur("diveprop", "zerolp", data, wManager, params, 0); });
		thpool.enqueue([&]() { runHeur("diveprop", "cliques2", data, wManager, params, 0); });
		thpool.wait();
		consoleInfo("Zero-obj LP heuristics done after {}s [sols={} feas={}]",
				gStopWatch().getElapsed(),
				solpool.getSols().size(),
				solpool.hasFeas());
	}
	else {
		consoleInfo("Could not solve zeroobj LP: skipping runs depending on it");
	}

	if (gStopWatch().getElapsed() >= params.timeLimit)  return;
	if (solpool.hasFeas())  return;

	/* Now try LP with original objective */
	data.xlp = solveLP(data.lp, params, 'S', true);
	if (gStopWatch().getElapsed() >= params.timeLimit)  return;
	if (!data.xlp.empty()) {
		// finally run those that depend on zeroobj solutions
		thpool.enqueue([&]() { runHeur("dfs", "lp", data, wManager, params, 0); });
		thpool.enqueue([&]() { runHeur("dive", "lp", data, wManager, params, 0); });
		thpool.enqueue([&]() { runHeur("diveprop", "lp", data, wManager, params, 0); });
		thpool.wait();
		consoleInfo("LP heuristics done after {}s [sols={} feas={}]",
				gStopWatch().getElapsed(),
				solpool.getSols().size(),
				solpool.hasFeas());
	}
	else {
		consoleInfo("Could not original LP: skipping runs depending on it");
	}
}


class MyApp : public App
{
protected:
	// helpers
	bool checkUsage()
	{
		if (args.input.size() < MIN_NUMBER_OF_INPUTS)
		{
			consoleError("usage: mh22 instance_file");
			consoleLog("");
			consoleLog("\tAdditional parameters (e.g., a time limit) can be passed to the command line as key=value pairs,");
			consoleLog("\tsee the list parameters printed to the console at the very beginning of an execution ([config] section).");
			consoleLog("\tExample: mh22 instance_file timeLimit=10 saveSol=1");
			return false;
		}
		return true;
	}
	void exec()
	{
		// read params
		Params params;
		params.readConfig();

		// log config
		consoleInfo("[config]");
		consoleLog("gitHash = {}", FIXPROPREP_GIT_HASH);
		consoleLog("probFile = {}", args.input[0]);
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
			// presolve made no reduction: just clone the original model
			premodel = model->clone();
			hasPresolvedModel = false;
		}
		double presolveTime = gStopWatch().getElapsed();
		consoleLog("presolvedProblem: #rows={} #cols={} #nnz={}", premodel->nrows(), premodel->ncols(), premodel->nnz());
		consoleInfo("Presolve done after {}s", presolveTime);
		assert( premodel );

		if (!params.paramFile.empty())  premodel->readParams(params.paramFile);
		if (params.savePresolved)  premodel->writeModel(getProbName(args.input[0]) + "_pre.mps");
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

		// Worker Data Manager
		bool forceStop = false;
		WorkerDataManager wManager{data, engine, forceStop};

		// run all heuristics
		params.enableOutput = false;
		if (data.nIntegers >= (0.5*data.mip.ncols))  params.maxConsecutiveInfeas = 1.0;
		runAll(data, wManager, params);

		// Print solution pool
		data.solpool.print();

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

		consoleInfo("[results]");
		consoleLog("presolveTime = {}", presolveTime);
		consoleLog("found = {}", (int)data.solpool.hasFeas());
		consoleLog("primalBound = {}", data.solpool.primalBound());
		consoleLog("dualBound = {}", data.dualBound);
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
		consoleLog("maxWork = {}", wManager.maxWork);
		consoleLog("time = {}", gStopWatch().getElapsed());
	}
};


int main (int argc, char const *argv[])
{
	MyApp theApp;
	theApp.parseArgsAndConfig(argc, argv);
	return theApp.run();
}
