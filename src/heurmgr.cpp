/**
 * @brief Heuristic manager
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2022 Domenico Salvagnin
 */

#include "fixproprep/heurmgr.h"
#include "fixproprep/dfs.h"
#include "fixproprep/strategies.h"


static void runDFS(WorkerDataPtr worker, Params params)
{
	const MIPData& data{worker->mipdata};
	const MIPInstance& mip{data.mip};
	PropagationEngine& engine{worker->engine};
	SolutionPool& pool{worker->solpool};
	const Domain& domain = engine.getDomain();
	int n = mip.ncols;
	assert( n == domain.ncols() );

	StaticStrategy strategy{data};

	// Refresh engines
	worker->engine.refresh();
	worker->repairEngine.refresh();

	// Create ranker and value chooser
	RankerPtr ranker = makeRanker(params.varSelect, params, data);
	ValuePtr chooser = makeValueChooser(params.valueSelect, params, data);
	assert( ranker );
	assert( chooser );
	strategy.setup(engine, ranker, chooser);

	PropagationEngine::iterator mark = engine.mark();
	dfsSearch(worker, params, strategy);

	// final repair if still infeasible
	if (params.finalRepair && !pool.hasFeas() && pool.hasSols()) {
		// load solution into engine
		SolutionPtr sol = pool.getSols()[0];
		assert( sol );
		for (int j = 0; j < n; j++)  engine.fix(j, sol->x[j]);
		double oldViol = engine.violation();
		assert( oldViol > 0.0 );

		// repair
		// consoleLog("Final repair attempt");
		params.maxSteps *= 5;
		RepairMIP repair(data, params);
		repair.walk(engine);
		double newViol = engine.violation();
		// collect work from final repair
		worker->work += engine.work;
		engine.work = 0;
		// consoleLog("Final repair: viol {} -> {}", oldViol, newViol);

		// add to pool
		if (newViol < oldViol) {
			consoleLog("Final repair: viol {} -> {}", oldViol, newViol);
			std::vector<double> x{domain.lbs().begin(), domain.lbs().end()};
			double objval = evalObj(mip, x);
			bool isFeas = isSolFeasible(mip, x, domain.feasTol);
			sol = makeFromSpan(mip, x, objval, isFeas, newViol);
			sol->timeFound = gStopWatch().getElapsed();
			sol->work = worker->work;
			sol->foundBy = fmt::format("{}_{}", params.run, params.strategy);
			pool.add(sol);
		}
	}

	engine.undo(mark);
}


void runHeur(
	const std::string& run,
	const std::string& strat,
	MIPData& data,
	WorkerDataManager& wManager,
	Params params,
	int trial)
{
	RunType runType = runTypes[run];
	Strategy strategy = strategies[strat];

	/* Set parameters according to run type and strategy */
	params.run               = run;
	params.propagate         = runType.propagate;
	params.repair            = runType.repair;
	params.backtrackOnInfeas = runType.backtrackOnInfeas;
	params.strategy          = strat;
	params.varSelect         = strategy.varSelect;
	params.valueSelect       = strategy.valueSelect;
	params.seed              = params.seed + 117*trial;
	if (params.maxNodes == -1) {
		params.maxNodes      = data.mip.ncols + 1;
	}
	// Run-specific tunings
	if (params.run == "walkmip") {
		params.maxConsecutiveInfeas = 1.0;
		params.maxNonImpr = 100;
		params.maxSteps = 10000;
		params.repairSearch = false;
		params.repairPropagate = false;
	}

	consoleLog("Heur {}-{} seed={}", run, strat, params.seed);

	// get worker data
	WorkerDataPtr worker = wManager.get();
	PropagationEngine::iterator mark = worker->engine.mark();

	// run DFS
	runDFS(worker, params);

	// release worker data to manager (automatically merges solpool into global one)
	worker->engine.undo(mark);
	wManager.release(worker);
}
