/**
 * @file dfs.h
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2020 Domenico Salvagnin
 */

#ifndef DFS_SEARCH_H
#define DFS_SEARCH_H

//#define DEBUG_LOG
#ifdef DEBUG_LOG
static int DEBUG_LEVEL = 3;
#endif //< DEBUG_LOG

#include "mip.h"
#include "propagation.h"
#include "repair.h"
#include <utils/app.h>
#include <utils/consolelog.h>
#include <utils/timer.h>
#include <utils/it_display.h>
#include <numeric>
#include <iostream>

using namespace dominiqs;


/** Branch data structure (simplified bound change) */
enum class BranchType { Variable = 0, Clique = 1 };

struct Branch
{
public:
	Branch() = default;
	Branch(BranchType t, int i, char s, double b) : type(t), index(i), sense(s), bound(b) {}
	BranchType type = BranchType::Variable;
	int index = -1; /**< column/clique index */
	char sense; /**< 'L','U','B' for lower, upper, both respectively */
	double bound; /**< new bound */
};


/* Node data structure */
struct Node
{
public:
	Node(Branch b, PropagationEngine::iterator _tp, size_t _d) : branch{b}, trailp(_tp), depth(_d) {}
	Branch branch;
	PropagationEngine::iterator trailp;
	size_t depth;
};


/* Policy Class to customize DFS behaviour */
class DFSStrategy
{
public:
	virtual ~DFSStrategy() {}
	virtual std::vector<Branch> branch(const PropagationEngine& engine, bool nodeInfeas, const Branch& oldBranch) = 0;
};


/* Helper class that applies a branch to the engine */
static inline void applyBranch(const MIPData& data, const Branch& branch, PropagationEngine& engine)
{
	assert( branch.index != -1 );

	if (branch.type == BranchType::Variable) {
		// standard variable branch
		if (branch.sense == 'U') {
			bool nodeInfeas = engine.changeUpperBound(branch.index, branch.bound);
			assert( !nodeInfeas );
		}
		else if (branch.sense == 'L') {
			bool nodeInfeas = engine.changeLowerBound(branch.index, branch.bound);
			assert( !nodeInfeas );
		}
		else {
			bool nodeInfeas = engine.changeLowerBound(branch.index, branch.bound);
			assert( !nodeInfeas );
			nodeInfeas = engine.changeUpperBound(branch.index, branch.bound);
			assert( !nodeInfeas );
		}
	}
	else {
		// clique cover zero branch
		assert(branch.sense == 'B');
		assert(branch.bound == 0.0);
		assert(branch.index >= 0);
		assert(branch.index < data.cliquecover.nCliques());
		bool nodeInfeas = false;
		for (int lit: data.cliquecover.getClique(branch.index)) {
			const auto [var,isPos] = varFromLit(lit, data.mip.ncols);
			if (isPos)  nodeInfeas = engine.changeUpperBound(var, 0.0);
			else        nodeInfeas = engine.changeLowerBound(var, 1.0);
			assert( !nodeInfeas );
		}
	}
}


/** Perform DFS on a given problem: customization of behaviour is provided via StrategyT */
template<typename StrategyT>
void dfsSearch(WorkerDataPtr worker, const Params& params, StrategyT&& strategy)
{
	const MIPData& data{worker->mipdata};
	const MIPInstance& mip{data.mip};
	PropagationEngine& engine{worker->engine};
	PropagationEngine& repairEngine{worker->repairEngine};
	SolutionPool& pool{worker->solpool};
	const Domain& domain = engine.getDomain();
	int n = mip.ncols;
	assert( n == domain.ncols() );
	MIPModelPtr lp = worker->lp;
	assert( lp );
	int consecutiveInfeas = 0;
	int maxConsecutiveInfeas = (int)(params.maxConsecutiveInfeas * n);
	int maxBacktracks = (int)(params.maxBacktracks * n);
	uint64_t maxWork = std::numeric_limits<uint64_t>::max();
	if (params.maxWorkRatio > 0.0) {
		maxWork = uint64_t(lp->nnz() * params.maxWorkRatio);
		maxWork = std::max(maxWork, (uint64_t)1'000'000);
		maxWork = std::min(maxWork, (uint64_t)1'000'000'000);
		consoleLog("maxWork = {}", maxWork);
	}

	/* Effort counting */
	engine.work = 0;
	repairEngine.work = 0;

	std::vector<Node> nodes;
	std::vector<int> allIdx(n);
	std::iota(allIdx.begin(), allIdx.end(), 0);

	RepairMIP repair(data, params);

	dominiqs::IterationDisplay display;
	if (params.enableOutput) {
		display.headerInterval = 10;
		display.iterationInterval = 1;
		display.addColumn("nodes", 0, 8);
		display.addColumn("open", 1, 8);
		display.addColumn("depth", 2, 8);
		display.addColumn("viol", 5, 8);
		display.addColumn("bt", 10, 8);
		display.addColumn("work", 18, 12);
		display.addColumn("time", 20, 10);
		display.printHeader(std::cout);
	}

	// push root node
	PropagationEngine::iterator start_mark = engine.mark();
	nodes.emplace_back(Branch{}, start_mark, 0);
	int nodecnt = 0;
	int nRepair = 0;
	int nBacktracks = 0;
	int lpSolved = 0;
	int numSolutions = 0;
	size_t maxDepth = 0;

	/* Small helpers */
	auto branch2str = [&](const Branch& br) {
		if (br.type == BranchType::Clique)  return fmt::format("sum(clq{}) <= 0", br.index);
		const char* sense;
		if (br.sense == 'B')       sense = "=";
		else if (br.sense == 'L')  sense = ">=";
		else                       sense = "<=";
		return fmt::format("{} {} {}", mip.cNames[br.index], sense, br.bound);
	};

	auto collectEffort = [&]() {
		worker->work += engine.work;
		worker->work += repairEngine.work;
		engine.work = 0;
		repairEngine.work = 0;
	};

	// DFS
	while (!nodes.empty()) {
		// pop node
		Node node = nodes.back();
		nodes.pop_back();
		// backtrack
		engine.undo(node.trailp);
		const Branch& branch = node.branch;
		// engine.debugChecks();
		// apply branch (if any)
		if (branch.index != -1) {
			consoleDebug(2, "Apply branching {}", branch2str(branch));
			applyBranch(data, branch, engine);
		}

		nodecnt++;
		maxDepth = std::max(maxDepth, node.depth);

		bool nodeInfeas = (!engine.violatedRows().empty());

		// propagate & repair
		if (params.propagate) {
			nodeInfeas = engine.propagate(false);
#ifdef DEBUG_LOG
			if (DEBUG_LEVEL >= 3)  printChangesSinceMark(engine, node.trailp);
#endif
			consoleDebug(2, "Bound changes after propagate: {} [infeas={}]", engine.mark() - node.trailp, nodeInfeas);
		}
		if (nodeInfeas && params.repair) {
			if (params.repairSearch)  repair.search(engine, repairEngine);
			else                      repair.walk(engine);
			nodeInfeas = (!engine.violatedRows().empty());
#ifdef DEBUG_LOG
			if (!nodeInfeas && (DEBUG_LEVEL >= 3))  printChangesSinceMark(engine, node.trailp);
#endif
			consoleDebug(2, "Bound changes after repair: {} [infeas={}]", engine.mark() - node.trailp, nodeInfeas);
			nRepair++;
		}

		if (params.enableOutput && ((nodecnt % params.displayInterval) == 0)) {
			collectEffort();
			display.resetIteration();
			display.set("nodes", nodecnt);
			display.set("open", nodes.size());
			display.set("depth", node.depth);
			display.set("viol", engine.violation());
			display.set("bt", nBacktracks);
			display.set("work", worker->work);
			display.set("time", gStopWatch().getElapsed());
			display.printIteration(std::cout);
		}

		if (!engine.violatedRows().empty()) {
			consecutiveInfeas++;
			if (params.backtrackOnInfeas) {
				consoleDebug(3, "Backtrack at depth {}", node.depth);
				nBacktracks++;
				continue;
			}
		}
		else {
			consecutiveInfeas = 0;
		}

		/* branch (other customization point) */
		const Domain& domain = engine.getDomain();
		std::vector<Branch> branches = strategy.branch(engine, nodeInfeas, branch);

		if (branches.empty()) {
			/* End of dive */

			/* If there are continuous variables, we need to compute their values by solving a LP */
			if (data.nContinuous && !nodeInfeas) {
				consoleDebug(3, "Leaf: solving LP");
				// copy bounds from domain to LP
				lp->lbs(allIdx.size(), allIdx.data(), domain.lbs().data());
				lp->ubs(allIdx.size(), allIdx.data(), domain.ubs().data());

				// solve LP
				lp->lpopt('D');
				lpSolved++;

				if (lp->isPrimalFeas()) {
					// looks like we have found a solution
					std::vector<double> x(n);
					consoleDebug(3, "Solution from fixed LP relaxation");
					lp->sol(x.data());
					// copy solution back to engine
					for (int j = 0; j < domain.ncols(); j++) {
						if (!dominiqs::equal(domain.lb(j), domain.ub(j), FEASTOL)) {
							engine.fix(j, x[j]);
						}
						assert( dominiqs::equal(domain.lb(j), domain.ub(j), FEASTOL) );
					}
					assert( engine.violatedRows().empty() );
				}
				else {
					consoleDebug(3, "LP relaxation infeasible");
					nodeInfeas = true;
					nBacktracks++;
				}
			}

			/* Apply 1-opt if feasible*/
			if (!nodeInfeas) {
				repair.oneOpt(engine);
				assert( engine.violatedRows().empty() );
			}
			else if (data.nContinuous) {
				// The solution got infeasible when solving the LP for the continuous variables
				// Hence, the engine does no yet have a complete solution. Make sure it does.
				for (int j = 0; j < domain.ncols(); j++) {
					if (!dominiqs::equal(domain.lb(j), domain.ub(j), FEASTOL)) {
						engine.fix(j, domain.lb(j));
					}
					assert( dominiqs::equal(domain.lb(j), domain.ub(j), FEASTOL) );
				}
			}

			/* Check all variables are fixed */
			for (int j = 0; j < domain.ncols(); j++) {
				assert( dominiqs::equal(domain.lb(j), domain.ub(j), FEASTOL) );
			}

			/* Add solution to pool no matter what */
			std::vector<double> x{domain.lbs().begin(), domain.lbs().end()};
			double objval = evalObj(mip, x);
			bool isFeas = isSolFeasible(mip, x, domain.feasTol);
			SolutionPtr sol = makeFromSpan(mip, x, objval, isFeas, engine.violation());
			sol->timeFound = gStopWatch().getElapsed();
			sol->foundBy = fmt::format("{}_{}", params.run, params.strategy);
			collectEffort();
			sol->work = worker->work;
			pool.add(sol);
			if (isFeas)  numSolutions++;
		}
		else {
			consoleDebug(3, "{}-way branch at depth {}", branches.size(), node.depth);
			// add them in reverse order (this is a stack after all)
			for (auto itr = branches.rbegin(); itr != branches.rend(); ++itr) {
				nodes.emplace_back(*itr, engine.mark(), node.depth+1);
			}
		}

		/* Collect Effort */
		collectEffort();

		/* Termination criteria */
		if ( (nodecnt >= params.maxNodes) ||
			 (lpSolved >= params.maxLpSolved) ||
			 (numSolutions >= params.maxSolutions) ||
			 (worker->work >= maxWork) ||
			 (nBacktracks >= maxBacktracks) ||
			 (consecutiveInfeas >= maxConsecutiveInfeas) ||
			 (gStopWatch().getElapsed() >= params.timeLimit) ||
			 worker->forceStop ||
			 dominiqs::UserBreak ) {
			consoleDebug(2, "Limits reached");
			break;
		}
	}

	// cleanup
	engine.undo(start_mark);
	collectEffort();

	if (params.enableOutput) {
		consoleLog("{} nodes processed: nrepair={} nbacktracks={} work={}", nodecnt, nRepair, nBacktracks, worker->work);
	}

	// make sure we restore the correct bounds on the LP object
	lp->lbs(allIdx.size(), allIdx.data(), mip.lb.data());
	lp->ubs(allIdx.size(), allIdx.data(), mip.ub.data());
}
#endif /* DFS_SEARCH_H */
