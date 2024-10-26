/**
 * @brief MIP related data structures and methods
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#include "fixproprep/mip.h"
#include <utils/consolelog.h>
#include <utils/fileconfig.h>
#include <utils/timer.h>
#include <utils/maths.h>

using namespace dominiqs;


SolutionPtr makeFromSpan(const MIPInstance& mip, std::span<const double> x, double objval, bool isFeas, double violation)
{
	SolutionPtr sol = std::make_shared<Solution>();
	sol->x.insert(sol->x.end(), x.begin(), x.end());
	sol->objval = objval;
	assert( equal(sol->objval, evalObj(mip, sol->x)) );
	sol->isFeas = isFeas;
	assert( sol->isFeas == isSolFeasible(mip, sol->x, FEASTOL) );
	sol->violation = violation;
	assert( (bool)(sol->isFeas) == isNull(sol->violation) );

	return sol;
}


void SolutionPool::add(SolutionPtr sol)
{
	if (!sol)  return;

	// avoid adding duplicates
	auto end = pool.end();
	auto compEq = [&](const SolutionPtr& other) {
		return (*sol) == (*other);
	};
	auto itr = std::find_if(pool.begin(), pool.end(), compEq);
	if (itr != end)  return;

	// feasible solution are kept sorted by objective and infeasible ones by violation
	// feasible solutions always come before infeasible ones
	auto comp = [&](const SolutionPtr& sol1, const SolutionPtr& sol2) {
		if (sol1->isFeas == sol2->isFeas) {
			if (sol1->isFeas)  return (objSense*sol1->objval < objSense*sol2->objval);
			else               return (sol1->violation < sol2->violation);
		}
		return (sol1->isFeas);
	};
	pool.insert( 
		std::upper_bound(pool.begin(), pool.end(), sol, comp),
		sol);

	// 20 solutions should be enough for our purposes
	if (pool.size() > 20)  pool.resize(20);
}


static double solDistance(std::span<const double> x1, std::span<const double> x2)
{
	double ret = 0.0;
	assert(x1.size() == x2.size());
	for (int j = 0; j < x1.size(); j++) {
		ret += fabs(x1[j] - x2[j]);
	}
	return ret;
}


void SolutionPool::print() const
{
	if (pool.empty()) {
		consoleInfo("Empty solution pool");
		return;
	}

	consoleInfo("Solution pool: {} solutions", pool.size());
	consoleLog("{:>8}{:>15}{:>15}{:>7}{:>12}{:>12}{:>12}  {}", "n", "Objective", "Violation", "Feas", "L1 dist", "Time", "Work", "FoundBy");

	SolutionPtr first = pool[0];
	for (int k = 0; k < pool.size(); k++) {
		SolutionPtr sol = pool[k];
		consoleLog("{:>8}{:>15.2f}{:>15.4f}{:>7}{:>12.2f}{:>12.2f}{:>12}  {}",
			k, sol->objval, sol->violation, sol->isFeas, solDistance(first->x, sol->x), sol->timeFound, sol->work, sol->foundBy);
	}
}


/* sort variables in a row by type ('B','I','C') and within each
 * subset by non-increasing absolute value of coefficients.
 */
static void normalizeRow(MIPInstance& mip, int* vars, double* coefs, int count)
{
	/* First bucket sort by type */
	int numBin = 0;
	int numInt = 0;
	for (int k = 0; k < count; k++) {
		if (mip.xtype[vars[k]] == 'B') numBin++;
		else if (mip.xtype[vars[k]] == 'I') numInt++;
	}

	int startBin = 0;
	int startInt = numBin;
	int startCont = numBin + numInt;
	std::vector<std::pair<int,double>> row(count);
	for (int k = 0; k < count; k++) {
		if (mip.xtype[vars[k]] == 'B')      row[startBin++] = { vars[k], coefs[k] };
		else if (mip.xtype[vars[k]] == 'I') row[startInt++] = { vars[k], coefs[k] };
		else                                row[startCont++] = { vars[k], coefs[k] };
	}
	assert(startBin == numBin);
	assert(startInt == (numBin+numInt));

	/* Now sort by abs(coefficient) within each type */
	auto cmp = [](const auto& e1, const auto& e2) {
		return (fabs(e1.second) > fabs(e2.second));
	};
	std::stable_sort(row.begin(), row.begin() + numBin, cmp);
	std::stable_sort(row.begin() + numBin, row.begin() + numBin + numInt, cmp);
	std::stable_sort(row.begin() + numBin + numInt, row.end(), cmp);

	/* Copy over original row */
	for (int k = 0; k < count; k++) {
		vars[k] = row[k].first;
		coefs[k] = row[k].second;
	}

	/* Make sure everything is as expected */
	for (int k = 0; k < (numBin-1); k++) {
		assert(mip.xtype[vars[k]] == 'B');
		assert(fabs(coefs[k]) >= fabs(coefs[k+1]));
	}
	for (int k = numBin; k < (numBin+numInt-1); k++) {
		assert(mip.xtype[vars[k]] == 'I');
		assert(fabs(coefs[k]) >= fabs(coefs[k+1]));
	}
	for (int k = numBin+numInt; k < count-1; k++) {
		assert(mip.xtype[vars[k]] == 'C');
		assert(fabs(coefs[k]) >= fabs(coefs[k+1]));
	}
}


void normalizeRows(MIPInstance& mip)
{
	int m = mip.nRows;
	for (int i = 0; i < m; i++) {
		int* vars = &(mip.rows.ind[mip.rows.beg[i]]);
		double* coefs = &(mip.rows.val[mip.rows.beg[i]]);
		normalizeRow(mip, vars, coefs, mip.rows.cnt[i]);
	}
}


#define READ_PARAM(name)  name = gConfig().get(""#name, defaults.name)

void Params::readConfig()
{
	Params defaults;
	READ_PARAM(seed);
	READ_PARAM(mipPresolve);
	READ_PARAM(saveSol);
	READ_PARAM(savePresolved);
	READ_PARAM(timeLimit);
	READ_PARAM(threads);
	READ_PARAM(paramFile);
	READ_PARAM(maxTries);
	READ_PARAM(run);
	READ_PARAM(strategy);
	READ_PARAM(propagate);
	READ_PARAM(repair);
	READ_PARAM(backtrackOnInfeas);
	READ_PARAM(maxConsecutiveInfeas);
	READ_PARAM(maxBacktracks);
	READ_PARAM(maxNodes);
	READ_PARAM(maxLpSolved);
	READ_PARAM(maxSolutions);
	READ_PARAM(maxWorkRatio);
	READ_PARAM(p);
	READ_PARAM(maxNonImpr);
	READ_PARAM(maxSteps);
	READ_PARAM(finalRepair);
	READ_PARAM(repairSearch);
	READ_PARAM(repairPropagate);
	READ_PARAM(enableOutput);
	READ_PARAM(displayInterval);
}


#define LOG_PARAM(name)  consoleLog(#name" = {}", name)

void Params::logToConsole()
{
	LOG_PARAM(seed);
	LOG_PARAM(mipPresolve);
	LOG_PARAM(saveSol);
	LOG_PARAM(savePresolved);
	LOG_PARAM(timeLimit);
	LOG_PARAM(threads);
	LOG_PARAM(paramFile);
	LOG_PARAM(maxTries);
	LOG_PARAM(run);
	LOG_PARAM(strategy);
	LOG_PARAM(propagate);
	LOG_PARAM(repair);
	LOG_PARAM(backtrackOnInfeas);
	LOG_PARAM(maxConsecutiveInfeas);
	LOG_PARAM(maxBacktracks);
	LOG_PARAM(maxNodes);
	LOG_PARAM(maxLpSolved);
	LOG_PARAM(maxSolutions);
	LOG_PARAM(maxWorkRatio);
	LOG_PARAM(p);
	LOG_PARAM(maxNonImpr);
	LOG_PARAM(maxSteps);
	LOG_PARAM(finalRepair);
	LOG_PARAM(repairSearch);
	LOG_PARAM(repairPropagate);
	LOG_PARAM(enableOutput);
	LOG_PARAM(displayInterval);
}


// classify a single row
RowClass classifyRow(SparseVector::view_type row, char sense, double rlb, double rub, std::span<const char> xtype)
{
	assert( row.size() ); //< no empty rows please!
	// compute stats for this row
	int posbincnt = 0;
	int negbincnt = 0;
	int othercnt = 0;
	double minAbs = INFTY;
	double maxAbs = 0.0;
	for (const auto& [j,v]: row) {
		if (xtype[j] != 'B') {
			othercnt++;
		}
		else {
			double absV = fabs(v);
			minAbs = std::min(minAbs, absV);
			maxAbs = std::max(maxAbs, absV);
			if (v > 0.0)  posbincnt++;
			else          negbincnt++;
		}
	}

	// use them for classification
	RowClass ret = GENERIC;
	// TODO: better management of ranged rows
	if (sense == 'R')  return ret;

	if (othercnt) {
		// non all-binary row
		if ((othercnt == 1) && ((posbincnt+negbincnt)==1)) {
			if (sense == 'E')  ret = DOUBLE_AGGR;
			else               ret= VBOUND;
		}
		else {
			ret = GENERIC;
		}
	}
	else {
		// all binary row
		if (greaterThan(maxAbs, minAbs, ZEROTOL)) {
			// knapsack types
			if (sense == 'E')  ret = KNAPSACK_EQ;
			else               ret = KNAPSACK;
		}
		else {
			// uniform binary row: (generalized) setcover/setpartition/setpacking/cardinality constraints
			assert(equal(minAbs, maxAbs, ZEROTOL));
			double factor = minAbs;
			if (sense == 'L') {
				if (equal(factor*(1 - negbincnt), rub, ZEROTOL))     ret = CLIQUE;
				else if (equal(factor*(posbincnt-1), rub, ZEROTOL))  ret = SETCOVER;
				else                                                 ret = CARD;
			}
			else if (sense == 'G') {
				if (equal(factor*(1 - negbincnt), rlb, ZEROTOL))     ret = SETCOVER;
				else if (equal(factor*(posbincnt-1), rlb, ZEROTOL))  ret = CLIQUE;
				else                                                 ret = CARD;
			}
			else if (sense == 'E') {
				assert(rlb == rub);
				if (equal(factor*(1 - negbincnt), rub, ZEROTOL))     ret = CLIQUE_EQ;
				else if (equal(factor*(posbincnt-1), rub, ZEROTOL))  ret = CLIQUE_EQ_N;
				else                                                 ret = CARD_EQ;
			}
			else {
				assert(false);
			}
		}
	}
	return ret;
}


// classify rows depending on structure
void rowClassification(MIPData& data)
{
	const MIPInstance& mip = data.mip;
	int m = mip.nRows;
	int n = mip.ncols;
	data.rclass.resize(m);

	// classify rows
	for (int i = 0; i < m; i++) {
		data.rclass[i] = classifyRow(mip.rows[i], mip.sense[i], mip.rlb[i], mip.rub[i], mip.xtype);
	}

	// print stats
	std::vector<int> counts(RowClass::NCLASSES, 0);
	for (int i = 0; i < m; i++)  counts[data.rclass[i]]++;
	consoleLog("Row classification:");
	for (int rc = 0; rc < (int)RowClass::NCLASSES; rc++)  consoleLog("{}: {}", rClassName((RowClass)rc), counts[rc]);
	consoleLog("");
}


// compute variable locks
void computeColLocks(MIPData& data)
{
	const MIPInstance& mip = data.mip;
	int n = mip.ncols;
	int m = mip.nRows;
	data.uplocks.resize(n);
	data.dnlocks.resize(n);
	std::fill(data.uplocks.begin(), data.uplocks.end(), 0);
	std::fill(data.dnlocks.begin(), data.dnlocks.end(), 0);

	for (int i = 0; i < m; i++) {
		const auto& row = mip.rows[i];
		for (const auto& [j,v]: row) {
			if (mip.sense[i] == 'E') {
				data.uplocks[j]++;
				data.dnlocks[j]++;
			}
			else {
				double mult = (mip.sense[i] == 'L') ? 1.0 : -1.0;
				if ((mult*v) > 0.0)  data.uplocks[j]++;
				else                 data.dnlocks[j]++;
			}
		}
	}
}


void colStats(MIPData& data)
{
	const MIPInstance& mip = data.mip;
	int numBin = 0;
	int numInt = 0;
	int numCont = 0;
	int numSingletons = 0;
	double minLB = INFTY; //< for integer variables only
	double maxUB = 0.0; //< for integer variables only
	double maxDomain = 0.0; //< for integer variables only
	int dnRoundable = 0;
	int upRoundable = 0;
	int objSupport = 0;

	// get variable data
	int n = mip.ncols;

	// collects stats
	for (int j = 0; j < n; j++) {
		if (isNotNull(mip.obj[j])) objSupport++;

		const auto& col = mip.cols[j];
		if (col.size() <= 1)  numSingletons++;

		if (!data.uplocks[j])  upRoundable++;
		if (!data.dnlocks[j])  dnRoundable++;

		if (mip.xtype[j] == 'B') {
			numBin++;
		}
		else if (mip.xtype[j] == 'I') {
			numInt++;
			minLB = std::min(minLB, mip.lb[j]);
			maxUB = std::max(maxUB, mip.ub[j]);
			maxDomain = std::max(maxDomain, mip.ub[j]-mip.lb[j]);
		}
		else if (mip.xtype[j] == 'C') {
			numCont++;
		}
		else {
			assert(false);
		}
	}

	if (!numInt) minLB = 0.0;

	data.nBinaries = numBin;
	data.nIntegers = numInt;
	data.nContinuous = numCont;
	data.numSingletons = numSingletons;
	data.objSupport = objSupport;

	// print stats
	consoleLog("Col classification:");
	consoleLog("BINARIES: {}", numBin);
	consoleLog("INTEGERS: {} [{},{}] |{}|", numInt, minLB, maxUB, maxDomain);
	consoleLog("CONTINUOUS: {}", numCont);
	consoleLog("SINGLETONS: {}", numSingletons);
	consoleLog("OBJSUPPORT: {}", objSupport);
	consoleLog("UPROUNDABLE: {}", upRoundable);
	consoleLog("DNROUNDABLE: {}", dnRoundable);
	consoleLog("");
}


/* Extract cliques from instance */
void constructCliquetable(MIPData& data) {
	const MIPInstance& mip = data.mip;
	data.cliquetable.setNcols(mip.ncols);
	assert(data.cliquetable.nCliques() == 0);
	assert(data.cliquetable.nNonzeros() == 0);

	int m = mip.nRows;
	int n = mip.ncols;
	std::vector<int> clique;

	for (int i = 0; i < m; i++) {
		if ((data.rclass[i] == CLIQUE) ||
			(data.rclass[i] == CLIQUE_EQ) ||
			(data.rclass[i] == CLIQUE_EQ_N)) {
			assert(mip.sense[i] != 'R');
			// This is a clique row straight from the matrix
			clique.clear();
			double mult = 1.0;
			// We need to figure out whether this is a negated clique or not
			if (data.rclass[i] == CLIQUE) {
				mult = (mip.sense[i] == 'G') ? -1.0 : +1.0;
			}
			else {
				mult = (data.rclass[i] == CLIQUE_EQ) ? 1.0 : -1.0;
			}
			const auto& row = mip.rows[i];
			for (const auto& [j,v]: row)
			{
				if (mult*v > 0.0)  clique.push_back(posLit(j,n));
				else               clique.push_back(negLit(j,n));
			}
			data.cliquetable.add(clique, mip.sense[i] == 'E');
		}
		/* TODO: we could extract cliques from knapsack rows, for example... */
	}

	data.cliquetable.constructLitWiseRepr();

	consoleLog("Cliquetable: {} cliques and {} nonzeros", data.cliquetable.nCliques(), data.cliquetable.nNonzeros());
}


/* Extract implications from instance */
void constructImpltable(MIPData& data) {
	const MIPInstance& mip = data.mip;
	data.impltable.setNcols(mip.ncols);
	assert(data.impltable.nImpls() == 0);

	int m = mip.nRows;
	int n = mip.ncols;

	for (int i = 0; i < m; i++) {
		if (data.rclass[i] == VBOUND) {
			const auto& row = mip.rows[i];
			assert(mip.sense[i] != 'E');
			assert(mip.sense[i] != 'R');
			assert(row.size() == 2);
			const int* vars = row.idx();
			const double* coefs = row.coef();
			assert(mip.xtype[vars[0]] == 'B');
			assert(mip.xtype[vars[1]] != 'B');
			double bincoef = coefs[0];
			double othercoef = coefs[1];
			char sense = mip.sense[i];
			bincoef /= othercoef;
			/* flip sense if we divided by a negative number */
			if (othercoef < 0.0)  sense = (sense == 'L') ? 'G' : 'L';
			if (sense == 'L') {
				double rhs = mip.rub[i] / othercoef;
				// we get an implied upper bound
				double bound = mip.ub[vars[1]];
				// for x = 0 the new upper bound is rhs
				if (lessThan(rhs, bound, ZEROTOL))  data.impltable.add(vars[0], false, vars[1], true, rhs);
				// for x = 1 the new upper bound is rhs-bincoef
				if (lessThan(rhs-bincoef, bound, ZEROTOL))  data.impltable.add(vars[0], true, vars[1], true, rhs-bincoef);
			}
			else {
				double rhs = mip.rub[i] / othercoef;
				// we get an implied lower bound
				double bound = mip.lb[vars[1]];
				// for x = 0 the new lower bound is rhs
				if (greaterThan(rhs, bound, ZEROTOL))  data.impltable.add(vars[0], false, vars[1], false, rhs);
				// for x = 1 the new lower bound is rhs-bincoef
				if (greaterThan(rhs-bincoef, bound, ZEROTOL))  data.impltable.add(vars[0], true, vars[1], false, rhs-bincoef);
			}
		}
		/* TODO: what about double_aggr and other constraints? */
	}

	data.impltable.sort();

	consoleLog("Impltable: {} implications", data.impltable.nImpls());
}


// Construct a cliquecover over the binary variables
void constructCliqueCover(MIPData& data)
{
	// clique covers
	std::vector<int> binaries;
	for (int j = 0; j < data.mip.ncols; j++) {
		if (data.mip.xtype[j] == 'B')  binaries.push_back(j);
	}
	data.cliquecover = greedyCliqueCover(data.cliquetable, binaries, false);

	consoleLog("Clique cover: {} cliques, {} / {}", data.cliquecover.nCliques(), data.cliquecover.nCovered(), binaries.size());
}


// Compute a trivial dual bound based on variable bounds and (possibly) a clique cover
double trivialDualBound(const MIPData& data)
{
	const MIPInstance& mip = data.mip;
	double dualBound = mip.objOffset;

	for (int j = 0; j < mip.ncols; j++) {
		if (mip.obj[j] * mip.objSense > 0.0)  dualBound += mip.obj[j] * mip.lb[j];
		else                                  dualBound += mip.obj[j] * mip.ub[j];
	}

	return dualBound;
}


// Solve LP relaxation
std::vector<double> solveLP(MIPModelPtr model, const Params& params, char method, bool enableOutput)
{
	std::vector<double> x;
	model->logging(enableOutput);
	model->handleCtrlC(true);
	model->seed(params.seed);
	model->dblParam(DblParam::TimeLimit, std::max(params.timeLimit-gStopWatch().getElapsed(), 0.0));
	model->intParam(IntParam::Threads, params.threads);
	model->lpopt(method);
	model->handleCtrlC(false);
	model->logging(false);
	if (model->isPrimalFeas()) {
		x.resize(model->ncols());
		model->sol(x.data());
	}
	return x;
}


// Init MIP Data from a MIP model
MIPData::MIPData(MIPModelPtr model)
{
	mip = extract(model);
	solpool.setObjSense(mip.objSense);
	dualBound = -INFTY * mip.objSense;

	// normalize rows
	normalizeRows(mip);

	// row classification
	rowClassification(*this);

	// variable locks
	computeColLocks(*this);
	colStats(*this);

	// global tables
	constructCliquetable(*this);
	constructImpltable(*this);
	constructCliqueCover(*this);
}


void printLocks(const MIPData& data)
{
	fmt::print("Locks:\n");
	for (int j = 0; j < data.mip.ncols; j++) {
		fmt::print("{}: up = {} dn = {}\n", data.mip.cNames[j], data.uplocks[j], data.dnlocks[j]);
	}
}


void printCliqueTable(const MIPData& data)
{
	const CliqueTable& ct = data.cliquetable;
	fmt::print("Clique table: {} cliques {} nonzeros\n", ct.nCliques(), ct.nNonzeros());
	for (int cl = 0; cl < ct.nCliques(); cl++) {
		for (int lit: ct.getClique(cl)) {
			const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
			fmt::print("{}{} ", isPos ? '+' : '-', data.mip.cNames[j]);
		}
		fmt::print("{} 1\n", ct.cliqueIsEqual(cl) ? ">=" : "==");
	}
}


void printCliqueCover(const MIPData& data)
{
	const CliqueCover& cc = data.cliquecover;
	fmt::print("Clique cover: {} cliques {} / {}\n", cc.nCliques(), cc.nCovered(), data.nBinaries);
	for (int cl = 0; cl < cc.nCliques(); cl++) {
		for (int lit: cc.getClique(cl)) {
			const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
			fmt::print("{}{} ", isPos ? '+' : '-', data.mip.cNames[j]);
		}
		fmt::print("\n");
	}
}
