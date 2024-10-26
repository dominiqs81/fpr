/**
 * @brief MIP related data structures and methods
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#ifndef MIP_H
#define MIP_H

#include "cliquetable.h"
#include "impltable.h"
#include "cliquecover.h"
#include <utils/maths.h>
#include <utils/mipmodel.h>
#include <utils/mip.h>
#include <vector>

const double INFTY = 1e20;
const double FEASTOL = 1e-6;
const double ZEROTOL = 1e-9;


/* Stores a complete MIP solution (feasible or not) */
struct Solution
{
public:
	std::vector<double> x;
	double objval;
	bool isFeas;
	double violation; //< if infeasible, reports the maximum violation
	double timeFound;
	std::string foundBy; //< which algorithm found this solution
	uint64_t work; //< deterministic effort counter
	// equality operator
	bool operator==(const Solution& other) const {
		if (!dominiqs::equal(objval, other.objval))  return false;
		if (!dominiqs::equal(violation, other.violation))  return false;
		if (isFeas != other.isFeas)  return false;
		if (x.size() != other.x.size())  return false;
		for (int j = 0; j < x.size(); j++) {
			if (!dominiqs::equal(x[j], other.x[j]))  return false;
		}
		return true;
	}
};
using SolutionPtr = std::shared_ptr<Solution>;


/* Construction a solution from a span range */
SolutionPtr makeFromSpan(const dominiqs::MIPInstance& mip, std::span<const double> x, double objval, bool isFeas = true, double violation = 0.0);

/* Stores a pool of MIP solutions (feasible or not).
 *
 * If not empty, the first solution is always the best one
 */
class SolutionPool
{
public:
	void setObjSense(double sense)
	{
		assert( (sense == 1.0) || (sense == -1.0) );
		objSense = sense;
	}
	void add(SolutionPtr sol);
	const std::vector<SolutionPtr>& getSols() const { return pool; }
	bool hasFeas() const
	{
		return ((!pool.empty()) && pool[0]->isFeas);
	}
	bool hasSols() const
	{
		return (!pool.empty());
	}
	SolutionPtr getIncumbent() const
	{
		return hasFeas() ? pool[0] : SolutionPtr();
	}
	double primalBound() const
	{
		return hasFeas() ? pool[0]->objval : objSense*INFTY;
	}
	double minViolation() const
	{
		assert( !pool.empty() );
		return pool[0]->violation;
	}
	void merge(SolutionPool& other)
	{
		for (SolutionPtr sol: other.pool)  add(sol);
		other.pool.clear();
		assert( !other.hasSols() );
	}
	void print() const;
protected:
	double objSense = 1.0;
	std::vector<SolutionPtr> pool;
};


/* sort variables within each row by type ('B','I','C') and within each
 * subset by non-increasing absolute value of coefficients.
 */
void normalizeRows(dominiqs::MIPInstance& mip);


/* Row classification */
enum RowClass {
	CLIQUE = 0,
	CLIQUE_EQ = 1,
	CLIQUE_EQ_N = 2,
	SETCOVER = 3,
	CARD = 4,
	CARD_EQ = 5,
	DOUBLE_AGGR = 6,
	VBOUND = 7,
	KNAPSACK = 8,
	KNAPSACK_EQ = 9,
	GENERIC = 10,
	NCLASSES
};

/* Get name of a row type */
inline std::string rClassName(RowClass rc)
{
	switch (rc)
	{
		case CLIQUE:      return "CLIQUE";
		case CLIQUE_EQ:   return "CLIQUE_EQ";
		case CLIQUE_EQ_N: return "CLIQUE_EQ_NEGATED";
		case SETCOVER:    return "SETCOVER";
		case CARD:        return "CARD";
		case CARD_EQ:     return "CARD_EQ";
		case VBOUND:      return "VBOUND";
		case DOUBLE_AGGR: return "DOUBLE_AGGR";
		case KNAPSACK:    return "KNAPSACK";
		case KNAPSACK_EQ: return "KNAPSACK_EQ";
		case GENERIC:     return "GENERIC";
		default:          return "(unknown)";
	}
	return "(unknown)";
}


/* Parameters */
struct Params
{
public:
	// global params
	uint64_t seed = 20220201;
	bool mipPresolve = true;
	bool saveSol = false;
	bool savePresolved = false;
	double timeLimit = 600;
	int threads = 4;
	int maxTries = 100;
	std::string paramFile;
	// dfs params
	bool propagate = true;
	bool repair = false;
	bool backtrackOnInfeas = true;
	double maxConsecutiveInfeas = 0.1; //< node limit as fraction of variables
	double maxBacktracks = 0.2; //< again as a fraction of variables
	// dfs limits
	int maxNodes = -1;
	int maxLpSolved = 10;
	int maxSolutions = 1;
	double maxWorkRatio = 100.0;
	// repair params
	double p = 0.75; //< random walk probability
	int maxNonImpr = 10;
	// repair limits
	int maxSteps = 100;
	bool finalRepair = true;
	bool repairSearch = true;
	bool repairPropagate = true;
	// strategies
	std::string run = "dfs";
	std::string strategy = "goodobj";
	std::string varSelect = "type";
	std::string valueSelect = "random";
	// output
	bool enableOutput = true;
	int displayInterval = 500;
	// read/log
	void readConfig();
	void logToConsole();
};


/* Global data structures on a MIP instance */
struct MIPData
{
	MIPData(MIPModelPtr model);
	// MIP instance
	dominiqs::MIPInstance mip;
	// global structures
	std::vector<RowClass> rclass; //< row classification
	std::vector<int> uplocks; //< up locks
	std::vector<int> dnlocks; //< down locks
	CliqueTable cliquetable;
	ImplTable impltable;
	CliqueCover cliquecover;
	// stats
	int nBinaries;
	int nIntegers;
	int nContinuous;
	int numSingletons;
	int objSupport;
	// solution and bounds
	double dualBound;
	SolutionPool solpool;
	// relaxations
	MIPModelPtr lp; //< LP relaxation model
	std::vector<double> zeroobj_corepoint; //< zero-obj corepoint
	std::vector<double> zeroobj_xlp; //< vertex of the zero-obj LP relxation
	std::vector<double> corepoint; //< optimal LP face corepoint
	std::vector<double> xlp; //< optimal LP relaxation vertex
};


/* classify rows depending on structure
 *
 * @note: In general, there can ambiguity for uniform binary rows:
 *        For example, x1 + x2 <= 1 can be either classified as a clique constraint
 *        or as a set covering constraint on the negated literals
 *        (1-x1) + (1-x2) >= 1.
 *        Similarly, x1 + x2 + x3 <= 1 is a cardinality constraint, but also
 *        a set covering on negated literals (1-x1) + (1-x2) + (1-x3) >= 1.
 */
RowClass classifyRow(dominiqs::SparseVector::view_type row, char sense, double rlb, double rub, std::span<const char> xtype);

void rowClassification(MIPData& data);

// compute variable locks
void computeColLocks(MIPData& data);

// column statistics
void colStats(MIPData& data);

// Extract cliques from instance
void constructCliquetable(MIPData& data);

// Extract implications from instance
void constructImpltable(MIPData& data);

// Construct a cliquecover over the binary variables
void constructCliqueCover(MIPData& data);

// Compute a trivial dual bound based on variable bounds and (possibly) a clique cover
double trivialDualBound(const MIPData& data);

// Solve LP relaxation
std::vector<double> solveLP(MIPModelPtr model, const Params& params, char method, bool enableOutput);


/* Debugging utilities */
void printLocks(const MIPData& data);
void printCliqueTable(const MIPData& data);
void printCliqueCover(const MIPData& data);

#endif /* MIP_H */
