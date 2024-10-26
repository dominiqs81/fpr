/**
 * @brief MIP instance data
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2023 Domenico Salvagnin
 */

#ifndef MIP_INST_H
#define MIP_INST_H

#include <vector>
#include "maths.h"
#include "mipmodel.h"

namespace dominiqs {

/* MIP instance data */
struct MIPInstance
{
	int ncols = 0;
	int nRows = 0;
	// obj
	double objSense = 1.0;
	double objOffset = 0.0;
	std::vector<double> obj;
	// col data
	std::vector<char> xtype;
	std::vector<double> lb;
	std::vector<double> ub;
	dominiqs::SparseMatrix cols;
	// row data
	std::vector<char> sense;
	std::vector<double> rlb;
	std::vector<double> rub;
	dominiqs::SparseMatrix rows;
	// names (for debugging and output)
	std::vector<std::string> rNames;
	std::vector<std::string> cNames;
};


/* Extract instance data from a MIP model */
MIPInstance extract(MIPModelPtr model);

/* Checks whether a given solution vector x is feasible */
bool isSolFeasible(const MIPInstance& mip, std::span<const double> x, double feasTol);

/* Evaluate the objective value of a given solution vector */
double evalObj(const MIPInstance& mip, std::span<const double> x);

/* Evaluate gap between two bounds (normalized to [0,1], see primal integral paper */
inline double evalGap(double primalBound, double dualBound, double epsZero = 1e-9)
{
	if (dominiqs::isNull(dualBound, epsZero) && dominiqs::isNull(primalBound, epsZero))  return 0.0;
	if ((primalBound * dualBound) < 0.0)  return 1.0;
	return fabs(primalBound - dualBound) / std::max(fabs(primalBound), fabs(dualBound));
}

/* Evaluate violation of constraint from activity bounds */
inline double rowViol(double minAct, double maxAct, char sense, double rLB, double rUB)
{
	double lbviol = (sense != 'L') ? std::max(rLB-maxAct, 0.0) : 0.0;
	double ubviol = (sense != 'G') ? std::max(minAct-rUB, 0.0) : 0.0;
	return std::max(lbviol, ubviol);
}

} //< namespace dominiqs

#endif /* MIP_INST_H */
