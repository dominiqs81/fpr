/**
 * @file main.cpp
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#include "fixproprep/mip.h"
#include "utils/mip.h"
#include <utils/app.h>
#include <utils/floats.h>
#include <utils/cpxmodel.h>
#include <utils/consolelog.h>
#include <utils/path.h>
#include <iostream>
#include <fstream>

using namespace dominiqs;

static const unsigned int MIN_NUMBER_OF_INPUTS = 2;
static const double FEAS_TOL = 1e-5;
static const double INT_TOL = 1e-4;


std::vector<double> readSolution(const MIPInstance& mip, const std::string& solFile)
{
	// create varName -> index map
	std::map<std::string, int> name2idx;
	for (int j = 0; j < mip.ncols; j++)  name2idx[mip.cNames[j]] = j;

	// read solution from file
	std::vector<double> x(mip.ncols, 0.0);
	std::ifstream in(solFile);
	for (int j = 0; j < mip.ncols; j++) {
		std::string var;
		double value;
		in >> var >> value;
		x[name2idx[var]] = value;
	}

	return x;
}


bool checkFeasible(const MIPInstance& mip, std::span<const double> x)
{
	assert( x.size() >= mip.ncols );
	double boundViol = 0.0;
	double intViol = 0.0;
	double rowViol = 0.0;

	// first check bounds and integrality
	for (int j = 0; j < mip.ncols; j++) {
		if (lessThan(x[j], mip.lb[j], FEAS_TOL)) {
			consoleError("Lower bound violation on {}: {} < {}", mip.cNames[j], x[j], mip.lb[j]);
			boundViol = std::max(boundViol, mip.lb[j] - x[j]);
		}
		if (greaterThan(x[j], mip.ub[j], FEAS_TOL)) {
			consoleError("Upper bound violation on {}: {} > {}", mip.cNames[j], x[j], mip.ub[j]);
			boundViol = std::max(boundViol, x[j] - mip.ub[j]);
		}
		if ((mip.xtype[j] != 'C') && !isInteger(x[j], INT_TOL)) {
			consoleError("Integrality violation on {}: {}", mip.cNames[j], x[j]);
			intViol = std::max(intViol, integralityViolation(x[j], INT_TOL)); 
		}
	}

	// check for constraints
	for (int i = 0; i < mip.nRows; i++) {
		// compute violation
		double act = dotProduct(mip.rows[i], x.data());
		double viol = dominiqs::rowViol(act, act, mip.sense[i], mip.rlb[i], mip.rub[i]);

		if (isPositive(viol, FEAS_TOL)) {
			consoleError("Row {} violated by {}", mip.rNames[i], viol);
			rowViol = std::max(rowViol, viol); 
		}
	}

	consoleLog("Max bound violation  : {}", boundViol);
	consoleLog("Max integer violation: {}", intViol);
	consoleLog("Max row violation    : {}", rowViol);

	return (boundViol <= FEAS_TOL) && (intViol <= INT_TOL) && (rowViol <= FEAS_TOL);
}


class MyApp : public App
{
protected:
	bool checkUsage()
	{
		if (args.input.size() < MIN_NUMBER_OF_INPUTS)
		{
			std::cout << "usage: solchecker instance_file sol_file" << std::endl;
			return false;
		}
		return true;
	}
	void exec()
	{
		// read problem
		MIPModelPtr model = MIPModelPtr(new CPXModel());
		model->readModel(args.input[0]);
		consoleLog("originalProblem: #rows={} #cols={} #nnz={}", model->nrows(), model->ncols(), model->nnz());
		MIPInstance mip = extract(model);

		// read solution
		std::vector<double> x = readSolution(mip, args.input[1]);


		// check feasiblity
		bool isFeas = checkFeasible(mip, x);
		if (!isFeas) {
			consoleLog("Solution infeasible");
		}
		else {
			consoleLog("Solution feasible with objective: {}", evalObj(mip, x));
		}
	}
};


int main (int argc, char const *argv[])
{
	MyApp theApp;
	theApp.parseArgsAndConfig(argc, argv);
	return theApp.run();
}
