/**
 * @brief MIP instance data
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#include "utils/mip.h"

namespace dominiqs {

MIPInstance extract(MIPModelPtr model)
{
	MIPInstance mip;
	mip.ncols = model->ncols();
	mip.nRows = model->nrows();
	// get obj
	mip.objSense = (model->objSense() == ObjSense::MIN) ? 1.0 : -1.0;
	mip.objOffset = model->objOffset();
	mip.obj.resize(mip.ncols);
	model->objcoefs(mip.obj.data());
	// get col data
	mip.xtype.resize(mip.ncols);
	model->ctypes(mip.xtype.data());
	mip.lb.resize(mip.ncols);
	model->lbs(mip.lb.data());
	mip.ub.resize(mip.ncols);
	model->ubs(mip.ub.data());
	model->cols(mip.cols);
	// get row data
	mip.sense.resize(mip.nRows);
	model->sense(mip.sense.data());
	mip.rlb.resize(mip.nRows);
	mip.rub.resize(mip.nRows);
	model->rowbounds(mip.rlb.data(), mip.rub.data());
	model->rows(mip.rows);
	// names
	model->rowNames(mip.rNames);
	model->colNames(mip.cNames);
	return mip;
}


/* Checks whether a given solution vector x is feasible */
bool isSolFeasible(const MIPInstance& mip, std::span<const double> x, double feasTol)
{
	assert( x.size() >= mip.ncols );
	// first check bounds and integrality
	for (int j = 0; j < mip.ncols; j++) {
		if (lessThan(x[j], mip.lb[j], feasTol))
			return false;
		if (greaterThan(x[j], mip.ub[j], feasTol))
			return false;
		if ((mip.xtype[j] != 'C') && !isInteger(x[j], feasTol))
			return false;
	}

	// check for constraints
	for (int i = 0; i < mip.nRows; i++) {
		// compute violation
		double act = dotProduct(mip.rows[i], x.data());
		double viol = rowViol(act, act, mip.sense[i], mip.rlb[i], mip.rub[i]);

		if (isPositive(viol, feasTol)) {
			return false;
		}
	}

	return true;
}


/* Evaluate the objective value of a given solution vector */
double evalObj(const MIPInstance& mip, std::span<const double> x)
{
	assert( (int)x.size() >= mip.ncols );
	return dotProduct(mip.obj.data(), x.data(), mip.ncols) + mip.objOffset;
}

} //< namespace dominiqs
