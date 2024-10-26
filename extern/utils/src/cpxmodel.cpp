/**
 * @file cpxmodel.cpp
 * @brief Implementation of MIPModelI for CPLEX
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2019
 */

#include "utils/cpxmodel.h"
#include <stdexcept>
#include <signal.h>
#include <cstring>

int CPXModel_UserBreak = 0;

static void userSignalBreak(int signum)
{
	CPXModel_UserBreak = 1;
}


static void throwCplexError(CPXCENVptr env, int status)
{
	const unsigned int BUF_SIZE = 4096;
	char errmsg[BUF_SIZE];
	CPXgeterrorstring(env, status, errmsg);
	int trailer = std::strlen(errmsg) - 1;
	if (trailer >= 0) errmsg[trailer] = '\0';
	throw std::runtime_error(errmsg);
}


/* Make a call to a Cplex API function checking its return status */
template<typename Func, typename... Args>
void CPX_CALL(Func cpxfunc, CPXENVptr env, Args&&... args)
{
	int status = cpxfunc(env, std::forward<Args>(args)...);
	if (status)  throwCplexError(env, status);
}


CPXModel::CPXModel()
{
	int status = 0;

	env = CPXopenCPLEX(&status);
	if (status)  throwCplexError(nullptr, status);

	lp = CPXcreateprob(env, &status, "");
	if (status)
	{
		CPXcloseCPLEX(&env);
		throwCplexError(nullptr, status);
	}
}


CPXModel::CPXModel(CPXENVptr _env, CPXLPptr _lp, bool _ownEnv, bool _ownLP) : env(_env), lp(_lp), ownEnv(_ownEnv), ownLP(_ownLP)
{
	assert(env && lp);
}


CPXModel::~CPXModel()
{
	if (restoreSignalHandler) handleCtrlC(false);
	if (ownLP) CPXfreeprob(env, &lp);
	if (ownEnv) CPXcloseCPLEX(&env);
}


/* Read/Write */
void CPXModel::readModel(const std::string& filename)
{
	assert(env && lp);
	CPX_CALL(CPXreadcopyprob, env, lp, filename.c_str(), nullptr);
}


void CPXModel::readParams(const std::string& filename)
{
	assert(env);
	CPX_CALL(CPXreadcopyparam, env, filename.c_str());
}


void CPXModel::writeModel(const std::string& filename, const std::string& format) const
{
	assert(env && lp);
	CPX_CALL(CPXwriteprob, env, lp, filename.c_str(), nullptr);
}


void CPXModel::writeSol(const std::string& filename) const
{
	assert(env && lp);
	CPX_CALL(CPXsolwrite, env, lp, filename.c_str());
}


/* Solve */
void CPXModel::lpopt(char method)
{
	assert(env && lp);
	switch(method)
	{
		case 'S': CPX_CALL(CPXlpopt, env, lp); break;
		case 'P': CPX_CALL(CPXprimopt, env, lp); break;
		case 'D': CPX_CALL(CPXdualopt, env, lp); break;
		case 'B': CPX_CALL(CPXbaropt, env, lp); break;
		default: throw std::runtime_error("Unexpected method for lpopt");
	}
}


void CPXModel::mipopt()
{
	assert(env && lp);
	CPX_CALL(CPXmipopt, env, lp);
}


void CPXModel::presolve()
{
	assert(env && lp);
	CPX_CALL(CPXpresolve, env, lp, CPX_ALG_NONE);
}


void CPXModel::postsolve()
{
	assert(env && lp);
	// no-op in CPLEX
}


std::vector<double> CPXModel::postsolveSolution(const std::vector<double>& preX) const
{
	// uncrush solution
	int n = ncols();
	std::vector<double> origX(n, 0.0);
	CPX_CALL(CPXuncrushx, env, lp, &origX[0], &preX[0]);
	return origX;
}


/* Get solution */
double CPXModel::objval() const
{
	assert(env && lp);
	double ret;
	CPX_CALL(CPXgetobjval, env, lp, &ret);
	return ret;
}


void CPXModel::sol(double* x, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	CPX_CALL(CPXgetx, env, lp, x, first, last);
}


bool CPXModel::isPrimalFeas() const
{
	assert(env && lp);
	int primalFeas = 0;
	CPX_CALL(CPXsolninfo, env, lp, nullptr, nullptr, &primalFeas, nullptr);
	return (primalFeas > 0);
}


/* Parameters */
void CPXModel::handleCtrlC(bool flag)
{
	if (flag)
	{
		CPXModel_UserBreak = 0;
		previousHandler = ::signal(SIGINT, userSignalBreak);
		restoreSignalHandler = true;
		CPX_CALL( CPXsetterminate, env, &CPXModel_UserBreak );
	}
	else
	{
		if (restoreSignalHandler)
		{
			::signal(SIGINT, previousHandler);
			restoreSignalHandler = false;
			CPX_CALL( CPXsetterminate, env, nullptr );
		}
	}
}


bool CPXModel::aborted() const
{
	return CPXModel_UserBreak;
}


void CPXModel::seed(int seed)
{
	assert(env);
	CPX_CALL(CPXsetintparam, env, CPX_PARAM_RANDOMSEED, seed);
}


void CPXModel::logging(bool log)
{
	assert(env);
	if (log)  CPX_CALL(CPXsetintparam, env, CPX_PARAM_SCRIND, CPX_ON);
	else      CPX_CALL(CPXsetintparam, env, CPX_PARAM_SCRIND, CPX_OFF);
}


int CPXModel::intParam(IntParam which) const
{
	assert(env);
	int value;

	switch(which)
	{
		case IntParam::Threads:
			CPX_CALL(CPXgetintparam, env, CPX_PARAM_THREADS, &value);
			break;
		case IntParam::SolutionLimit:
			CPX_CALL(CPXgetintparam, env, CPX_PARAM_INTSOLLIM, &value);
			break;
		case IntParam::NodeLimit:
			CPX_CALL(CPXgetintparam, env, CPX_PARAM_NODELIM, &value);
			break;
		case IntParam::IterLimit:
			CPX_CALL(CPXgetintparam, env, CPX_PARAM_ITLIM, &value);
			break;
		default:
			throw std::runtime_error("Unknown integer parameter");
	}

	return (int)value;
}


void CPXModel::intParam(IntParam which, int value)
{
	assert(env);

	switch(which)
	{
		case IntParam::Threads:
			CPX_CALL(CPXsetintparam, env, CPX_PARAM_THREADS, value);
			break;
		case IntParam::SolutionLimit:
			CPX_CALL(CPXsetintparam, env, CPX_PARAM_INTSOLLIM, value);
			break;
		case IntParam::NodeLimit:
			CPX_CALL(CPXsetintparam, env, CPX_PARAM_NODELIM, value);
			break;
		case IntParam::IterLimit:
			CPX_CALL(CPXsetintparam, env, CPX_PARAM_ITLIM, value);
			break;
		default:
			throw std::runtime_error("Unknown integer parameter");
	}
}


double CPXModel::dblParam(DblParam which) const
{
	assert(env);
	double value;

	switch(which)
	{
		case DblParam::TimeLimit:
			CPX_CALL(CPXgetdblparam, env, CPX_PARAM_TILIM, &value);
			break;
		case DblParam::FeasibilityTolerance:
			CPX_CALL(CPXgetdblparam, env, CPX_PARAM_EPRHS, &value);
			break;
		case DblParam::IntegralityTolerance:
			CPX_CALL(CPXgetdblparam, env, CPX_PARAM_EPINT, &value);
			break;
		default:
			throw std::runtime_error("Unknown double parameter");
	}

	return value;
}


void CPXModel::dblParam(DblParam which, double value)
{
	assert(env);
	switch(which)
	{
		case DblParam::TimeLimit:
			CPX_CALL(CPXsetdblparam, env, CPX_PARAM_TILIM, value);
			break;
		case DblParam::FeasibilityTolerance:
			CPX_CALL(CPXsetdblparam, env, CPX_PARAM_EPRHS, value);
			break;
		case DblParam::IntegralityTolerance:
			CPX_CALL(CPXsetdblparam, env, CPX_PARAM_EPINT, value);
			break;
		default:
			throw std::runtime_error("Unknown double parameter");
	}
}


int CPXModel::intAttr(IntAttr which) const
{
	assert(env && lp);
	int value;

	switch(which)
	{
		case IntAttr::Nodes:
			value = CPXgetnodecnt(env, lp);
			break;
		case IntAttr::NodesLeft:
			value = CPXgetnodeleftcnt(env, lp);
			break;
		case IntAttr::BarrierIterations:
			value = CPXgetbaritcnt(env, lp);
			break;
		case IntAttr::SimplexIterations:
			value = CPXgetitcnt(env, lp);
			break;
		default:
			throw std::runtime_error("Unknown integer attribute");
	}

	return value;
}


double CPXModel::dblAttr(DblAttr which) const
{
	assert(env && lp);
	double value;

	switch(which)
	{
		case DblAttr::MIPDualBound:
			CPX_CALL(CPXgetbestobjval, env, lp, &value);
			break;
		default:
			throw std::runtime_error("Unknown double attribute");
	}

	return value;
}


void CPXModel::intParamInternal(int which, int value)
{
	assert(env);
	CPX_CALL(CPXsetintparam, env, which, value);
}


void CPXModel::dblParamInternal(int which, double value)
{
	assert(env);
	CPX_CALL(CPXsetdblparam, env, which, value);
}


/* Access model data */
int CPXModel::nrows() const
{
	assert(env && lp);
	return CPXgetnumrows(env, lp);
}


int CPXModel::ncols() const
{
	assert(env && lp);
	return CPXgetnumcols(env, lp);
}


int CPXModel::nnz() const
{
	int nnz;
	int tmp = 0;
	CPXgetrows(env, lp, &tmp, nullptr, nullptr, nullptr, 0, &nnz, 0, nrows()-1);
	assert(nnz <= 0);
	return -nnz;
}


double CPXModel::objOffset() const
{
	assert(env && lp);
	double objOffset = 0.0;
	CPX_CALL(CPXgetobjoffset, env, lp, &objOffset);
	return objOffset;
}


ObjSense CPXModel::objSense() const
{
	assert(env && lp);
	int cpxobjsen = CPXgetobjsen(env, lp);
	return (cpxobjsen > 0) ? ObjSense::MIN : ObjSense::MAX;
}


void CPXModel::lbs(double* lb, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	CPX_CALL(CPXgetlb, env, lp, lb, first, last);
}


void CPXModel::ubs(double* ub, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	CPX_CALL(CPXgetub, env, lp, ub, first, last);
}


void CPXModel::objcoefs(double* obj, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	CPX_CALL(CPXgetobj, env, lp, obj, first, last);
}


void CPXModel::ctypes(char* ctype, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	CPX_CALL(CPXgetctype, env, lp, ctype, first, last);
}


void CPXModel::sense(char* sense, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	CPX_CALL(CPXgetsense, env, lp, sense, first, last);
}


void CPXModel::rhs(double* rhs, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	CPX_CALL(CPXgetrhs, env, lp, rhs, first, last);
}


void CPXModel::rowbounds(double* rowLB, double* rowUB, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	int count = last - first + 1;
	std::vector<double> rhs(count);
	std::vector<double> rngval(count);
	std::vector<char> sense(count);
	CPX_CALL(CPXgetrhs, env, lp, rhs.data(), first, last);
	CPX_CALL(CPXgetrngval, env, lp, rngval.data(), first, last);
	CPX_CALL(CPXgetsense, env, lp, sense.data(), first, last);
	// CPLEX assumes a ranged row to be in the interval [rhs, rhs+rngval]
	for (int i = 0; i < count; i++) {
		switch(sense[i]) {
			case 'L':
				assert(rngval[i] == 0.0);
				rowLB[first+i] = -CPX_INFBOUND;
				rowUB[first+i] = rhs[i];
				break;
			case 'E':
				assert(rngval[i] == 0.0);
				rowLB[first+i] = rhs[i];
				rowUB[first+i] = rhs[i];
				break;
			case 'G':
				assert(rngval[i] == 0.0);
				rowLB[first+i] = rhs[i];
				rowUB[first+i] = +CPX_INFBOUND;
				break;
			case 'R':
				assert(rngval[i] >= 0.0);
				rowLB[first+i] = rhs[i];
				rowUB[first+i] = rhs[i] + rngval[i];
				break;
			default:
				throw std::runtime_error("Unknown row sense");
		}
	}
}


void CPXModel::row(int ridx, dominiqs::SparseVector& row, char& sense, double& rhs, double& rngval) const
{
	assert(env && lp);
	assert((ridx >= 0) && (ridx < nrows()));
	// get row nz
	int tmp = 0;
	int size;
	CPXgetrows(env, lp, &tmp, &tmp, 0, 0, 0, &size, ridx, ridx);
	// get coef
	size = -size;
	if (size)
	{
		row.resize(size);
		CPX_CALL(CPXgetrows, env, lp, &tmp, &tmp, row.idx(), row.coef(), size, &tmp, ridx, ridx);
	}
	else
	{
		row.clear();
	}
	// here to correctly handle empty constraints
	// get rhs
	CPX_CALL(CPXgetrhs, env, lp, &rhs, ridx, ridx);
	// get sense
	CPX_CALL(CPXgetsense, env, lp, &sense, ridx, ridx);
	// get rngval
	CPX_CALL(CPXgetrngval, env, lp, &rngval, ridx, ridx);
	// CPLEX treats ranged rows considering the constraint satisfied
	// if the linear expression is in the range [rhs, rhs+rngval].
	// However, we interpret ranged rows differently, and the allowed
	// range is [rhs-rngval,rhs] (in both cases, rngval >= 0).
	// So we might have to update the rhs
	if (sense == 'R')
	{
		assert(rngval >= 0.0);
		rhs += rngval;
	}
}


void CPXModel::rows(dominiqs::SparseMatrix& matrix) const
{
	assert(env && lp);
	// get nnz
	int tmp = 0;
	int m = nrows();
	int size;
	matrix.k = m;
	matrix.U = ncols();
	matrix.beg.resize(m);
	CPXgetrows(env, lp, &tmp, matrix.beg.data(), nullptr, nullptr, 0, &size, 0, m - 1);
	size = -size;
	assert(size >= 0);
	matrix.nnz = size;
	// get coefs
	if (size)
	{
		matrix.ind.resize(size);
		matrix.val.resize(size);
		CPX_CALL(CPXgetrows, env, lp, &tmp,
				matrix.beg.data(), matrix.ind.data(), matrix.val.data(),
				size, &tmp, 0, m - 1);
		// fill up cnt
		matrix.cnt.resize(m);
		for (int i = 0; i < (m-1); i++)
		{
			matrix.cnt[i] = matrix.beg[i+1] - matrix.beg[i];
		}
		matrix.cnt[m-1] = matrix.nnz - matrix.beg[m-1];
	}
	else
	{
		matrix.cnt.clear();
		matrix.ind.clear();
		matrix.val.clear();
	}
}


void CPXModel::col(int cidx, dominiqs::SparseVector& col, char& type, double& lb, double& ub, double& obj) const
{
	assert(env && lp);
	assert((cidx >= 0) && (cidx < ncols()));
	// get col nz
	int tmp = 0;
	int size;
	CPXgetcols(env, lp, &tmp, &tmp, 0, 0, 0, &size, cidx, cidx);
	// get coefs
	size = -size;
	if (size)
	{
		col.resize(size);
		CPX_CALL(CPXgetcols, env, lp, &tmp, &tmp, col.idx(), col.coef(), size, &tmp, cidx, cidx);
	}
	else
	{
		col.clear();
	}
	// here to correctly handle empty vars
	// get bounds
	CPX_CALL(CPXgetlb, env, lp, &lb, cidx, cidx );
	CPX_CALL(CPXgetub, env, lp, &ub, cidx, cidx );
	// get obj
	CPX_CALL(CPXgetobj, env, lp, &obj, cidx, cidx );
	// get type
	int status = CPXgetctype(env, lp, &type, cidx, cidx);
	if (status) type = 'C'; // cannot call CPXgetctype on an LP
}


void CPXModel::cols(dominiqs::SparseMatrix& matrix) const
{
	assert(env && lp);
	// get nnz
	int tmp = 0;
	int n = ncols();
	int size;
	matrix.beg.resize(n);
	matrix.k = n;
	matrix.U = nrows();
	CPXgetcols(env, lp, &tmp, matrix.beg.data(), nullptr, nullptr, 0, &size, 0, n - 1);
	size = -size;
	assert(size >= 0);
	matrix.nnz = size;
	// get coefs
	if (size)
	{
		matrix.ind.resize(size);
		matrix.val.resize(size);
		CPX_CALL(CPXgetcols, env, lp, &tmp,
				matrix.beg.data(), matrix.ind.data(), matrix.val.data(),
				size, &tmp, 0, n - 1);
		// fill up cnt
		matrix.cnt.resize(n);
		for (int j = 0; j < (n-1); j++)
		{
			matrix.cnt[j] = matrix.beg[j+1] - matrix.beg[j];
		}
		matrix.cnt[n-1] = matrix.nnz - matrix.beg[n-1];
	}
	else
	{
		matrix.cnt.clear();
		matrix.ind.clear();
		matrix.val.clear();
	}
}


void CPXModel::colNames(std::vector<std::string>& names, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	names.clear();
	int count = last - first + 1;
	std::vector<char> buffer;
	std::vector<char*> cnames(count, nullptr);
	int surplus;
	CPXgetcolname(env, lp, &cnames[0], nullptr, 0, &surplus, first, last);
	if (surplus)
	{
		buffer.resize(-surplus);
		CPX_CALL(CPXgetcolname, env, lp, &cnames[0], &buffer[0], buffer.size(), &surplus, first, last);
		for (int i = 0; i < count; i++) names.push_back(std::string(cnames[i]));
	}
	else
	{
		// no names
		for (int i = 0; i < count; i++) names.push_back("");
	}
}


void CPXModel::rowNames(std::vector<std::string>& names, int first, int last) const
{
	assert(env && lp);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	names.clear();
	int count = last - first + 1;
	std::vector<char> buffer;
	std::vector<char*> rnames(count, 0);
	int surplus;
	CPXgetrowname(env, lp, &rnames[0], nullptr, 0, &surplus, first, last);
	if (surplus)
	{
		buffer.resize(-surplus);
		CPX_CALL(CPXgetrowname, env, lp, &rnames[0], &buffer[0], buffer.size(), &surplus, first, last);
		for (int i = 0; i < count; i++) names.push_back(std::string(rnames[i]));
	}
	else
	{
		// no names
		for (int i = 0; i < count; i++) names.push_back("");
	}
}


/* Data modifications */
void CPXModel::addEmptyCol(const std::string& name, char ctype, double lb, double ub, double obj)
{
	assert(env && lp);
	char* cname = (char*)(name.c_str());
	const char* ctypeptr = (ctype == 'C') ? nullptr : &ctype; //< do not risk turning the model into a MIP
	CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, ctypeptr, &cname);
}


void CPXModel::addCol(const std::string& name, const int* idx, const double* val, int cnt, char ctype, double lb, double ub, double obj)
{
	assert(env && lp);
	int matbeg = 0;
	char* cname = (char*)(name.c_str());
	if (cnt > 0)
	{
		assert(idx && val);
		CPX_CALL(CPXaddcols, env, lp, 1, cnt, &obj, &matbeg, idx, val, &lb, &ub, &cname);
		if (ctype != 'C')
		{
			int last = ncols() - 1;
			CPX_CALL(CPXchgctype, env, lp, 1, &last, &ctype);
		}
	}
	else CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, &ctype, &cname);
}


void CPXModel::addRow(const std::string& name, const int* idx, const double* val, int cnt, char sense, double rhs, double rngval)
{
	assert(env && lp);
	int matbeg = 0;
	char* rname = (char*)(name.c_str());
	if (sense == 'R')
	{
		assert(rngval >= 0.0);
		// for ranged rows, we assue [rhs-rngval,rhs] while CPLEX uses [rhs, rhs+rngval]
		rhs -= rngval;
	}
	CPX_CALL(CPXaddrows, env, lp, 0, 1, cnt, &rhs, &sense, &matbeg, idx, val, nullptr, &rname);
	if (sense == 'R')
	{
		assert(rngval >= 0.0);
		int ridx = nrows()-1;
		assert(ridx >= 0);
		CPX_CALL(CPXchgrngval, env, lp, 1, &ridx, &rngval);
	}
}


void CPXModel::delRow(int ridx)
{
	assert(env && lp);
	CPX_CALL(CPXdelrows, env, lp, ridx, ridx);
}


void CPXModel::delCol(int cidx)
{
	assert(env && lp);
	CPX_CALL(CPXdelcols, env, lp, cidx, cidx);
}


void CPXModel::delRows(int first, int last)
{
	assert(env && lp);
	assert((first >= 0) && (first < nrows()));
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	CPX_CALL(CPXdelrows, env, lp, first, last);
}


void CPXModel::delCols(int first, int last)
{
	assert(env && lp);
	assert((first >= 0) && (first < ncols()));
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	CPX_CALL(CPXdelcols, env, lp, first, last);
}


void CPXModel::objSense(ObjSense objsen)
{
	assert(env && lp);
	CPXchgobjsen(env, lp, static_cast<int>(objsen));
}


void CPXModel::objOffset(double val)
{
	assert(env && lp);
	CPX_CALL(CPXchgobjoffset, env, lp, val);
}


void CPXModel::lb(int cidx, double val)
{
	assert(env && lp);
	assert((cidx >= 0) && (cidx < ncols()));
	char lu = 'L';
	CPX_CALL(CPXchgbds, env, lp, 1, &cidx, &lu, &val);
}


void CPXModel::lbs(int cnt, const int* cols, const double* values)
{
	assert(env && lp);
	std::vector<char> lu(cnt, 'L');
	CPX_CALL(CPXchgbds, env, lp, cnt, cols, &lu[0], values);
}


void CPXModel::ub(int cidx, double val)
{
	assert(env && lp);
	assert((cidx >= 0) && (cidx < ncols()));
	char lu = 'U';
	CPX_CALL(CPXchgbds, env, lp, 1, &cidx, &lu, &val);
}


void CPXModel::ubs(int cnt, const int* cols, const double* values)
{
	assert(env && lp);
	std::vector<char> lu(cnt, 'U');
	CPX_CALL(CPXchgbds, env, lp, cnt, cols, &lu[0], values);
}


void CPXModel::fixCol(int cidx, double val)
{
	assert(env && lp);
	assert((cidx >= 0) && (cidx < ncols()));
	char lu = 'B';
	CPX_CALL(CPXchgbds, env, lp, 1, &cidx, &lu, &val);
}


void CPXModel::objcoef(int cidx, double val)
{
	assert(env && lp);
	assert((cidx >= 0) && (cidx < ncols()));
	CPX_CALL(CPXchgobj, env, lp, 1, &cidx, &val);
}


void CPXModel::objcoefs(int cnt, const int* cols, const double* values)
{
	assert(env && lp);
	CPX_CALL(CPXchgobj, env, lp, cnt, cols, values);
}


void CPXModel::ctype(int cidx, char val)
{
	assert(env && lp);
	assert((cidx >= 0) && (cidx < ncols()));
	assert((val == 'B') || (val == 'I') || (val == 'C'));
	CPX_CALL(CPXchgctype, env, lp, 1, &cidx, &val);
}


void CPXModel::ctypes(int cnt, const int* cols, const char* values)
{
	assert(env && lp);
	CPX_CALL(CPXchgctype, env, lp, cnt, cols, values);
}


void CPXModel::switchToLP()
{
	assert(env && lp);
	CPX_CALL(CPXchgprobtype, env, lp, CPXPROB_LP);
}


/* Private interface */
CPXModel* CPXModel::clone_impl() const
{
	/* We have to create a new environment as well: one of the main reasons
	 * for cloning is to allow for parallel processing of clones, which is
	 * not possible if the env is shared among clones.
	 */
	assert(env && lp);
	int status = 0;

	CPXENVptr newEnv = CPXopenCPLEX(&status);
	if (status)  throwCplexError(nullptr, status);

	CPXLPptr cloned = CPXcloneprob(newEnv, lp, &status);
	if (status) {
		CPXcloseCPLEX(&newEnv);
		throwCplexError(nullptr, status);
	}

	/* For a proper clone, we need to copy nondefault parameters over.
	 * There is not simple API function for that, so we do it by hand...
	 */
	int nchanged = 0;
	int surplus = 0;
	CPXgetchgparam(env, &nchanged, nullptr, 0, &surplus);
	if (surplus < 0) {
		std::vector<int> params(-surplus);
		int ivalue;
		double dvalue;
		char svalue[CPX_STR_PARAM_MAX];
		CPXLONG lvalue;
		CPX_CALL( CPXgetchgparam, env, &nchanged, params.data(), -surplus, &surplus );
		assert(surplus == 0);
		assert(nchanged == params.size());
		for (int param: params) {
			int paramType = CPX_PARAMTYPE_NONE;
			CPX_CALL( CPXgetparamtype, env, param, &paramType );
			switch (paramType) {
				case CPX_PARAMTYPE_NONE:
					/* nothing to do */
					break;
				case CPX_PARAMTYPE_INT:
					CPX_CALL( CPXgetintparam, env, param, &ivalue );
					CPX_CALL( CPXsetintparam, newEnv, param, ivalue );
					break;
				case CPX_PARAMTYPE_DOUBLE:
					CPX_CALL( CPXgetdblparam, env, param, &dvalue );
					CPX_CALL( CPXsetdblparam, newEnv, param, dvalue );
					break;
				case CPX_PARAMTYPE_STRING:
					CPX_CALL( CPXgetstrparam, env, param, svalue );
					CPX_CALL( CPXsetstrparam, newEnv, param, svalue );
					break;
				case CPX_PARAMTYPE_LONG:
					CPX_CALL( CPXgetlongparam, env, param, &lvalue );
					CPX_CALL( CPXsetlongparam, newEnv, param, lvalue );
					break;
				default:
					throw std::runtime_error("Unexpected parameter type");
			}
		}
	}

	return new CPXModel(newEnv, cloned, true, true);
}


CPXModel* CPXModel::presolvedmodel_impl()
{
	int preStat;
	CPX_CALL(CPXgetprestat, env, lp, &preStat, nullptr, nullptr, nullptr, nullptr);

	if (preStat == 2) // reduced to empty problem
	{
		throw std::runtime_error("Problem too easy (solved in presolve stage)");
	}
	else if (preStat == 0) // no presolve
	{
		// nothing to do
		return nullptr;
	}
	else // return clone of presolved problem
	{
		CPXCLPptr redlp = nullptr;
		CPX_CALL(CPXgetredlp, env, lp, &redlp);
		int status = 0;
		CPXLPptr cloned = CPXcloneprob(env, redlp, &status);
		if (status)  throwCplexError(env, status);
		CPXModel* premodel = new CPXModel(env, cloned, false, true);
		return premodel;
	}
	assert( false );
	return nullptr;
}
