/**
 * @file xprsmodel.cpp
 * @brief Implementation of MIPModelI for XPRESS
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * 2019
 */

#include "utils/xprsmodel.h"
#include <signal.h>
#include <cstring>
#include <climits>
#include <numeric>
#include <fmt/format.h>


XPRSprob probSignalHandler = nullptr;

static void userSignalBreak(int signum)
{
	if (probSignalHandler)  XPRSinterrupt(probSignalHandler, XPRS_STOP_CTRLC);
}


static void throwXpressError(XPRSprob prob)
{
	const unsigned int BUF_SIZE = 512;
	char errmsg[BUF_SIZE];
	XPRSgetlasterror(prob, errmsg);
	throw std::runtime_error(errmsg);
}


/* Make a call to a XPRESS API function checking its return status */
template<typename Func, typename... Args>
void XPRS_CALL(Func xprsfunc, XPRSprob prob, Args&&... args)
{
	int status = xprsfunc(prob, std::forward<Args>(args)...);
	if (status)  throwXpressError(prob);
}


XPRSModel::XPRSModel()
{
	int status = XPRSinit(NULL);
	if (status)  throw std::runtime_error("Cannot initialize XPRESS library");

	status = XPRScreateprob(&prob);
	if (status)
	{
		XPRSfree();
		throw std::runtime_error("Cannot create XPRESS problem");
	}

	/* immediately disable log */
	XPRS_CALL(XPRSsetintcontrol, prob, XPRS_OUTPUTLOG, 0);

	// Oddly enough, this is not yet an empty problem we can add rows/cols to...
	int beg = 0;
	status = XPRSloadlp(prob, "prob", 0, 0, nullptr, nullptr, nullptr, nullptr,
						&beg, nullptr, nullptr, nullptr, nullptr, nullptr);
	if (status)
	{
		XPRSdestroyprob(prob);
		XPRSfree();
		throw std::runtime_error("Cannot create XPRESS problem");
	}
}


XPRSModel::XPRSModel(XPRSprob _prob, bool _ownProb) : prob(_prob), ownProb(_ownProb)
{
	assert(prob);
}


XPRSModel::~XPRSModel()
{
	if (restoreSignalHandler) handleCtrlC(false);
	if (ownProb)
	{
		XPRSdestroyprob(prob);
		XPRSfree();
	}
}


/* Read/Write */
void XPRSModel::readModel(const std::string& filename)
{
	assert(prob);
	XPRS_CALL(XPRSreadprob, prob, filename.c_str(), "");
}


void XPRSModel::readParams(const std::string& filename)
{
	assert(prob);
	throw std::runtime_error("Reading params from file unsupported in XPRESS");
}


void XPRSModel::writeModel(const std::string& filename, const std::string& format) const
{
	assert(prob);
	XPRS_CALL(XPRSwriteprob, prob, filename.c_str(), format.c_str());
}


void XPRSModel::writeSol(const std::string& filename) const
{
	assert(prob);
	XPRS_CALL(XPRSwritebinsol, prob, filename.c_str(), "");
}


/* Solve */
void XPRSModel::lpopt(char method)
{
	assert(prob);
	switch(method)
	{
		case 'S': XPRS_CALL(XPRSlpoptimize, prob, "pdn"); break;
		case 'P': XPRS_CALL(XPRSlpoptimize, prob, "p"); break;
		case 'D': XPRS_CALL(XPRSlpoptimize, prob, "d"); break;
		case 'B': XPRS_CALL(XPRSlpoptimize, prob, "b"); break;
		default: throw std::runtime_error("Unexpected method for lpopt");
	}
}


void XPRSModel::mipopt()
{
	assert(prob);
	XPRS_CALL(XPRSmipoptimize, prob, "");
}


void XPRSModel::presolve()
{
	assert(prob);

	// presolve ~=~ call mipopt with an iteration limit of zero
	// disable reductions that can introduce strange global entities like
	// semicontinuous or ranged rows
	int oldLpIterLimit;
	XPRS_CALL(XPRSgetintcontrol, prob, XPRS_LPITERLIMIT, &oldLpIterLimit);
	XPRS_CALL(XPRSsetintcontrol, prob, XPRS_LPITERLIMIT, 0);
	int presolveOps;
	XPRS_CALL(XPRSgetintcontrol, prob, XPRS_PRESOLVEOPS, &presolveOps);
	int oldPresolveOps = presolveOps;
	presolveOps |= (1 << 10); //< No semi-continuous variable detection
	presolveOps |= (1 << 15); //< No integer variable and SOS detection
	XPRS_CALL(XPRSsetintcontrol, prob, XPRS_PRESOLVEOPS, presolveOps);
	XPRS_CALL(XPRSmipoptimize, prob, "d");
	// restore controls
	XPRS_CALL(XPRSsetintcontrol, prob, XPRS_LPITERLIMIT, oldLpIterLimit);
	XPRS_CALL(XPRSsetintcontrol, prob, XPRS_PRESOLVEOPS, oldPresolveOps);
}


void XPRSModel::postsolve()
{
	assert(prob);
	XPRS_CALL(XPRSpostsolve, prob);
}


std::vector<double> XPRSModel::postsolveSolution(const std::vector<double>& preX) const
{
	// XPRESS does not provide an API to uncrush a single solution
	// so we copy the current problem (which is supposed to be in a presolved state)
	// fix all variables to the values in preX and the optimize (should be fast)
	// to get a solution in the original space
	std::unique_ptr<XPRSModel> copy(clone_impl());
	int n = copy->ncols();
	assert( n == (int)preX.size());
	std::vector<char> xtype(n);
	copy->ctypes(&xtype[0]);
	for (int j = 0; j < n; j++)
	{
		if (xtype[j] != 'C')  copy->fixCol(j, preX[j]);
	}
	copy->mipopt();
	int mipStatus = 0;
	XPRS_CALL(XPRSgetintattrib, copy->prob, XPRS_MIPSTATUS, &mipStatus);
	if (mipStatus != XPRS_MIP_OPTIMAL)
	{
		throw std::runtime_error(fmt::format("Unexpected mipstatus {} in postsolve", mipStatus));
	}
	n = copy->ncols(); //< after a successful lpopt we are no longer in presolve stage!
	std::vector<double> x(n);
	copy->sol(&x[0]);
	return x;
}


/* Get solution */
double XPRSModel::objval() const
{
	assert(prob);
	double ret;
	XPRS_CALL(XPRSgetdblattrib, prob, XPRS_LPOBJVAL, &ret);
	return ret;
}


void XPRSModel::sol(double* x, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));

	int nglents = 0;
	int nsos = 0;
	XPRS_CALL(XPRSgetglobal, prob, &nglents, &nsos, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
	int isMIP = (nglents > 0) || (nsos > 0);

	// XPRESS allows only querying the full solution: if this is also
	// what the user want, we don't need to allocate any temporary memory
	if ((first == 0) && (last == (ncols()-1)))
	{
		if (isMIP)  XPRS_CALL(XPRSgetmipsol, prob, x, nullptr);
		else        XPRS_CALL(XPRSgetlpsol, prob, x, nullptr, nullptr, nullptr);
	}
	else
	{
		// Otherwise we need to store the solution into a temporary buffer
		// and then copy the subsequence the user asked for
		std::vector<double> sol(ncols());
		if (isMIP)  XPRS_CALL(XPRSgetmipsol, prob, &sol[0], nullptr);
		else        XPRS_CALL(XPRSgetlpsol, prob, &sol[0], nullptr, nullptr, nullptr);
		std::copy(sol.begin() + first, sol.begin() + last + 1, x); 
	}
}


bool XPRSModel::isPrimalFeas() const
{
	assert(prob);
	int primalFeas = 0;
	int lpstat = 0;
	XPRS_CALL(XPRSgetintattrib, prob, XPRS_LPSTATUS, &lpstat);
	return (lpstat == XPRS_LP_OPTIMAL);
}


/* Parameters */
void XPRSModel::handleCtrlC(bool flag)
{
	if (flag)
	{
		probSignalHandler = prob;
		previousHandler = ::signal(SIGINT, userSignalBreak);
		restoreSignalHandler = true;
	}
	else
	{
		if (restoreSignalHandler)
		{
			::signal(SIGINT, previousHandler);
			restoreSignalHandler = false;
			probSignalHandler = nullptr;
		}
	}
}


bool XPRSModel::aborted() const
{
	int stopstatus;
	XPRS_CALL(XPRSgetintattrib, prob, XPRS_STOPSTATUS, &stopstatus);
	return (stopstatus == XPRS_STOP_CTRLC);
}


void XPRSModel::seed(int seed)
{
	assert(prob);
	XPRS_CALL(XPRSsetintcontrol, prob, XPRS_RANDOMSEED, seed);
}


void XPRSModel::logging(bool log)
{
	assert(prob);
	if (log)  XPRS_CALL(XPRSsetintcontrol, prob, XPRS_OUTPUTLOG, 1);
	else      XPRS_CALL(XPRSsetintcontrol, prob, XPRS_OUTPUTLOG, 0);
}


int XPRSModel::intParam(IntParam which) const
{
	assert(prob);
	int value;

	switch(which)
	{
		case IntParam::Threads:
			XPRS_CALL(XPRSgetintcontrol, prob, XPRS_THREADS, &value);
			break;
		case IntParam::SolutionLimit:
			XPRS_CALL(XPRSgetintcontrol, prob, XPRS_MAXMIPSOL, &value);
			break;
		case IntParam::NodeLimit:
			XPRS_CALL(XPRSgetintcontrol, prob, XPRS_MAXNODE, &value);
			break;
		case IntParam::IterLimit:
			XPRS_CALL(XPRSgetintcontrol, prob, XPRS_LPITERLIMIT, &value);
			break;
		default:
			throw std::runtime_error("Unknown integer parameter");
	}

	return (int)value;
}


void XPRSModel::intParam(IntParam which, int value)
{
	assert(prob);

	switch(which)
	{
		case IntParam::Threads:
			XPRS_CALL(XPRSsetintcontrol, prob, XPRS_THREADS, value);
			break;
		case IntParam::SolutionLimit:
			XPRS_CALL(XPRSsetintcontrol, prob, XPRS_MAXMIPSOL, value);
			break;
		case IntParam::NodeLimit:
			XPRS_CALL(XPRSsetintcontrol, prob, XPRS_MAXNODE, value);
			break;
		case IntParam::IterLimit:
			XPRS_CALL(XPRSsetintcontrol, prob, XPRS_LPITERLIMIT, value);
			break;
		default:
			throw std::runtime_error("Unknown integer parameter");
	}
}


double XPRSModel::dblParam(DblParam which) const
{
	assert(prob);
	double value;

	switch(which)
	{
		case DblParam::TimeLimit:
			int xprmaxtime;
			XPRS_CALL(XPRSgetintcontrol, prob, XPRS_MAXTIME, &xprmaxtime);
			value = (double)abs(xprmaxtime);
			break;
		case DblParam::FeasibilityTolerance:
			XPRS_CALL(XPRSgetdblcontrol, prob, XPRS_FEASTOL, &value);
			break;
		case DblParam::IntegralityTolerance:
			XPRS_CALL(XPRSgetdblcontrol, prob, XPRS_MIPTOL, &value);
			break;
		default:
			throw std::runtime_error("Unknown double parameter");
	}

	return value;
}


void XPRSModel::dblParam(DblParam which, double value)
{
	assert(prob);
	switch(which)
	{
		case DblParam::TimeLimit:
			int xprmaxtime;
			if (value >= INT_MAX)  xprmaxtime = 0; //< 0 means no time limit
			else                   xprmaxtime = (int) floor(value);
			XPRS_CALL(XPRSsetintcontrol, prob, XPRS_MAXTIME, -xprmaxtime);
			break;
		case DblParam::FeasibilityTolerance:
			XPRS_CALL(XPRSsetdblcontrol, prob, XPRS_FEASTOL, value);
			break;
		case DblParam::IntegralityTolerance:
			XPRS_CALL(XPRSsetdblcontrol, prob, XPRS_MIPTOL, value);
			break;
		default:
			throw std::runtime_error("Unknown double parameter");
	}
}


int XPRSModel::intAttr(IntAttr which) const
{
	assert(prob);
	int value;

	switch(which)
	{
		case IntAttr::Nodes:
			XPRS_CALL(XPRSgetintattrib, prob, XPRS_NODES, &value);
			break;
		case IntAttr::NodesLeft:
			XPRS_CALL(XPRSgetintattrib, prob, XPRS_ACTIVENODES, &value);
			break;
		case IntAttr::BarrierIterations:
			XPRS_CALL(XPRSgetintattrib, prob, XPRS_BARITER, &value);
			break;
		case IntAttr::SimplexIterations:
			XPRS_CALL(XPRSgetintattrib, prob, XPRS_SIMPLEXITER, &value);
			break;
		default:
			throw std::runtime_error("Unknown integer attribute");
	}

	return value;
}


double XPRSModel::dblAttr(DblAttr which) const
{
	assert(prob);
	double value;

	switch(which)
	{
		case DblAttr::MIPDualBound:
			XPRS_CALL(XPRSgetdblattrib, prob, XPRS_BESTBOUND, &value);
			break;
		default:
			throw std::runtime_error("Unknown double attribute");
	}

	return value;
}


void XPRSModel::intParamInternal(int which, int value)
{
	assert(prob);
	XPRS_CALL(XPRSsetintcontrol, prob, which, value);
}


void XPRSModel::dblParamInternal(int which, double value)
{
	assert(prob);
	XPRS_CALL(XPRSsetdblcontrol, prob, which, value);
}


/* Access model data */
int XPRSModel::nrows() const
{
	assert(prob);
	int ret;
	XPRS_CALL(XPRSgetintattrib, prob, XPRS_ROWS, &ret);
	return ret;
}


int XPRSModel::ncols() const
{
	assert(prob);
	int ret;
	XPRS_CALL(XPRSgetintattrib, prob, XPRS_COLS, &ret);
	return ret;
}


int XPRSModel::nnz() const
{
	assert(prob);
	int ret;
	XPRS_CALL(XPRSgetintattrib, prob, XPRS_ELEMS, &ret);
	return ret;
}


double XPRSModel::objOffset() const
{
	assert(prob);
	double ret;
	XPRS_CALL(XPRSgetdblattrib, prob, XPRS_OBJRHS, &ret);
	return ret;
}


ObjSense XPRSModel::objSense() const
{
	assert(prob);
	double objsen;
	XPRS_CALL(XPRSgetdblattrib, prob, XPRS_OBJSENSE, &objsen);
	return (objsen > 0) ? ObjSense::MIN : ObjSense::MAX;
}


void XPRSModel::lbs(double* lb, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	XPRS_CALL(XPRSgetlb, prob, lb, first, last);
}


void XPRSModel::ubs(double* ub, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	XPRS_CALL(XPRSgetub, prob, ub, first, last);
}


void XPRSModel::objcoefs(double* obj, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	XPRS_CALL(XPRSgetobj, prob, obj, first, last);
}


void XPRSModel::ctypes(char* ctype, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	XPRS_CALL(XPRSgetcoltype, prob, ctype, first, last);
}


void XPRSModel::sense(char* sense, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	XPRS_CALL(XPRSgetrowtype, prob, sense, first, last);
}


void XPRSModel::rhs(double* rhs, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	XPRS_CALL(XPRSgetrhs, prob, rhs, first, last);
}


void XPRSModel::rowbounds(double* rowLB, double* rowUB, int first, int last) const
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
	XPRS_CALL(XPRSgetrhs, prob, rhs.data(), first, last);
	XPRS_CALL(XPRSgetrngval, prob, rngval.data(), first, last);
	XPRS_CALL(XPRSgetrhsrange, prob, sense.data(), first, last);
	// XPRESS assumes a ranged row to be in the interval [rhs-rngval, rhs]
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
				rowLB[first+i] = rhs[i] - rngval[i];
				rowUB[first+i] = rhs[i];
				break;
			default:
				throw std::runtime_error("Unknown row sense");
		}
	}
}


void XPRSModel::row(int ridx, dominiqs::SparseVector& row, char& sense, double& rhs, double& rngval) const
{
	assert(prob);
	assert((ridx >= 0) && (ridx < nrows()));
	// get row nz
	int size;
	XPRSgetrows(prob, nullptr, nullptr, nullptr, 0, &size, ridx, ridx);
	// get coef
	if (size)
	{
		row.resize(size);
		int matbeg = 0;
		int rcnt = 0;
		XPRS_CALL(XPRSgetrows, prob, &matbeg, row.idx(), row.coef(), size, &rcnt, ridx, ridx);
		assert(rcnt == size);
	}
	else
	{
		row.clear();
	}
	// here to correctly handle empty constraints
	// get rhs
	XPRS_CALL(XPRSgetrhs, prob, &rhs, ridx, ridx);
	// get sense
	XPRS_CALL(XPRSgetrowtype, prob, &sense, ridx, ridx);
	// get rngval
	XPRS_CALL(XPRSgetrhsrange, prob, &rngval, ridx, ridx);
	assert(rngval >= 0.0);
}


void XPRSModel::rows(dominiqs::SparseMatrix& matrix) const
{
	assert(prob);
	// get nnz
	int m = nrows();
	int size;
	matrix.k = m;
	matrix.U = ncols();
	matrix.beg.resize(m);
	XPRSgetrows(prob, nullptr, nullptr, nullptr, 0, &size, 0, m - 1);
	assert(size >= 0);
	matrix.nnz = size;
	// get actual data
	if (size)
	{
		matrix.ind.resize(size);
		matrix.val.resize(size);
		int nnz;
		XPRS_CALL(XPRSgetrows, prob,
				matrix.beg.data(), matrix.ind.data(), matrix.val.data(),
				size, &nnz, 0, m - 1);
		assert(nnz == size);
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


void XPRSModel::col(int cidx, dominiqs::SparseVector& col, char& type, double& lb, double& ub, double& obj) const
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	// get col nz
	int size;
	XPRSgetcols(prob, nullptr, nullptr, nullptr, 0, &size, cidx, cidx);
	// get coef
	if (size)
	{
		col.resize(size);
		int matbeg = 0;
		int cnt = 0;
		XPRS_CALL(XPRSgetcols, prob, &matbeg, col.idx(), col.coef(), size, &cnt, cidx, cidx);
		assert(cnt == size);
	}
	else
	{
		col.clear();
	}
	// here to correctly handle empty vars
	// get bounds
	XPRS_CALL(XPRSgetlb, prob, &lb, cidx, cidx );
	XPRS_CALL(XPRSgetub, prob, &ub, cidx, cidx );
	// get obj
	XPRS_CALL(XPRSgetobj, prob, &obj, cidx, cidx );
	// get type
	int status = XPRSgetcoltype(prob, &type, cidx, cidx);
	if (status) type = 'C'; // cannot call XPRSgetcoltype on an LP
}


void XPRSModel::cols(dominiqs::SparseMatrix& matrix) const
{
	assert(prob);
	// get nnz
	int n = ncols();
	int size;
	matrix.beg.resize(n);
	matrix.k = n;
	matrix.U = nrows();
	XPRSgetcols(prob, nullptr, nullptr, nullptr, 0, &size, 0, n - 1);
	assert(size >= 0);
	matrix.nnz = size;
	// get actual data
	if (size)
	{
		matrix.ind.resize(size);
		matrix.val.resize(size);
		int nnz;
		XPRS_CALL(XPRSgetcols, prob,
				matrix.beg.data(), matrix.ind.data(), matrix.val.data(),
				size, &nnz, 0, n - 1);
		assert(nnz == size);
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


void XPRSModel::colNames(std::vector<std::string>& names, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	if (last == -1)  last = ncols()-1;
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	names.clear();
	int count = last - first + 1;
	// get number of bytes required to store the names
	int buffer_len;
    XPRSgetnamelist(prob, 2, nullptr, 0, &buffer_len, first, last);
	if (buffer_len)
	{
		// get names
		std::vector<char> buffer(buffer_len);
		XPRS_CALL(XPRSgetnamelist, prob, 2, &buffer[0], buffer_len, nullptr, first, last);
		int offset = 0;
		for (int i = 0; i < count; i++)
		{
			const char* name = &buffer[offset];
			names.emplace_back(name);
			offset += std::strlen(name) + 1;
		}
	}
	else
	{
		// no names
		for (int i = 0; i < count; i++) names.push_back("");
	}
}


void XPRSModel::rowNames(std::vector<std::string>& names, int first, int last) const
{
	assert(prob);
	assert((first >= 0) && (first < nrows()));
	if (last == -1)  last = nrows()-1;
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	names.clear();
	int count = last - first + 1;
	// get number of bytes required to store the names
	int buffer_len;
    XPRSgetnamelist(prob, 1, nullptr, 0, &buffer_len, first, last);
	if (buffer_len)
	{
		// get names
		std::vector<char> buffer(buffer_len);
		XPRS_CALL(XPRSgetnamelist, prob, 1, &buffer[0], buffer_len, nullptr, first, last);
		int offset = 0;
		for (int i = 0; i < count; i++)
		{
			const char* name = &buffer[offset];
			names.emplace_back(name);
			offset += std::strlen(name) + 1;
		}
	}
	else
	{
		// no names
		for (int i = 0; i < count; i++) names.push_back("");
	}
}


/* Data modifications */
void XPRSModel::addEmptyCol(const std::string& name, char ctype, double lb, double ub, double obj)
{
	assert(prob);
	int matbeg = 0;
	XPRS_CALL(XPRSaddcols, prob, 1, 0, &obj, &matbeg, nullptr, nullptr, &lb, &ub);

	int cidx = ncols()-1;
	if (ctype != 'C')
	{
		XPRS_CALL(XPRSchgcoltype, prob, 1, &cidx, &ctype);
	}

	char* cname = (char*)(name.c_str());
	XPRS_CALL(XPRSaddnames, prob, 2, cname, cidx, cidx);
}


void XPRSModel::addCol(const std::string& name, const int* idx, const double* val, int cnt, char ctype, double lb, double ub, double obj)
{
	assert(prob);
	assert(cnt && idx && val);

	int matbeg = 0;
	XPRS_CALL(XPRSaddcols, prob, 1, cnt, &obj, &matbeg, idx, val, &lb, &ub);

	int cidx = ncols()-1;
	if (ctype != 'C')
	{
		XPRS_CALL(XPRSchgcoltype, prob, 1, &cidx, &ctype);
	}

	char* cname = (char*)(name.c_str());
	XPRS_CALL(XPRSaddnames, prob, 2, cname, cidx, cidx);
}


void XPRSModel::addRow(const std::string& name, const int* idx, const double* val, int cnt, char sense, double rhs, double rngval)
{
	assert(prob);

	int matbeg = 0;
	XPRS_CALL(XPRSaddrows, prob, 1, cnt, &sense, &rhs, &rngval, &matbeg, idx, val);

	int ridx = nrows()-1;
	char* rname = (char*)(name.c_str());
	XPRS_CALL(XPRSaddnames, prob, 1, rname, ridx, ridx);
}


void XPRSModel::delRow(int ridx)
{
	assert(prob);
	assert((ridx >= 0) && (ridx < nrows()));
	XPRS_CALL(XPRSdelrows, prob, 1, &ridx);
}


void XPRSModel::delCol(int cidx)
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	XPRS_CALL(XPRSdelcols, prob, 1, &cidx);
}


void XPRSModel::delRows(int first, int last)
{
	assert(prob);
	assert((first >= 0) && (first < nrows()));
	assert((last >= 0) && (last < nrows()));
	assert(first <= last);
	int count = last - first + 1;
	std::vector<int> idx(count);
	std::iota(idx.begin(), idx.end(), first);
	XPRS_CALL(XPRSdelrows, prob, count, &idx[0]);
}


void XPRSModel::delCols(int first, int last)
{
	assert(prob);
	assert((first >= 0) && (first < ncols()));
	assert((last >= 0) && (last < ncols()));
	assert(first <= last);
	int count = last - first + 1;
	std::vector<int> idx(count);
	std::iota(idx.begin(), idx.end(), first);
	XPRS_CALL(XPRSdelcols, prob, count, &idx[0]);
}


void XPRSModel::objSense(ObjSense objsen)
{
	assert(prob);
	XPRS_CALL(XPRSchgobjsense, prob, (objsen == ObjSense::MIN) ? XPRS_OBJ_MINIMIZE : XPRS_OBJ_MAXIMIZE);
}


void XPRSModel::objOffset(double val)
{
	assert(prob);
	int idx = -1;
	val = -val;
	XPRS_CALL(XPRSchgobj, prob, 1, &idx, &val);
}


void XPRSModel::lb(int cidx, double val)
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	char lu = 'L';
	XPRS_CALL(XPRSchgbounds, prob, 1, &cidx, &lu, &val);
}


void XPRSModel::lbs(int cnt, const int* cols, const double* values)
{
	assert(prob);
	std::vector<char> lu(cnt, 'L');
	XPRS_CALL(XPRSchgbounds, prob, cnt, cols, &lu[0], values);
}


void XPRSModel::ub(int cidx, double val)
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	char lu = 'U';
	XPRS_CALL(XPRSchgbounds, prob, 1, &cidx, &lu, &val);
}


void XPRSModel::ubs(int cnt, const int* cols, const double* values)
{
	assert(prob);
	std::vector<char> lu(cnt, 'U');
	XPRS_CALL(XPRSchgbounds, prob, cnt, cols, &lu[0], values);
}


void XPRSModel::fixCol(int cidx, double val)
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	char lu = 'B';
	XPRS_CALL(XPRSchgbounds, prob, 1, &cidx, &lu, &val);
}


void XPRSModel::objcoef(int cidx, double val)
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	XPRS_CALL(XPRSchgobj, prob, 1, &cidx, &val);
}


void XPRSModel::objcoefs(int cnt, const int* cols, const double* values)
{
	assert(prob);
	XPRS_CALL(XPRSchgobj, prob, cnt, cols, values);
}


void XPRSModel::ctype(int cidx, char val)
{
	assert(prob);
	assert((cidx >= 0) && (cidx < ncols()));
	assert((val == 'B') || (val == 'I') || (val == 'C'));
	XPRS_CALL(XPRSchgcoltype, prob, 1, &cidx, &val);
}


void XPRSModel::ctypes(int cnt, const int* cols, const char* values)
{
	assert(prob);
	XPRS_CALL(XPRSchgcoltype, prob, cnt, cols, values);
}


void XPRSModel::switchToLP()
{
	int n = ncols();
	for (int j = 0; j < n; j++)  ctype(j, 'C');
}


/* Private interface */
XPRSModel* XPRSModel::clone_impl() const
{
	assert(prob);
	std::unique_ptr<XPRSModel> cloned(new XPRSModel());
	XPRS_CALL(XPRScopyprob, cloned->prob, prob, "cloned");
	XPRS_CALL(XPRScopycontrols, cloned->prob, prob);
	return cloned.release();
}


XPRSModel* XPRSModel::presolvedmodel_impl()
{
	assert(prob);
	// get mip status
	int mipStatus = 0;
	XPRS_CALL(XPRSgetintattrib, prob, XPRS_MIPSTATUS, &mipStatus);

	if ((mipStatus == XPRS_MIP_INFEAS) || (mipStatus == XPRS_MIP_OPTIMAL))
	{
		throw std::runtime_error("Problem too easy (solved in presolve stage)");
	}
	else // return clone of current presolved problem
	{
		// we cannot just call XPRScopyprob, as that remembers that this is a
		// presolved model, while we want a fresh copy. So we need to copy
		// the problem data piece by piece...
		int presolveState;
		XPRS_CALL(XPRSgetintattrib, prob, XPRS_PRESOLVESTATE, &presolveState);
		bool isPresolved = presolveState & ((1 << 1) /* LP presolve */ | (1 << 2) /* MIP presolve */);
		assert(isPresolved);

		// get data of the presolved model
		int n = ncols();
		int m = nrows();
		std::vector<double> xlb(n);
		std::vector<double> xub(n);
		std::vector<double> xobj(n);
		std::vector<char> xtype(n);
		std::vector<char> rsense(m);
		std::vector<double> rrhs(m);
		std::vector<double> rngval(m);
		dominiqs::SparseMatrix matrix;
		std::vector<std::string> xNames;
		std::vector<std::string> rNames;
		lbs(&xlb[0]);
		ubs(&xub[0]);
		objcoefs(&xobj[0]);
		ctypes(&xtype[0]);
		sense(&rsense[0]);
		rhs(&rrhs[0]);
		rows(matrix);
		XPRS_CALL(XPRSgetrhsrange, prob, &rngval[0], 0, m-1);
		rowNames(rNames);
		colNames(xNames);

		// create new problem
		std::unique_ptr<XPRSModel> premodel(new XPRSModel());

		// copy everything into it
		for (int j = 0; j < n; j++)
		{
			if ((xtype[j] != 'C') && (xtype[j] != 'I') && (xtype[j] != 'B'))  throw std::runtime_error("Unsupported variable type for FP");
			premodel->addEmptyCol(xNames[j], xtype[j], xlb[j], xub[j], xobj[j]);
		}
		for (int i = 0; i < m; i++)
		{
			int rcnt = matrix.beg[i+1] - matrix.beg[i];
			const int* idx = &matrix.ind[matrix.beg[i]];
			const double* val = &matrix.val[matrix.beg[i]];
			premodel->addRow(rNames[i], idx, val, rcnt, rsense[i], rrhs[i], rngval[i]);
		}
		// don't forget objective offset
		double objoff = objOffset();
		premodel->objOffset(objoff);
		return premodel.release();
	}
	assert( false );
	return nullptr;
}
