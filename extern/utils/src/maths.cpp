/**
 * @file maths.cpp
 * @brief Mathematical utilities Source
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 * 2008
 */

#include "utils/maths.h"
#include "utils/floats.h"
#include <cstring>

namespace dominiqs {

double Constraint::violation(const double* x) const
{
	double slack = rhs - dotProduct(row, x);
	if (sense == 'L') return -slack;
	if (sense == 'G') return slack;
	if (sense == 'R') return std::max(-slack, slack-range);
	return fabs(slack);
}

void accumulate(double* v, const double* w, int n, double lambda)
{
	for (int i = 0; i < n; ++i) v[i] += (lambda * w[i]);
}

void accumulate(double* v, const int* wIdx, const double* wCoef, int n, double lambda)
{
	for (int i = 0; i < n; ++i) v[wIdx[i]] += (lambda * wCoef[i]);
}

void scale(double* v, int n, double lambda)
{
	for (int i = 0; i < n; ++i) v[i] *= lambda;
}

double dotProduct(const double* x, const double* y, int n)
{
	double ans = 0.0;
	int i;
	if ( n >= 8 )
	{
		for ( i = 0; i < ( n >> 3 ); ++i, x += 8, y += 8 )
			ans += x[0] * y[0] + x[1] * y[1] +
				x[2] * y[2] + x[3] * y[3] +
				x[4] * y[4] + x[5] * y[5] +
				x[6] * y[6] + x[7] * y[7];
		n -= i << 3;
	}
	for ( i = 0; i < n; ++i ) ans += (x[i] * y[i]);
	return ans;
}

double dotProduct(const int* idx, const double* x, int n, const double* y)
{
	double ans = 0.0;
	for (int i = 0; i < n; i++) ans += (x[i] * y[idx[i]]);
	return ans;
}

bool disjoint(const double* x, const double* y, int n)
{
	int cnt = n;
	cnt = cnt >> 1;
	for (int i = cnt; i > 0 ; --i)
	{
		double p1 = fabs(x[0] * y[0]);
		double p2 = fabs(x[1] * y[1]);
		if (isPositive(p1 + p2)) return false;
		x += 2;
		y += 2;
	}
	if (n % 2) return (isPositive(fabs(x[0] * y[0])));
	return true;
}

double euclidianNorm(const double* v, int n)
{
	double ans = 0.0;
	for (int i = 0; i < n; i++) ans += (v[i] * v[i]);
	return sqrt(ans);
}

double euclidianDistance(const double* a1, const double* a2, int n)
{
	double ans = 0.0;
	for (int i = 0; i < n; i++) ans += ((a1[i] - a2[i]) * (a1[i] - a2[i]));
	return sqrt(ans);
}

int lexComp(const double* s1, const double* s2, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (lessThan(s1[i], s2[i])) return -1;
		if (greaterThan(s1[i], s2[i])) return 1;
	}
	return 0;
}


void SparseMatrix::add(SparseMatrix::view_type row)
{
	const int* idx = row.idx();
	const double* coef = row.coef();
	int size = row.size();
	beg.push_back((int)ind.size());
	cnt.push_back(size);
	ind.insert(ind.end(), idx, idx + size);
	val.insert(val.end(), coef, coef + size);
	k++;
	nnz += size;
}


SparseMatrix SparseMatrix::transpose() const
{
	SparseMatrix transposed;
	transposed.k = U;
	transposed.U = k;
	transposed.nnz = nnz;

	if (!U)  return transposed;

	/* compute value counts */
	transposed.cnt.resize(U);
	std::fill(transposed.cnt.begin(), transposed.cnt.end(), 0);
	for (int i = 0; i < k; i++)
	{
		for (const auto& [j,coef]: (*this)[i])
		{
			assert((j >= 0) && (j < U));
			transposed.cnt[j]++;
		}
	}

	/* compute beg */
	transposed.beg.resize(U);
	transposed.beg[0] = 0;
	for (int j = 1; j < U; j++)
	{
		transposed.beg[j] = transposed.beg[j-1] + transposed.cnt[j-1];
	}

	std::vector<int> start = transposed.beg;

	transposed.ind.resize(nnz);
	transposed.val.resize(nnz);

	/* bucket sort */
	for (int i = 0; i < k; i++)
	{
		for (const auto& [j,coef]: (*this)[i])
		{
			transposed.ind[start[j]] = i;
			transposed.val[start[j]] = coef;
			start[j]++;
		}
	}

	/* some checks */
	for (int j = 0; j < (U-1); j++)
	{
		assert(start[j] == transposed.beg[j+1]);
	}
	assert(start[U-1] == transposed.nnz);

	return transposed;
}

} // namespace dominiqs
