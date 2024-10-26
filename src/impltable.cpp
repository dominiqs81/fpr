#include "fixproprep/impltable.h"
#include <algorithm>


void ImplTable::clear()
{
	impls.clear();
	byImplied.clear();
}


void ImplTable::add(int binvar, bool value, int implied, bool isUpper, double bound)
{
	int lit = makeLit(binvar, n, value);
	int var = makeImplied(implied, n, isUpper);
	impls.emplace_back(lit, var, bound);
	sorted = false;
	byImpliedValid = false;
	byImpliedSorted = false;
}


void ImplTable::sort()
{
	// first of all, sort by binary variables (bucket sort)
	std::vector<size_t> count(2*n, 0);
	for (const auto& impl: impls)  count[impl.lit]++;

	firstBin.resize(2*n+1);
	std::fill(firstBin.begin(), firstBin.end(), 0);
	for (int lit = 1; lit <= 2*n; lit++)  firstBin[lit] = firstBin[lit-1] + count[lit-1];

	// use by implied vector as temporary storage (it will be rewritten anyway)
	byImplied.resize(nImpls());
	std::vector<size_t> start = firstBin;
	for (const auto& impl: impls)  byImplied[start[impl.lit]++] = impl;
	std::swap(byImplied, impls);

	for (int lit = 0; lit < 2*n; lit++)
	{
		assert(start[lit] == firstBin[lit+1]);
	}
	sorted = true;

	// now fill sorted representation by implied variable (again a bucket sort)
	std::fill(count.begin(), count.end(), 0);
	for (const auto& impl: impls)  count[impl.implied]++;

	firstImplied.resize(2*n+1);
	std::fill(firstImplied.begin(), firstImplied.end(), 0);
	for (int var = 1; var <= 2*n; var++)  firstImplied[var] = firstImplied[var-1] + count[var-1];
	start = firstImplied;
	for (const auto& impl: impls)  byImplied[start[impl.implied]++] = impl;

	for (int var = 0; var < 2*n; var++)
	{
		assert(start[var] == firstImplied[var+1]);
	}
	byImpliedValid = true;

	/* We store the implications in the "by implied" array sorted by bound.
	 * For example, if there are k implications that impose a upper bound u_i
	 * on variable y, those are sorted such that u1 <= u2 <= ... <= uk,
	 * i.e., from stronger to weaker (and similarly for lower bounds).
	 * This is to speed up backward propagation.
	 */
	for (int var = 0; var < 2*n; var++)
	{
		double mult = (var < n) ? +1.0 : -1.0;
		std::sort(&byImplied[firstImplied[var]], &byImplied[firstImplied[var+1]],
			[=](const Implication& x, const Implication& y)
		{
			return (mult*x.bound) < (mult*y.bound);
		});
	}

	for (int var = 0; var < 2*n; var++)
	{
		double mult = (var < n) ? +1.0 : -1.0;
		for (int i = firstImplied[var]+1; i < firstImplied[var+1]; i++)
		{
			assert((mult*byImplied[i-1].bound) <= (mult*byImplied[i].bound));
		}
	}
	byImpliedSorted = true;
}
