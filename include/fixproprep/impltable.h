/**
 * @file   impltable.h
 * @brief  Implication table
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */


#ifndef IMPLTABLE_H
#define IMPLTABLE_H

#include <span>
#include <cassert>
#include <vector>
#include "lit.h"


/** An implication of the form:
 *
 * (x = 0/1) -> (y >=/<= value)
 *
 * Encoding is as follows:
 *
 * Implication               | Lit | Implied
 * -----------------------------------------
 * (x[b] = 1) -> (y[j] >= l) |   b |     j+n
 * (x[b] = 0) -> (y[j] >= l) | b+n |     j+n
 * (x[b] = 1) -> (y[j] <= u) |   b |       j
 * (x[b] = 0) -> (y[j] <= u) | b+n |       j
 * -----------------------------------------
 *
 * So both lit and implied are in [0,2*ncols).
 *
 * In any case, the value stored in BOUND is the actual one (no negation there!).
 */

inline int makeImplied(int var, int ncols, bool isUpper) { return isUpper ? var : var + ncols; }
inline std::pair<int,bool> getImplied(int implied, int ncols)
{
	if (implied < ncols) return {implied, true};
	return {implied-ncols, false};
}

class Implication
{
public:
	Implication() : lit(-1), implied(-1), bound(0.0) {}
	Implication(int l, int y, double v) : lit(l), implied(y), bound(v) {}
	int lit;
	int implied;
	double bound;
};


/** Collection of implications
 *
 * This is used to store implications of the form:
 *
 * (binvar = 0/1) -> (implied >=/<= bound)
 *
 * where binvar is a binary variable and implied is a non-binary variable (implications between binary variables
 * are stored in the cliquetable!).
 *
 */

class ImplTable
{
public:
	/* Getters */
	inline int nImpls() const { return impls.size(); }
	inline int ncols() const { return n; }
	using view_type = std::span<const Implication>;
	inline view_type getImpls() const { return impls; }
	inline view_type getImplsByBin(int binvar, bool isPos) const
	{
		assert((binvar >= 0) && binvar < n);
		assert(sorted);
		int lit = makeLit(binvar, n, isPos);
		return view_type(&impls[firstBin[lit]], &impls[firstBin[lit+1]]); 
	}
	inline view_type getImplsByImplied(int implied, bool isUpper) const
	{
		assert((implied >= 0) && implied < n);
		assert(sorted);
		assert(byImpliedValid);
		assert(byImpliedSorted);
		int var = makeImplied(implied, n, isUpper);
		return view_type(&byImplied[firstImplied[var]], &byImplied[firstImplied[var+1]]); 
	}
	/* Modifiers */
	void clear();
	/* Set the number of columns */
	inline void setNcols(int _n)
	{
		n = _n;
		clear();
	}
	void add(int binvar, bool fixval, int implied, bool isUpper, double bound);
	void sort();
private:
	int n = 0;
	std::vector<Implication> impls;
	bool sorted = false;
	std::vector<size_t> firstBin; //< first implication for each binary
	std::vector<Implication> byImplied; //< copy by implied variable
	bool byImpliedValid = false;
	std::vector<size_t> firstImplied; //< first implication for each implied
	bool byImpliedSorted = false; //< are implications sorted by bound?
};

#endif /* IMPLTABLE_H */
