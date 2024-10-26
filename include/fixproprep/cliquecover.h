/**
 * @file cliquecover.h
 * @brief Clique Covers API
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef CLIQUECOVER_H
#define CLIQUECOVER_H

#include "vectorlist.h"
#include "cliquetable.h"
#include "lit.h"
#include <span>
#include <vector>


/** Data structure to store a clique cover for a set of columns.
 *
 * Column indices are in the range [0,ncols) for positive literals and [ncols,2ncols) for negative literals.
 *
 * NOTE: Equality cliques have to added to the cover _before_ non equality cliques!
 */
class CliqueCover
{
public:
	/* Getters */
	inline int nCliques() const { return cliques.nVectors(); }
	inline int ncols() const { return n; }
	inline int nCovered() const { return cliques.nNonzeros(); }
	using view_type = VectorList::view_type;
	inline view_type getClique(int cl) const
	{
		assert((cl >= 0) && (cl < cliques.nVectors()));
		return cliques[cl];
	}
	inline bool cliqueIsEqual(int cl) const
	{
		assert((cl >= 0) && (cl < cliques.nVectors()));
		return (cl < firstNonEq);
	}
	inline std::pair<int,bool> coveredBy(int var) const
	{
		assert((var >= 0) && (var < n));
		assert(hasVarView);
		return {var2clique[var], varNegated[var]};
	}
	inline bool isCovered(int var) const
	{
		assert((var >= 0) && (var < n));
		assert(hasVarView);
		return (var2clique[var] != -1);
	}
	/* Set the number of columns */
	inline void setNcols(int _n)
	{
		n = _n;
		cliques.setMaxValue(2*n);
	}
	/* Add a clique to the clique cover */
	void add(std::span<const int> clique, bool isEqual)
	{
		assert(!clique.empty());
		assert(!isEqual || (nCliques() == firstNonEq));
		cliques.add(clique);
		if (isEqual)  firstNonEq++;
		hasVarView = false;
	}
	/* Construct variable wise representation */
	void constructVarRepr()
	{
		var2clique.resize(n);
		varNegated.resize(n);
		std::fill(var2clique.begin(), var2clique.end(), -1);
		std::fill(varNegated.begin(), varNegated.end(), false);
		for (int cl = 0; cl < nCliques(); cl++) {
			for (int lit: cliques[cl]) {
				const auto [var,isPos] = varFromLit(lit, n);
				var2clique[var] = cl;
				varNegated[var] = !isPos;
			}
		}
		hasVarView = true;
	}
private:
	VectorList cliques; //< sorted list of cliques
	int n = 0; //< number of columns
	int firstNonEq = 0; //< index of first non-equality clique in the cover
	std::vector<int> var2clique; //< map variables to the clique index they are covered by
	std::vector<int> varNegated; //< is a given variable negated in the clique covering it?
	bool hasVarView = false; //< is the variable-wise repr valid?
	
};


/** Construct a greedy clique cover of a given set of variables */
CliqueCover greedyCliqueCover(const CliqueTable& table, const std::vector<int> vars, bool eqOnly = false);

#endif /* CLIQUECOVER_H */
