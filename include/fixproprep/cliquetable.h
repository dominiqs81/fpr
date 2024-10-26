/**
 * @file cliquetable.h
 * @brief Cliquetable
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef CLIQUETABLE_H
#define CLIQUETABLE_H

#include "vectorlist.h"
#include "lit.h"
#include <span>
#include <vector>


/** Collection of cliques
 *
 * Column indices are in the range [0,ncols) for positive literals and [ncols,2ncols) for negative literals.
 */

class CliqueTable
{
public:
	/* Getters */
	inline int nCliques() const { return cliques.nVectors(); }
	inline int ncols() const { return n; }
	inline int nNonzeros() const { return cliques.nNonzeros(); }
	using view_type = VectorList::view_type;
	inline view_type getClique(int cl) const
	{
		assert((cl >= 0) && (cl < cliques.nVectors()));
		return cliques[cl];
	}
	inline bool cliqueIsEqual(int cl) const
	{
		assert((cl >= 0) && (cl < cliques.nVectors()));
		return (type[cl] == 'E');
	}
	inline view_type getLit(int lit) const
	{
		assert(hasLitwise);
		assert(lit >= 0);
		assert(lit < 2*n);

		return literals[lit];
	}
	/* Set the number of columns (this invalidates the litwise repr) */
	inline void setNcols(int _n)
	{
		n = _n;
		cliques.setMaxValue(2*n);
		hasLitwise = false;
	}
	/* Add a clique to the cliquetable */
	void add(std::span<const int> clique, bool isEqual);
	/* Construct literal wise representation */
	void constructLitWiseRepr();
private:
	VectorList cliques; //< clique repr
	std::vector<char> type; //< 'E' for equality cliques, 'L' for the others
	VectorList literals; //< literal repr
	int n = 0; //< number of columns
	bool hasLitwise = false; //< is the literal wise repr valid?
};

#endif /* CLIQUETABLE_H */
