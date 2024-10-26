/**
 * @file cliquecover.cpp
 * @brief Clique Covers API
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#include "fixproprep/cliquecover.h"
#include <algorithm>
#include <iostream>


/** Construct a greedy cliques cover of a given set of variables */
CliqueCover greedyCliqueCover(const CliqueTable& table, const std::vector<int> vars, bool eqOnly)
{
	int n = table.ncols();
	CliqueCover cc;
	cc.setNcols(n);

	/*************************************
	 * First equality cliques
	 ************************************/

	// extract all equality clique indices
	std::vector<int> cliques;
	for (int cl = 0; cl < table.nCliques(); cl++) {
		if (table.cliqueIsEqual(cl))  cliques.push_back(cl);
	}

	// now process them one by one in the order they appear in the cliquetable
	std::vector<bool> covered(n, false);
	for (int cl: cliques) {
		// check if it has non empty intersection with already covered vars
		bool intersect = false;
		for (int lit: table.getClique(cl)) {
			const auto [var,isPos] = varFromLit(lit, n);
			if (covered[var]) {
				intersect = true;
				break;
			}
		}

		// if so, we need to skip it
		if (intersect)  continue;

		// otherwise we keep it and cover variables with it
		cc.add(table.getClique(cl), true);
		int lastAdded = cc.nCliques()-1;
		assert(0 <= lastAdded && lastAdded < cc.nCliques());
		for (int lit: cc.getClique(lastAdded)) {
			const auto [var,isPos] = varFromLit(lit, n);
			assert(!covered[var]);
			covered[var] = true;
		}
	}

	if (eqOnly) {
		cc.constructVarRepr();
		return cc;
	}

	/*************************************
	 * Then the rest
	 ************************************/

	// Count how much each clique covers the yet to be covered variables (either as is or negated)
	std::vector<int> cliqueCount(table.nCliques(), 0);
	for (int k = 0; k < (int)vars.size(); k++) {
		int var = vars[k];

		if (covered[var])  continue;

		for (int neg = 0; neg < 2; neg++) {
			int lit = makeLit(var, neg, n);
			for (int cl: table.getLit(lit))  cliqueCount[cl]++;
		}
	}

	// greedy step: assign each column to the largest clique covering it
	std::vector<int> var2clique(n, -1);
	for (int k = 0; k < (int)vars.size(); k++) {
		int var = vars[k];

		if (covered[var])  continue;

		int assignedClique = -1;
		int assignedCount = 0;

		for (int neg = 0; neg < 2; neg++) {
			int lit = makeLit(var, neg, n);
			for (int cl: table.getLit(lit)) {
				if (cliqueCount[cl] > assignedCount) {
					assignedClique = cl;
					assignedCount = cliqueCount[cl];
				}
			}
		}

		var2clique[var] = assignedClique;
	}

	// Recount how each clique (among the assigned ones only!) covers the variables
	std::fill(cliqueCount.begin(), cliqueCount.end(), 0);
	for (int k = 0; k < (int)vars.size(); k++) {
		int var = vars[k];

		if (var2clique[var] == -1)  continue;

		cliqueCount[var2clique[var]]++;
	}

	// Local adjustments to increase size of largest cover and possibly use fewer cliques
	for (int k = 0; k < (int)vars.size(); k++) {
		int var = vars[k];

		if (var2clique[var] == -1)  continue;

		// current assignment
		int assignedClique = var2clique[var];
		int assignedCount = cliqueCount[assignedClique];

		for (int neg = 0; neg < 2; neg++) {
			int lit = makeLit(var, neg, n);
			for (int cl: table.getLit(lit)) {
				if (cliqueCount[cl] >= assignedCount) {
					// move variable var from clique assignedClique to cl
					cliqueCount[cl]++;
					cliqueCount[assignedClique]--;
					assignedClique = cl;
					assignedCount = cliqueCount[cl];
				}
			}
		}

		var2clique[var] = assignedClique;
	}

	// Now collect all cliques that are actually in use in the cover (ignoring those that cover only one var)
	std::vector<bool> cliqueInUse(table.nCliques(), false);
	cliques.clear();
	for (int k = 0; k < (int)vars.size(); k++) {
		int var = vars[k];

		if (var2clique[var] == -1)  continue;

		int cl = var2clique[var];

		if (cliqueCount[cl] <= 1) {
			var2clique[var] = -1;
			cliqueCount[cl]--;
			continue;
		}

		if (!cliqueInUse[cl]) {
			cliques.push_back(cl);
			cliqueInUse[cl] = true;
		}
	}

	// Sort them by decreasing size
	std::stable_sort(cliques.begin(), cliques.end(), [&](int cl1, int cl2) {
		return (cliqueCount[cl1] > cliqueCount[cl2]);
	});

	// Finally add them to the cover
	std::vector<int> projectedClique;
	for (int cl: cliques) {
		projectedClique.clear();
		for (int lit: table.getClique(cl)) {
			const auto [var,isPos] = varFromLit(lit, n);
			if (var2clique[var] == cl)  projectedClique.push_back(lit);
		}
		assert(projectedClique.size() <= table.getClique(cl).size());
		assert(projectedClique.size() == cliqueCount[cl]);
		cc.add(projectedClique, false);
	}

#ifndef NDEBUG
	// Debugging checks
	std::fill(covered.begin(), covered.end(), false);
	int nCovered = 0;
	for (int cl = 0; cl < cc.nCliques(); cl++) {
		for (int lit: cc.getClique(cl)) {
			const auto [var,isPos] = varFromLit(lit, n);
			assert(!covered[var]);
			covered[var] = true;
			nCovered++;
		}
	}
	assert(nCovered == cc.nCovered());
	for (int cl = 1; cl < cc.nCliques(); cl++) {
		int prevEq = cc.cliqueIsEqual(cl-1);
		int thisEq = cc.cliqueIsEqual(cl);
		assert(prevEq >= thisEq);
	}
#endif

	cc.constructVarRepr();
	return cc;
}
