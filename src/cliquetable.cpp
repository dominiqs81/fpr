/**
 * @file cliquetable.h
 * @brief Cliquetable
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#include "fixproprep/cliquetable.h"
#include <fmt/format.h>


void CliqueTable::add(std::span<const int> clique, bool isEqual)
{
	assert(!clique.empty());
	cliques.add(clique);
	type.push_back(isEqual ? 'E' : 'L');
	hasLitwise = false;
}


void CliqueTable::constructLitWiseRepr()
{
	if (hasLitwise)
		return;

	literals = cliques.transpose();
	hasLitwise = true;
}
