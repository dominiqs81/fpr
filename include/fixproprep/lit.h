/**
 * @file   lit.h
 * @brief  Simple utilities for literals
 *
 * @author Domenico Salvagnin dominiqs@gmail.com
 */

#ifndef LIT_HELPERS_H
#define LIT_HELPERS_H

/* Helpers to convert from literals to column indices and back */
inline int posLit(int var, int ncols) { return var; }
inline int negLit(int var, int ncols) { return var+ncols; }
inline int makeLit(int var, int ncols, bool isPos) { return isPos ? posLit(var, ncols) : negLit(var, ncols); }
inline std::pair<int,bool> varFromLit(int lit, int ncols)
{
	if (lit < ncols) return {lit, true};
	return {lit-ncols, false};
}
inline int negateLit(int lit, int ncols)
{
	return (lit < ncols) ? (lit+ncols) : (lit-ncols);
}

#endif /* LIT_HELPERS_H */
