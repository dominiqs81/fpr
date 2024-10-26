/**
 * @file   Table-based propagators.h
 * @author Domenico Salvagnin
 */

#ifndef TABLE_PROPAGATORS_H
#define TABLE_PROPAGATORS_H

#include <queue>
#include "propagation.h"
#include "cliquetable.h"
#include "impltable.h"


/* Propagator for the cliques in a cliquetable.
 *
 * Note that equality cliques are propagated only as '<=' cliques.
 * The propagator is stateless.
 */
class CliquesPropagator : public PropagatorI
{
public:
	CliquesPropagator(const CliqueTable& ct) : PropagatorI{}, cliquetable(ct) {}
	std::string name() const override { return "CliquePropagator"; }
	PropagatorPtr clone() const override;
	void propagate(PropagationEngine& engine, PropagationEngine::iterator mark, bool initialProp) override;
private:
	const CliqueTable& cliquetable;
	void propagateOneLiteral(PropagationEngine& engine, int lit, bool& infeas);
};


/* Propagator for implication table 
 *
 * It propagates implications both forward and backward.
 */
class ImplPropagator : public PropagatorI
{
public:
	ImplPropagator(const ImplTable& impls) : PropagatorI{}, impltable(impls) {}
	std::string name() const override { return "ImplPropagator"; }
	PropagatorPtr clone() const override;
	void propagate(PropagationEngine& engine, PropagationEngine::iterator mark, bool initialProp) override;
private:
	const ImplTable& impltable;
	void forwardProp(PropagationEngine& engine, int binvar, bool value, bool& infeas);
	void backwardProp(PropagationEngine& engine, int impliedvar, bool isLower, double bound, bool& infeas);
};

#endif /* TABLE_PROPAGATORS_H */
