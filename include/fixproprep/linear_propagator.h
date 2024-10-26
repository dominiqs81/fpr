/**
 * @file   linear_propagator.h
 * @author Domenico Salvagnin
 */

#ifndef LINEAR_PROPAGATOR_H
#define LINEAR_PROPAGATOR_H

#include <queue>
#include "propagation.h"
#include "mip.h"


/* Propagates generic linear contraints */
class LinearPropagator : public PropagatorI
{
public:
	LinearPropagator(const MIPData& data);
	std::string name() const override { return "LinearPropagator"; }
	PropagatorPtr clone() const override;
	void init(const PropagationEngine& engine) override;
	void update(const PropagationEngine& engine, PropagationEngine::iterator mark) override;
	void undo(const PropagationEngine& engine, PropagationEngine::iterator mark) override;
	void propagate(PropagationEngine& engine, PropagationEngine::iterator mark, bool initialProp) override;
	void commit(const PropagationEngine& engine) override;
private:
	// Matrix data
	const MIPData& data;
	// State (note: activities are maintained by the engine itself!)
	PropagationEngine::iterator lastUpdated;
	struct State
	{
		double diameter = 0.0;
		bool infeas = false;
	};
	std::vector<State> states;
	std::vector<int> firstNonBin;
	// Helpers
	State computeState(const Domain& domain, dominiqs::SparseMatrix::view_type row) const;
	void propagateOneRowLessThan(PropagationEngine& engine, int i, double mult, double bound);
	void propagateOneRow(PropagationEngine& engine, int i);
};

#endif /* LINEAR_PROPAGATOR_H */
