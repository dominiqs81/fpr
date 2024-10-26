/**
 * @file   linear_propagator.c
 * @author Domenico Salvagnin
 */

#include "fixproprep/linear_propagator.h"
#include <utils/consolelog.h>

using namespace dominiqs;


LinearPropagator::LinearPropagator(const MIPData& _data) : PropagatorI{}, data{_data} {}


PropagatorPtr LinearPropagator::clone() const
{
	return PropagatorPtr{new LinearPropagator{*this}};
}


/* Compute state for a given row */
LinearPropagator::State LinearPropagator::computeState(const Domain& domain, SparseMatrix::view_type row) const
{
	State s{}; //< all fields initialized to zero
	int n = domain.ncols();

	for (const auto& [var,coef]: row) {
		assert((var >= 0) && (var < n));
		assert(coef != 0.0);
		double lb = domain.lb(var);
		assert( lb > -domain.infinity );
		double ub = domain.ub(var);
		assert( ub < +domain.infinity );

		/* update diameter */
		double spread = ub - lb;
		if (domain.type(var) == 'C')  spread *= (1.0-domain.minContRed);
		else                          spread = std::max(spread - domain.feasTol, 0.0);
		s.diameter = std::max(s.diameter, fabs(coef) * spread);
	}

	return s;
}


void LinearPropagator::init(const PropagationEngine& engine)
{
	const Domain& domain = engine.getDomain();
	int m = data.mip.nRows;
	states.resize(m);
	for (int i = 0; i < m; i++)  states[i] = computeState(domain, data.mip.rows[i]);

	/* Rows are normalized: find position of first non binary variable.
	 * This is not part of computeState as this piece of information nevers gets invalidated.
	 */
	firstNonBin.resize(m);
	for (int i = 0; i < m; i++) {
		const auto& row = data.mip.rows[i];
		const int* idx = row.idx();
		int cnt = row.size();
		int k = 0;

		while ((k < cnt) && domain.type(idx[k]) == 'B')  k++;
		firstNonBin[i] = k;
	}

	lastUpdated = engine.changesEnd();
}


void LinearPropagator::update(const PropagationEngine& engine, PropagationEngine::iterator mark)
{
	assert(lastUpdated <= mark);

	while (lastUpdated < mark) {
		const BoundChange& bdchg = engine.change(lastUpdated++);

		/* Shifts reset the infeasibility status! */
		if (bdchg.type != BoundChange::Type::SHIFT)  continue;

		for (const auto& [i,coef]: data.mip.cols[bdchg.var])  states[i].infeas = false;
	}
}


void LinearPropagator::undo(const PropagationEngine& engine, PropagationEngine::iterator mark)
{
	/* Nothing to do if we are backtracking to a position that is still
	 * in the future w.r.t. this propagator.
	 */
	if (lastUpdated <= mark)  return;

	/* update activities incrementally */
	const Domain& domain = engine.getDomain();
	while (lastUpdated != mark) {
		const BoundChange& bdchg = engine.change(--lastUpdated);
		int j = bdchg.var;

		for (const auto& [i,coef]: data.mip.cols[j])
		{
			State& s = states[i];

			/* Reset infeas status on backtrack */
			s.infeas = false;

			/* update diameter */
			if (bdchg.type != BoundChange::Type::SHIFT) {
				double spread = domain.ub(j) - domain.lb(j);
				if (domain.type(j) == 'C')  spread *= (1.0 - domain.minContRed);
				else                        spread = std::max(spread - domain.feasTol, 0.0);
				s.diameter = std::max(s.diameter, fabs(coef) * spread);
			}
		}
	}
}


/* Propagate row i as a <= constraint */
void LinearPropagator::propagateOneRowLessThan(PropagationEngine& engine, int i, double mult, double bound)
{
	State& s = states[i];

	assert((mult == 1.0) || (mult == -1.0));

	const Domain& domain = engine.getDomain();

	/* Perform simple checks first */
	double act = (mult > 0) ? engine.getMinAct(i) : engine.getMaxAct(i);
	double slack = mult*(bound - act);

	if (slack < -domain.feasTol) {
		// consoleLog("Row {} mult={} infeasible", data.mip.rNames[i], mult);
		s.infeas = true;
		return;
	}

	/* Make sure slack is non-negative */
	if (slack < 0.0)  slack = 0.0;

	if (lessEqualThan(s.diameter, slack, domain.feasTol)) {
		// consoleLog("Row {} mult={} cannot propagate: {} <= {}", data.mip.rNames[i], mult, s.diameter, slack);
		// consoleLog("minAct={} maxAct={}", s.minAct, s.maxAct);
		return;
	}

	// consoleLog("Tighten vars from row {} mult={}", data.mip.rNames[i], mult);

	/* Loop over the row, tighten variables and update diameter */
	bool infeas = false;
	const auto& row = data.mip.rows[i];
	const int* idx = row.idx();
	const double* coefs = row.coef();
	int size = row.size();
	double newDiameter = 0.0;

	/* First loop over binaries */
	for (int k = 0; k < firstNonBin[i]; k++) {
		int j = idx[k];
		double coef = mult*coefs[k];

		/* Stop if abs(coef) is below slack (rows are normalized!) */
		if (lessEqualThan(fabs(coef), slack, domain.feasTol)) {
			/* We take the diameter of this variable as an upper bound
			 * on the diameters of remaining binaries we skip.
			 */
			newDiameter = std::max(newDiameter, fabs(coef) * (1.0 - domain.feasTol)); 
			break;
		}

		/* Skip fixed binaries */
		if (domain.lb(j) == domain.ub(j))  continue;

		assert(fabs(coef) > slack);
		if (coef > 0.0) {
			/* Fix binary to zero */
			infeas = engine.changeUpperBound(j, 0.0);
			assert(!infeas);
		}
		else {
			/* Fix binary to one */
			infeas = engine.changeLowerBound(j, 1.0);
			assert(!infeas);
		}
	}
	/* Then the rest */
	for (int k = firstNonBin[i]; k < size; k++) {
		int j = idx[k];
		double coef = mult*coefs[k];

		/* Compute current diameter of this variable */
		double thisDiameter = domain.ub(j) - domain.lb(j);
		if (domain.type(j) == 'C')  thisDiameter *= (1.0-domain.minContRed);
		else                        thisDiameter = std::max(thisDiameter - domain.feasTol, 0.0);
		thisDiameter *= fabs(coef);

		/* Do not derive bounds from tiny coefficients */
		if (fabs(coef) <= domain.feasTol) {
			newDiameter = std::max(newDiameter, thisDiameter);
			continue;
		}

		/* Skip fixed vars */
		if ((domain.ub(j) - domain.lb(j)) <= domain.feasTol) {
			newDiameter = std::max(newDiameter, thisDiameter);
			continue;
		}

		bool tightened = false;
		if (coef > 0.0) {
			double delta = slack / coef;
			assert(delta >= 0.0);
			double newBound = domain.lb(j) + delta;
			if (domain.type(j) == 'I')  newBound = floorEps(newBound, domain.feasTol);
			if (domain.isNewUpperBoundAcceptable(j, newBound)) {
				infeas = engine.changeUpperBound(j, newBound);
				assert(!infeas);
				tightened = true;
			}
		}
		else if (coef < 0.0) {
			double delta = slack / coef;
			assert(delta <= 0.0);
			double newBound = domain.ub(j) + delta;
			if (domain.type(j) == 'I')  newBound = ceilEps(newBound, domain.feasTol);
			if (domain.isNewLowerBoundAcceptable(j, newBound)) {
				infeas = engine.changeLowerBound(j, newBound);
				assert(!infeas);
				tightened = true;
			}
		}

		if (tightened) {
			/* Recompute diameter of this variable as it was changed */
			thisDiameter = domain.ub(j) - domain.lb(j);
			if (domain.type(j) == 'C')  thisDiameter *= (1.0-domain.minContRed);
			else                        thisDiameter = std::max(thisDiameter - domain.feasTol, 0.0);
			thisDiameter *= fabs(coef);
		}
		newDiameter = std::max(newDiameter, thisDiameter);
	}
	engine.work += size;
	s.diameter = newDiameter;

	update(engine, engine.changesEnd());
}


/* propagate a single linear constraint */
void LinearPropagator::propagateOneRow(PropagationEngine& engine, int i)
{
	char sense = data.mip.sense[i];

	bool rowHasLB = (sense != 'L');
	bool rowHasUB = (sense != 'G');
	double rowLB = data.mip.rlb[i];
	double rowUB = data.mip.rub[i];

	const Domain& domain = engine.getDomain();
	State& s = states[i];
	if (s.infeas)  return;

	/* Is the whole row fixed? */
	double minAct = engine.getMinAct(i);
	double maxAct = engine.getMaxAct(i);
	if ((fabs(maxAct - minAct) <= domain.feasTol)) {
		/* If all variables in the row are fixed, then the diameter must be zero */
		states[i].diameter = 0.0;
		/* Check if the row became infeasible */
		if (rowHasLB && lessThan(maxAct, rowLB, domain.feasTol)) {
			// consoleLog("Row {} mult={} infeasible", data.mip.rNames[i], 1.0);
			states[i].infeas = true;
		}
		else if (rowHasUB && greaterThan(minAct, rowUB, domain.feasTol)) {
			// consoleLog("Row {} mult={} infeasible", data.mip.rNames[i], -1.0);
			states[i].infeas = true;
		}
		/* If not, then all variables must be fixed and there is nothing to do anyway */
		return;
	}

	if (rowHasUB)  propagateOneRowLessThan(engine, i, 1.0, rowUB);
	if (rowHasLB)  propagateOneRowLessThan(engine, i, -1.0, rowLB);
}


void LinearPropagator::propagate(PropagationEngine& engine, PropagationEngine::iterator mark, bool initialProp)
{
	/* Update state */
	update(engine, engine.changesEnd());

	std::queue<int> qrows;
	if (initialProp) {
		/* Enqueue all rows */
		for (int i = 0; i < data.mip.nRows; i++)  qrows.push(i);
	}
	else {
		/* Loop over list of bound changes */
		assert(lastPropagated <= mark);
		PropagationEngine::iterator itr = lastPropagated;
		while (itr != mark) {
			const auto& bdchg = engine.change(itr++);
			int j = bdchg.var;
			for (const auto& [i,coef]: data.mip.cols[j]) {
				const State& s = states[i];

				if (s.infeas)  continue;

				qrows.push(i);
			}
		}
	}

	while (!qrows.empty()) {
		int i = qrows.front(); qrows.pop();
		propagateOneRow(engine, i);
		if (states[i].infeas) {
			// consoleLog("Linear propagator infeasiblity on row {}", i);
			break;
		}
	}
}


void LinearPropagator::commit(const PropagationEngine& engine)
{
	PropagationEngine::iterator end = engine.changesEnd();
	update(engine, end);
	lastUpdated = engine.changesBegin();
}
