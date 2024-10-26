/**
 * @file   propagation.c
 * @brief  Constraint Propagation API
 * @author Domenico Salvagnin
 */

//#define DEBUG_LOG
#ifdef DEBUG_LOG
static int DEBUG_LEVEL = 3;
#endif //< DEBUG_LOG

#include <cassert>
#include "fixproprep/propagation.h"
#include <utils/floats.h>
#include <utils/consolelog.h>

using namespace dominiqs;


void Domain::init(std::span<const double> _lb, std::span<const double> _ub, std::span<const char> _xtype)
{
	assert(_lb.size() == _ub.size());
	assert(_ub.size() == _xtype.size());
	xlb.resize(_lb.size());
	xub.resize(_ub.size());
	xtype.resize(_xtype.size());
	std::copy(_lb.begin(), _lb.end(), xlb.begin());
	std::copy(_ub.begin(), _ub.end(), xub.begin());
	std::copy(_xtype.begin(), _xtype.end(), xtype.begin());
	for (int j = 0; j < ncols(); j++) {
		/* Make sure that all binaries have proper bounds */
		assert(lb(j) <= ub(j));
		if (type(j) == 'B') {
			assert((lb(j) == 0.0) || (lb(j) == 1.0));
			assert((ub(j) == 0.0) || (ub(j) == 1.0));
		}
		/* Heuristically clip infinite bounds to some (not too large) value */
		if (xlb[j] <= -infinity)
			xlb[j] = -infinity / 1000.0;
		if (xub[j] >= +infinity)
			xub[j] = +infinity / 1000.0;
	}
}


bool Domain::isNewLowerBoundAcceptable(int var, double newBound) const
{
	double currLB = lb(var);
	double currUB = ub(var);
	assert(currLB <= currUB);

	if (fabs(newBound) >= infinity)  return false;

	double delta = newBound - currLB;
	if (delta <= feasTol)  return false;

	if (type(var) == 'C') {
		double thr = minContRed*(currUB - currLB);
		if (delta <= thr)  return false;
	}

	return true;
}


bool Domain::isNewUpperBoundAcceptable(int var, double newBound) const
{
	double currLB = lb(var);
	double currUB = ub(var);
	assert(currLB <= currUB);

	if (fabs(newBound) >= infinity)  return false;

	double delta = currUB - newBound;
	if (delta <= feasTol)  return false;

	if (type(var) == 'C') {
		double thr = minContRed*(currUB - currLB);
		if (delta <= thr)  return false;
	}

	return true;
}


PropagationEngine::PropagationEngine(const MIPData& _data) : data(_data),
	minAct(data.mip.nRows),
	maxAct(data.mip.nRows),
	violated(data.mip.nRows) {}


PropagationEngine::PropagationEngine(const PropagationEngine& other) :
	work(0),
	data(other.data),
	domain(other.domain),
	stack(other.stack),
	minAct(other.minAct),
	maxAct(other.maxAct),
	violated(other.violated),
	totViol(other.totViol)
{
	for (PropagatorPtr prop: other.propagators) {
		propagators.push_back(prop->clone());
	}
}


PropagatorPtr PropagationEngine::getPropagator(const std::string& name) const
{
	for (PropagatorPtr prop: propagators) {
		if (prop->name() == name)  return prop;
	}
	return PropagatorPtr{};
}


static inline bool needsRecomputation(double oldValue, double newValue)
{
	return (fabs(newValue) < (fabs(oldValue) * 1e-3));
}


void PropagationEngine::init(std::span<const double> _lb, std::span<const double> _ub, std::span<const char> _xtype)
{
	/* Initialize domain */
	assert(_lb.size() == _ub.size());
	assert(_ub.size() == _xtype.size());
	domain.init(_lb, _ub, _xtype);
	domain.cNames = data.mip.cNames;
	stack.clear();

	/* Compute activities and violations from scratch */
	refresh();

	/* Initialize propagators */
	for (const auto& prop: propagators)  prop->init(*this);
}


bool PropagationEngine::changeLowerBound(int var, double newBound)
{
	/* collect the current bounds of the column */
	double oldLowerBound = domain.lb(var);
	double oldUpperBound = domain.ub(var);
	assert(oldLowerBound <= oldUpperBound + domain.zeroTol);
	assert(oldLowerBound <= newBound + domain.zeroTol);

	/* detect infeasibilities */
	if (greaterThan(newBound, oldUpperBound, domain.feasTol)) {
		consoleDebug(2, "{} >= {} inconsistent with ub = {}", data.mip.cNames[var], newBound, oldUpperBound);
		return true;
	}
	else if (greaterEqualThan(newBound, oldUpperBound, domain.zeroTol)) {
		/* if the bounds are nearly identical, make them the same */
		newBound = oldUpperBound;
	}

	/* nothing to do if nothing changed */
	if (equal(newBound, oldLowerBound, domain.zeroTol))  return false;

	consoleDebug(3, "{} >= {}", data.mip.cNames[var], newBound);

	/* push bound change to the stack */
	stack.emplace_back(var, BoundChange::Type::LOWER, newBound, oldLowerBound);

	/* actually modify bound */
	domain.changeLowerBound(var, newBound);

	/* update activities */
	double delta = newBound - oldLowerBound;

	double deltaViol = 0.0;
	for (const auto& [i,coef]: data.mip.cols[var]) {
		double oldViol = rowViol(i);
		bool recomp = false;
		if (coef > 0.0) {
			double oldAct = minAct[i];
			minAct[i] += (coef * delta);
			recomp |= needsRecomputation(oldAct, minAct[i]);
		}
		else {
			double oldAct = maxAct[i];
			maxAct[i] += (coef * delta);
			recomp |= needsRecomputation(oldAct, maxAct[i]);
		}
		if (recomp)  recomputeRowActivity(i);
		double newViol = rowViol(i);
		if (newViol > domain.feasTol)  violated.add(i);
		else                           violated.remove(i);
		deltaViol += (newViol - oldViol);
#ifdef DEBUG_EXPENSIVE
		debugCheckRow(i);
#endif
	}
	work += data.mip.cols[var].size();
	double oldTotViol = totViol;
	totViol += deltaViol;
	if (needsRecomputation(oldTotViol, totViol))  recomputeTotViol();

	return false;
}


bool PropagationEngine::changeUpperBound(int var, double newBound)
{
	/* collect the current bounds of the column */
	double oldLowerBound = domain.lb(var);
	double oldUpperBound = domain.ub(var);
	assert(oldLowerBound <= oldUpperBound + domain.zeroTol);
	assert(newBound <= oldUpperBound + domain.zeroTol);

	/* detect infeasibilities */
	if (lessThan(newBound, oldLowerBound, domain.feasTol)) {
		consoleDebug(2, "{} <= {} inconsistent with lb = {}", data.mip.cNames[var], newBound, oldLowerBound);
		return true;
	}
	else if (lessEqualThan(newBound, oldLowerBound, domain.zeroTol)) {
		/* if the bounds are nearly identical, make them the same */
		newBound = oldLowerBound;
	}

	/* nothing to do if nothing changed */
	if (equal(newBound, oldUpperBound, domain.zeroTol))  return false;

	consoleDebug(3, "{} <= {}", data.mip.cNames[var], newBound);

	/* push bound change to the stack */
	stack.emplace_back(var, BoundChange::Type::UPPER, newBound, oldUpperBound);

	/* actually modify bound */
	domain.changeUpperBound(var, newBound);

	/* update activities */
	double delta = newBound - oldUpperBound;

	double deltaViol = 0.0;
	for (const auto& [i,coef]: data.mip.cols[var]) {
		double oldViol = rowViol(i);
		bool recomp = false;
		if (coef < 0.0) {
			double oldAct = minAct[i];
			minAct[i] += (coef * delta);
			recomp |= needsRecomputation(oldAct, minAct[i]);
		}
		else {
			double oldAct = maxAct[i];
			maxAct[i] += (coef * delta);
			recomp |= needsRecomputation(oldAct, maxAct[i]);
		}
		if (recomp)  recomputeRowActivity(i);
		double newViol = rowViol(i);
		if (newViol > domain.feasTol)  violated.add(i);
		else                           violated.remove(i);
		deltaViol += (newViol - oldViol);
#ifdef DEBUG_EXPENSIVE
		debugCheckRow(i);
#endif
	}
	work += data.mip.cols[var].size();
	double oldTotViol = totViol;
	totViol += deltaViol;
	if (needsRecomputation(oldTotViol, totViol))  recomputeTotViol();

	return false;
}


bool PropagationEngine::fix(int var, double value)
{
	bool infeas = changeLowerBound(var, value);
	if (infeas) return infeas;
	infeas = changeUpperBound(var, value);
	return infeas;
}


void PropagationEngine::shift(int var, double delta)
{
	if (delta == 0.0)  return;

	/* update domain */
	double oldLowerBound = domain.lb(var);
	double oldUpperBound = domain.ub(var);

	if (oldLowerBound == oldUpperBound) {
		consoleDebug(3, "{} {} -> {}", data.mip.cNames[var], oldUpperBound, oldUpperBound + delta);
	}
	else {
		consoleDebug(3, "{} [{},{}] -> [{},{}]", data.mip.cNames[var], oldLowerBound, oldUpperBound,
			oldLowerBound + delta, oldUpperBound + delta);
	}

	/* push bound change to the stack */
	stack.emplace_back(var, BoundChange::Type::SHIFT, delta, 0.0);

	/* actually modify bounds */
	domain.shift(var, delta);

	double deltaViol = 0.0;
	for (const auto& [i,coef]: data.mip.cols[var]) {
		double oldViol = rowViol(i);
		bool recomp = false;
		double oldMinAct = minAct[i];
		minAct[i] += (coef * delta);
		recomp |= needsRecomputation(oldMinAct, minAct[i]);
		double oldMaxAct = maxAct[i];
		maxAct[i] += (coef * delta);
		recomp |= needsRecomputation(oldMaxAct, maxAct[i]);
		if (recomp)  recomputeRowActivity(i);
		double newViol = rowViol(i);
		if (newViol > domain.feasTol)  violated.add(i);
		else                           violated.remove(i);
		deltaViol += (newViol - oldViol);
#ifdef DEBUG_EXPENSIVE
		debugCheckRow(i);
#endif
	}
	work += data.mip.cols[var].size();
	double oldTotViol = totViol;
	totViol += deltaViol;
	if (needsRecomputation(oldTotViol, totViol))  recomputeTotViol();
}


void PropagationEngine::recomputeRowActivity(int i)
{
	computeActivity(i, minAct[i], maxAct[i]);
	work += data.mip.rows[i].size();
}


void PropagationEngine::refresh()
{
	/* Evaluate current mininimum and maximum activities for each row */
	int m = data.mip.nRows;
	for (int i = 0; i < m; i++) {
		recomputeRowActivity(i);
	}

	/* Compute current violation and set of violated constraints */
	recomputeViolation();
}


bool PropagationEngine::propagate(bool initialProp)
{
	/* Nothing to do if we are already infeasible */
	if (!violated.empty())  return true;

	int nPasses = 0;
	bool infeas = false;

	while (!infeas) {
		/* stop if we are taking too many passes */
		if (nPasses >= maxPasses)  break;

		iterator prevEnd = changesEnd();

		/* One pass through all propagators */
		for (const auto& prop: propagators) {
			/* Skip disable propagators */
			if (!prop->enabled())  continue;

			/* Run each propagator until it reaches its own fixpoint or infeasibility */
			bool initialFirstPass = initialProp;

			while (true) {
				iterator currEnd = changesEnd();
				assert(currEnd >= prop->lastPropagated);
				bool hasChanges = (prop->lastPropagated != currEnd);

				/* Fixpoint reached */
				if (!hasChanges && !initialFirstPass)  break;

				/* Propagate */
				prop->propagate(*this, currEnd, initialFirstPass);
				prop->lastPropagated = currEnd;
				initialFirstPass = false;

				/* Stop on infeasibility */
				if (!violated.empty()) {
					// consoleLog("{} infeasible", prop->name());
					infeas = true;
					break;
				}
			}

			/* No need to run the next propagators */
			if (infeas)  break;
		}

		/* Stop if there was no additional reduction in the last pass */
		if (changesEnd() == prevEnd)  break;

		/* Stop if infeasibility was detected */
		if (infeas)  break;

		nPasses++;
	}

	return infeas;
}


bool PropagationEngine::directImplications()
{
	iterator currEnd = changesEnd();

	/* One pass through all propagators */
	for (const auto& prop: propagators) {
		/* Skip disable propagators */
		if (!prop->enabled())  continue;

		/* Run each propagator once */
		assert(currEnd >= prop->lastPropagated);

		/* Propagate */
		prop->propagate(*this, currEnd, false);
		prop->lastPropagated = currEnd;

		/* Stop on infeasibility */
		if (!violated.empty()) {
			// consoleLog("{} infeasible", prop->name());
			return true;
		}
	}

	return false;
}


void PropagationEngine::undo(iterator mark)
{
	/* Backtrack domain one change at time, updating activites */
	iterator end = changesEnd();
	assert( end >= mark );
	size_t dist = end - mark;

	for (size_t i = 0; i < dist; i++) {
		const BoundChange& bdchg = change(end - (i+1));

		/* Undo change to domain, but without changing the stack yet */
		switch(bdchg.type) {
			case BoundChange::Type::SHIFT:
				domain.shift(bdchg.var, -bdchg.value);
				break;
			case BoundChange::Type::UPPER:
				domain.changeUpperBound(bdchg.var, bdchg.oldValue);
				break;
			case BoundChange::Type::LOWER:
				domain.changeLowerBound(bdchg.var, bdchg.oldValue);
				break;
			default:
				assert( false );
		}

		/* Undo changes to activities */
		double delta = bdchg.value - bdchg.oldValue;
		int var = bdchg.var;

		double deltaViol = 0.0;
		for (const auto& [i,coef]: data.mip.cols[var]) {
			double oldViol = rowViol(i);
			double oldMinAct = minAct[i];
			double oldMaxAct = maxAct[i];
			switch (bdchg.type) {
				case BoundChange::Type::SHIFT:
					minAct[i] -= (coef * delta);
					maxAct[i] -= (coef * delta);
					break;
				case BoundChange::Type::UPPER:
					if (coef > 0.0)  maxAct[i] -= (coef * delta);
					else             minAct[i] -= (coef * delta);
					break;
				case BoundChange::Type::LOWER:
					if (coef > 0.0)  minAct[i] -= (coef * delta);
					else             maxAct[i] -= (coef * delta);
					break;
				default:
					assert( false );
			}
			if ( needsRecomputation(oldMinAct, minAct[i]) ||
				 needsRecomputation(oldMaxAct, maxAct[i]) ) {
				recomputeRowActivity(i);
			}
			double newViol = rowViol(i);
			if (newViol > domain.feasTol)  violated.add(i);
			else                           violated.remove(i);
			deltaViol += (newViol - oldViol);
#ifdef DEBUG_EXPENSIVE
			debugCheckRow(i);
#endif
		}
		work += data.mip.cols[var].size();
		totViol += deltaViol;
	}

	/* Backtrack propagators' states */
	for (const auto& prop: propagators) {
		prop->undo(*this, mark);
		if (prop->lastPropagated > mark) {
			prop->lastPropagated = mark;
		}
	}

	/* Now we can resize the stack */
	stack.resize(mark);
	assert( changesEnd() == mark );
}


void PropagationEngine::commit()
{
	for (const auto& prop: propagators) {
		prop->commit(*this);
		prop->lastPropagated = changesBegin();
	}

	stack.clear();
}


void PropagationEngine::enableAll()
{
	for (PropagatorPtr prop: propagators)  prop->enable();
}


void PropagationEngine::disableAll()
{
	for (PropagatorPtr prop: propagators)  prop->disable();
}


/* Compute current violation and set of violated constraints */
void PropagationEngine::recomputeViolation()
{
	int m = data.mip.nRows;
	violated.clear();
	totViol = 0.0;
	for (int i = 0; i < m; i++) {
		double viol = rowViol(i);
		assert( viol >= 0.0 );

		totViol += viol;
		if (viol > domain.feasTol)  violated.add(i);
	}
}


/* Recompute totViol only from the current set of violated rows */
void PropagationEngine::recomputeTotViol()
{
	totViol = 0.0;
	for (int i: violated.data()) {
		double viol = rowViol(i);
		assert( viol >= 0.0 );
		totViol += viol;
	}
}


void PropagationEngine::computeActivity(int i, double& minA, double& maxA) const
{
	minA = 0.0;
	maxA = 0.0;
	for (const auto& [var,coef]: data.mip.rows[i]) {
		assert((var >= 0) && (var < data.mip.ncols));
		assert(coef != 0.0);
		double lb = domain.lb(var);
		assert( lb > -domain.infinity );
		double ub = domain.ub(var);
		assert( ub < +domain.infinity );

		if (coef > 0.0) {
			minA += (coef * lb);
			maxA += (coef * ub);
		}
		else {
			minA += (coef * ub);
			maxA += (coef * lb);
		}
	}
}


void PropagationEngine::debugCheckRow(int i) const
{
	double minA;
	double maxA;
	computeActivity(i, minA, maxA);
	assert( relEqual(minA, minAct[i], domain.zeroTol) );
	assert( relEqual(maxA, maxAct[i], domain.zeroTol) );
	double viol = rowViol(i);
	assert( (viol > domain.feasTol) == violated.has(i) );
}


void PropagationEngine::debugChecks() const
{
	int m = data.mip.nRows;
	for (int i = 0; i < m; i++)  debugCheckRow(i);
}


void printChangesSinceMark(const PropagationEngine& engine, PropagationEngine::iterator mark)
{
	PropagationEngine::iterator end = engine.changesEnd();
	const Domain& domain = engine.getDomain();
	while (mark != end) {
		const BoundChange& bdchg = engine.change(mark++);
		int var = bdchg.var;
		assert( (0 <= var) && (var < domain.ncols()) );
		switch (bdchg.type) {
			case (BoundChange::Type::SHIFT):
				consoleLog("[{},{}] {} shifted by {}", var, domain.type(var), domain.cNames[var], bdchg.value);
				break;
			case (BoundChange::Type::UPPER):
				consoleLog("[{},{}] {} <= {}", var, domain.type(var), domain.cNames[var], bdchg.value);
				break;
			case (BoundChange::Type::LOWER):
				consoleLog("[{},{}] {} >= {}", var, domain.type(var), domain.cNames[var], bdchg.value);
				break;
			default:
				assert( false );
				break;
		}
	}
}
