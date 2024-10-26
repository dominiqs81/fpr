/**
 * @file   propagation.h
 * @brief  Constraint Propagation API
 * @author Domenico Salvagnin
 */

#ifndef PROPAGATION_H
#define PROPAGATION_H

#include <vector>
#include <string>
#include <span>
#include <memory>
#include "fixproprep/mip.h"
#include "fixproprep/index_set.h"


/** @brief domain data structure.
 *
 * This keeps a copy of the domain (i.e., bounds and type) of variables.
 * It is on this copy that propagation operates.
 */

class Domain
{
public:
	/* Initialize a domain from given bounds and variable types */
	void init(std::span<const double> _lb, std::span<const double> _ub, std::span<const char> _xtype);
	/* Simple getters */
	int ncols() const { return (int)xtype.size(); }
	inline double lb(int var) const { return xlb[var]; }
	inline double ub(int var) const { return xub[var]; }
	inline char type(int var) const { return xtype[var]; }
	inline std::span<const double> lbs() const { return xlb; }
	inline std::span<const double> ubs() const { return xub; }
	/* Checks if new candidate bounds would be accepted given current bounds and parameters */
	bool isNewLowerBoundAcceptable(int var, double newBound) const;
	bool isNewUpperBoundAcceptable(int var, double newBound) const;
	/* Change a bound on a variable */
	void changeLowerBound(int var, double newBound) { xlb[var] = newBound; }
	void changeUpperBound(int var, double newBound) { xub[var] = newBound; }
	/* Shift a variable's domain by some amount */
	void shift(int var, double delta) { xlb[var] += delta; xub[var] += delta; }
	/* Domain (numerical parameters) */
	double feasTol = 1e-6;
	double zeroTol = 1e-9;
	double infinity = 1e8;
	double minContRed = 0.1;
	/* Variable names for debugging */
	std::vector<std::string> cNames;
private:
	std::vector<double> xlb;
	std::vector<double> xub;
	std::vector<char> xtype;
};


/* Forward declaration of PropagatorI */
class PropagatorI;
using PropagatorPtr = std::shared_ptr<PropagatorI>;


/** Bound change data structure */
struct BoundChange
{
public:
	enum class Type : uint8_t { SHIFT=0, UPPER, LOWER };
	BoundChange() = default;
	BoundChange(int j, Type t, double v, double old) : var(j), type(t), value(v), oldValue(old) {}
	Type type; /**< bound change type */
	int var; /**< column index */
	double value; /**< new value (delta for shifts) */
	double oldValue; /**< old value (unused for shifts) */
};


/* Propagation Engine API */
class PropagationEngine
{
public:
	PropagationEngine(const MIPData& _data);
	PropagationEngine(const PropagationEngine& other);
	/* Add a propagator to the engine */
	void add(PropagatorPtr prop) { propagators.push_back(prop); }
	/* Get a propagator by name */
	PropagatorPtr getPropagator(const std::string& name) const;
	/* Initialize domain and propagators */
	void init(std::span<const double> _lb, std::span<const double> _ub, std::span<const char> _xtype);
	/* Get a (const) reference to current domain */
	const Domain& getDomain() const { return domain; }
	/* Get a (const) reference to current MIP data */
	const MIPData& getMIPData() const { return data; }
	/* Domain changes */
	bool changeLowerBound(int var, double newBound);
	bool changeUpperBound(int var, double newBound);
	bool fix(int var, double value);
	void shift(int var, double delta);
	/* Actitivies */
	double getMinAct(int i) const { return minAct[i]; }
	double getMaxAct(int i) const { return maxAct[i]; }
	void recomputeRowActivity(int i);
	void refresh();
	/* Violation */
	const IndexSet<int>& violatedRows() const { return violated; }
	double violation() const { return totViol; }
	inline double rowViol(int i) const
	{
		return dominiqs::rowViol(minAct[i], maxAct[i], data.mip.sense[i], data.mip.rlb[i], data.mip.rub[i]);
	}
	inline bool isRowRedundant(int i) const
	{
		char sense = data.mip.sense[i];
		double rlb = data.mip.rlb[i];
		double rub = data.mip.rub[i];
		assert(rlb <= rub);
		switch(sense) {
			case 'L':
				return dominiqs::lessEqualThan(maxAct[i], rub, FEASTOL);
			case 'G':
				return dominiqs::greaterEqualThan(minAct[i], rlb, FEASTOL);
			case 'E':
				assert(rlb == rub);
				return (dominiqs::equal(minAct[i], maxAct[i], FEASTOL) &&
					dominiqs::equal(minAct[i], rlb, FEASTOL));
			case 'R':
				return (dominiqs::lessEqualThan(maxAct[i], rub, FEASTOL) &&
					dominiqs::greaterEqualThan(minAct[i], rlb, FEASTOL));
		}
		return false;
	}
	/* Propagate until fixpoint or infeasibility */
	bool propagate(bool initialProp);
	/* Lightweight propagation: only direct implications of the current changes */
	bool directImplications();
	/* API for keeping track of changes to the domain and backtrack */
	/* We cannot use a regular std::vector::iterator here, as that can be invalidated
	 * if memory is reallocated on insert. So we use the actual size instead.
	 * The price to pay is that we always need the domain reference as well,
	 * hence it is not an iterator in the STL sense...
	 */
	using iterator = size_t;
	inline iterator changesBegin() const { return 0; }
	inline iterator changesEnd() const { return stack.size(); }
	const BoundChange& change(iterator itr) const { return stack[itr]; }
	inline iterator mark() const { return changesEnd(); }
	void undo(iterator mark);
	const BoundChange& lastChange() const { return stack.back(); }
	/* Commit current set of changes */
	void commit();
	/* Enable/disable propagators */
	void enableAll();
	void disableAll();
	/* Engine parameters */
	int maxPasses = 100;
	void debugChecks() const;
	/* Effort counter */
	uint64_t work = 0;
protected:
	const MIPData& data;
	Domain domain;
	std::vector<PropagatorPtr> propagators;
	std::vector<BoundChange> stack;
	// state
	std::vector<double> minAct;
	std::vector<double> maxAct;
	IndexSet<int> violated;
	double totViol = 0.0;
	// helpers
	void recomputeViolation();
	void recomputeTotViol();
	void computeActivity(int i, double& minAct, double& maxAct) const;
private:
	void debugCheckRow(int i) const;
};


/* Propagator base class */
class PropagatorI
{
public:
	virtual ~PropagatorI() {}
	virtual std::string name() const = 0;
	virtual std::shared_ptr<PropagatorI> clone() const = 0;
	virtual void init(const PropagationEngine& engine) {}
	virtual void update(const PropagationEngine& engine, PropagationEngine::iterator mark) {}
	virtual void undo(const PropagationEngine& engine, PropagationEngine::iterator mark) {}
	virtual void propagate(PropagationEngine& engine, PropagationEngine::iterator mark, bool initialProp) = 0;
	virtual void commit(const PropagationEngine& engine) {}
	bool enabled() const { return _enabled; }
	void enable() { _enabled = true; }
	void disable() { _enabled = false; }
public:
	PropagationEngine::iterator lastPropagated;
	bool _enabled = true;
};



/* Debugging aids */
void printChangesSinceMark(const PropagationEngine& engine, PropagationEngine::iterator mark);

#endif /* PROPAGATION_H */
