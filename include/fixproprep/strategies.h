/**
 * @file strategies.h
 * @brief
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2022 Domenico Salvagnin
 */

#ifndef STRATEGIES_H
#define STRATEGIES_H

#include "mip.h"
#include "dfs.h"
#include <random>
#include <memory>


/* Policy to sort variables
 *
 * Optionally, a ranker can also choose the preferred value of variables!
 * */
class Ranker
{
public:
	Ranker(uint64_t _seed = 0) : seed{_seed} {}
	virtual ~Ranker() {}
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) = 0;
	virtual bool isChooser() const { return false; }
	virtual double choose(const MIPData& data, const PropagationEngine& engine, int var) { throw std::runtime_error("Not supported"); }
protected:
	uint64_t seed;
};

using RankerPtr = std::shared_ptr<Ranker>;


/* Policy to choose a value for a given variable */
class ValueChooser
{
public:
	ValueChooser(uint64_t _seed = 0) : seed{_seed} {}
	virtual ~ValueChooser() {}
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) = 0;
protected:
	uint64_t seed;
};

using ValuePtr = std::shared_ptr<ValueChooser>;


/* Policy to rank and pick preferred value for variables in a clique.
 * A separate function is provided for the single variable case.
 */
class CliqueStrategy
{
public:
	CliqueStrategy(uint64_t _seed = 0) : seed{_seed} {}
	virtual void operator()(const MIPData& data, const PropagationEngine& engine,
			std::span<const int> clique, bool isEq,
			std::span<int> vars, std::span<double> values) = 0;
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) = 0;
protected:
	uint64_t seed;
};

using CliqueStrategyPtr = std::shared_ptr<CliqueStrategy>;


/******** Clique Strategies ********/
class LeftRightC : public CliqueStrategy
{
public:
	virtual void operator()(const MIPData& data, const PropagationEngine& engine,
			std::span<const int> clique, bool isEq,
			std::span<int> vars, std::span<double> values) override
	{
		assert( clique.size() == vars.size() );
		assert( (int)values.size() == data.mip.ncols );
		const Domain& domain = engine.getDomain();
		size_t count = 0;
		for (int lit: clique) {
			const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
			vars[count++] = j;
			values[j] = isPos ? domain.ub(j) : domain.lb(j);
		}
		assert( clique.size() == count );
	}
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		return engine.getDomain().ub(var);
	}
};


class ByObjC : public CliqueStrategy
{
public:
	ByObjC(double _mult = 1.0) : CliqueStrategy{0}, mult{_mult} {}
	virtual void operator()(const MIPData& data, const PropagationEngine& engine,
			std::span<const int> clique, bool isEq,
			std::span<int> vars, std::span<double> values) override
	{
		assert( mult == 1.0 || mult == -1.0 );
		assert( clique.size() == vars.size() );
		assert( (int)values.size() == data.mip.ncols );
		const Domain& domain = engine.getDomain();
		// delta of a variable in the clique w.r.t. the value achieving minimum clique activity level
		// e.g. for a clique x1 + (1-x2) + x3 <= 1, the values achieving min clique activity are
		// (0,1,0) and the deltas are (+1,-1,+1)
		delta.resize(data.mip.ncols);
		size_t count = 0;
		for (int lit: clique) {
			const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
			vars[count++] = j;
			values[j] = isPos ? domain.ub(j) : domain.lb(j);
			delta[j] = isPos ? +1.0 : -1.0;
		}
		// sort variable by their objective delta w.r.t. base objective level
		// i.e., the objective when all variables minimize the clique activity (see above)
		std::stable_sort(vars.begin(), vars.end(), [&](int v1, int v2) {
			double score1 = mult * data.mip.objSense * delta[v1] * data.mip.obj[v1];
			double score2 = mult * data.mip.objSense * delta[v2] * data.mip.obj[v2];
			return (score1 < score2);
		});
		assert( clique.size() == count );
	}
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		assert( mult == 1.0 || mult == -1.0 );
		const Domain& domain = engine.getDomain();
		return (mult * data.mip.objSense * data.mip.obj[var] >= 0.0) ? domain.lb(var) : domain.ub(var);
	}
private:
	double mult;
	std::vector<double> delta;
};


inline CliqueStrategyPtr makeCliqueStrategy(const std::string& name)
{
	// select ranker
	if (name == "LR") {
		return CliqueStrategyPtr{new LeftRightC()};
	}
	else if (name == "badobj") {
		return CliqueStrategyPtr{new ByObjC(-1.0)};
	}
	else if (name == "goodobj") {
		return CliqueStrategyPtr{new ByObjC(+1.0)};
	}
	else {
		assert( false );
	}
	return CliqueStrategyPtr{};
}


/************* Rankers *************/
class LeftRight : public Ranker
{
public:
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		assert( data.mip.ncols == (data.nBinaries + data.nIntegers + data.nContinuous) );
		std::vector<int> sorted;
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] != 'C') sorted.push_back(j);
		}
		return sorted;
	}
};


class ByCliqueStrategy : public Ranker
{
public:
	ByCliqueStrategy(uint64_t seed, const std::string& csName) : Ranker{seed}
	{
		cs = makeCliqueStrategy(csName);
	}
	virtual bool isChooser() const override { return true; }
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		const Domain& domain = engine.getDomain();
		const CliqueCover& cc = data.cliquecover;
		std::mt19937_64 rndgen;
		rndgen.seed(seed);
		std::vector<int> sorted;
		std::vector<bool> added(data.mip.ncols, false);
		values.resize(data.mip.ncols);

		for (int cl = 0; cl < cc.nCliques(); cl++) {
			auto clique = cc.getClique(cl);
			bool isEq = cc.cliqueIsEqual(cl);
			std::vector<int> vars(clique.size());
			(*cs)(data, engine, clique, isEq, vars, values);

			for (int j: vars) {
				sorted.push_back(j);
				added[j] = true;
			}
		}

		// now add remaining non-continuous vars
		for (int j = 0; j < data.mip.ncols; j++) {
			if ((domain.type(j) != 'C') && !added[j]) {
				sorted.push_back(j);
				values[j] = (*cs)(data, engine, j);
				added[j] = true;
			}
		}

		return sorted;
	}
	virtual double choose(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		return values[var];
	}
protected:
	CliqueStrategyPtr cs;
	std::vector<double> values;
};


class ByType : public Ranker
{
public:
	ByType(bool _useCliqueCover = false) : useCliqueCover{_useCliqueCover} {}
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		assert( data.mip.ncols == (data.nBinaries + data.nIntegers + data.nContinuous) );
		const Domain& domain = engine.getDomain();

		// bucket sort by type (note that the sorting is stable)
		int startBin = 0;
		int startInt = startBin + data.nBinaries;
		std::vector<int> sorted(data.nBinaries + data.nIntegers);
		std::vector<bool> added(data.mip.ncols, false);
		if (useCliqueCover && data.cliquecover.nCliques()) {
			// first add variables as they appear in the clique cover
			for (int cl = 0; cl < data.cliquecover.nCliques(); cl++) {
				for (int lit: data.cliquecover.getClique(cl)) {
					const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
					if (dominiqs::equal(domain.lb(j), domain.ub(j)))  continue;
					assert( !added[j] );
					sorted[startBin++] = j;
					added[j] = true;
				}
			}
		}

		// then the rest
		for (int j = 0; j < data.mip.ncols; j++) {
			if (added[j])  continue;
			if (data.mip.xtype[j] == 'B')       sorted[startBin++] = j;
			else if (data.mip.xtype[j] == 'I')  sorted[startInt++] = j;
			added[j] = true;
		}
		assert( startBin == data.nBinaries );
		assert( startInt == (data.nBinaries  + data.nIntegers) );
		return sorted;
	}
private:
	bool useCliqueCover;
};


class ByIntegrality : public Ranker
{
public:
	ByIntegrality(const std::vector<double>& _xref) : xref(_xref) {}
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		assert( data.mip.ncols == (data.nBinaries + data.nIntegers + data.nContinuous) );
		assert( !xref.empty() );
		// always start with a bucket sort by type
		int startBin = 0;
		int startInt = startBin + data.nBinaries;
		std::vector<int> sorted(data.nBinaries + data.nIntegers);
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] == 'B')       sorted[startBin++] = j;
			else if (data.mip.xtype[j] == 'I')  sorted[startInt++] = j;
		}
		assert( startBin == data.nBinaries );
		assert( startInt == (data.nBinaries  + data.nIntegers) );

		// now we can assign a score and sort (stably!) within each bucket
		std::vector<double> score(data.mip.ncols);
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] == 'C')  continue;
			score[j] = dominiqs::integralityViolation(xref[j]);
		}
		std::stable_sort(sorted.begin(), sorted.begin() + data.nBinaries, [&](int v1, int v2) {
			return (score[v1] < score[v2]);
		});
		std::stable_sort(sorted.begin() + data.nBinaries, sorted.end(), [&](int v1, int v2) {
			return (score[v1] < score[v2]);
		});
		return sorted;
	}
protected:
	std::vector<double> xref;
};


class ByObj : public Ranker
{
public:
	ByObj(double _mult) : mult(_mult) {}
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		assert( mult == 1.0 || mult == -1.0 );
		assert( data.mip.ncols == (data.nBinaries + data.nIntegers + data.nContinuous) );
		// always start with a bucket sort by type
		int startBin = 0;
		int startInt = startBin + data.nBinaries;
		std::vector<int> sorted(data.nBinaries + data.nIntegers);
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] == 'B')       sorted[startBin++] = j;
			else if (data.mip.xtype[j] == 'I')  sorted[startInt++] = j;
		}
		assert( startBin == data.nBinaries );
		assert( startInt == (data.nBinaries  + data.nIntegers) );

		// now we can sort (stably!) within each bucket
		std::stable_sort(sorted.begin(), sorted.begin() + data.nBinaries, [&](int v1, int v2) {
			double score1 = mult * data.mip.objSense * data.mip.obj[v1];
			double score2 = mult * data.mip.objSense * data.mip.obj[v2];
			return (score1 < score2);
		});
		std::stable_sort(sorted.begin() + data.nBinaries, sorted.end(), [&](int v1, int v2) {
			double score1 = mult * data.mip.objSense * data.mip.obj[v1];
			double score2 = mult * data.mip.objSense * data.mip.obj[v2];
			return (score1 < score2);
		});
		return sorted;
	}
protected:
	double mult;
};


class RandomOrder : public Ranker
{
public:
	RandomOrder(uint64_t seed) : Ranker{seed} {}
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		// bucket sort by type
		std::vector<int> sorted(data.nBinaries + data.nIntegers);
		int startBin = 0;
		int startInt = startBin + data.nBinaries;
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] == 'B')       sorted[startBin++] = j;
			else if (data.mip.xtype[j] == 'I')  sorted[startInt++] = j;
		}
		assert( startBin == data.nBinaries );
		assert( startInt == (data.nBinaries  + data.nIntegers) );

		// random shuffle within each bucket
		std::mt19937_64 rndgen;
		rndgen.seed(seed);
		std::shuffle(sorted.begin(), sorted.begin() + data.nBinaries, rndgen);
		std::shuffle(sorted.begin() + data.nBinaries, sorted.end(), rndgen);
		return sorted;
	}
};


class Locks : public Ranker
{
public:
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		// always start with a bucket sort by type
		int startBin = 0;
		int startInt = startBin + data.nBinaries;
		std::vector<int> sorted(data.nBinaries + data.nIntegers);
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] == 'B')       sorted[startBin++] = j;
			else if (data.mip.xtype[j] == 'I')  sorted[startInt++] = j;
		}
		assert( startBin == data.nBinaries );
		assert( startInt == (data.nBinaries  + data.nIntegers) );

		// compute score based on locks
		std::vector<int> score(data.mip.ncols);
		int sentinel = 2*data.mip.nRows;
		for (int j = 0; j < data.mip.ncols; j++) {
			if (data.mip.xtype[j] == 'C')  continue;
			// both bounds are finite
			if (data.uplocks[j] > data.dnlocks[j])  score[j] = -data.dnlocks[j];
			else                                    score[j] = -data.uplocks[j];
		}

		// sort by increasing score within each bucket
		std::stable_sort(sorted.begin(), sorted.begin() + data.nBinaries, [&](int v1, int v2) {
			return (score[v1] < score[v2]);
		});
		std::stable_sort(sorted.begin() + data.nBinaries, sorted.end(), [&](int v1, int v2) {
			return (score[v1] < score[v2]);
		});

		return sorted;
	}
};


class ByCliqueCover : public Ranker
{
public:
	ByCliqueCover(uint64_t seed, bool useCorepoint) : Ranker{seed}, corepoint{useCorepoint} {}
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		const Domain& domain = engine.getDomain();
		std::mt19937_64 rndgen;
		rndgen.seed(seed);
		std::vector<int> sorted;
		std::vector<bool> added(data.mip.ncols, false);
		std::vector<double> weights(data.mip.ncols, false);

		for (int cl = 0; cl < data.cliquecover.nCliques(); cl++) {
			std::vector<int> vars;
			for (int lit: data.cliquecover.getClique(cl)) {
				const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
				if (dominiqs::equal(domain.lb(j), domain.ub(j)))  continue;
				vars.push_back(j);
				double value = corepoint ? data.zeroobj_corepoint[j] : data.zeroobj_xlp[j];
				if (isPos)  weights[j] = value;
				else        weights[j] = 1.0 - value;
			}

			// now "randomly sort" by weight
			// This is based on a O(n logn) algorithm for weighted random sampling without replacement
			// as described in:
			// https://blog.taboola.com/going-old-school-designing-algorithms-fast-weighted-sampling-production
			if (vars.size() > 1) {
				std::uniform_real_distribution<double> d;
				for (int j: vars) {
					double w = weights[j];
					weights[j] = std::log(d(rndgen)) / w;
				}
				std::stable_sort(vars.begin(), vars.end(), [&](int i, int j){
					return (weights[i] > weights[j]);
				});
			}

			// add them to sorted list
			sorted.insert(sorted.end(), vars.begin(), vars.end());
		}

		// now add remaining non-continuous vars
		for (int j: sorted)  added[j] = true;

		for (int j = 0; j < data.mip.ncols; j++) {
			if ((domain.type(j) != 'C') && !added[j])  sorted.push_back(j);
		}

		return sorted;
	}
protected:
	bool corepoint;
};


class ByCliques : public Ranker
{
public:
	virtual bool isChooser() const override { return true; }
	virtual std::vector<int> operator()(const MIPData& data, const PropagationEngine& engine) override
	{
		const Domain& domain = engine.getDomain();
		std::mt19937_64 rndgen;
		rndgen.seed(seed);
		std::vector<int> sorted;
		std::vector<bool> added(data.mip.ncols, false);
		values.resize(data.mip.ncols);
		const CliqueTable& cliquetable = data.cliquetable;
		int nCliques = cliquetable.nCliques();
		assert( nCliques > 0 );

		for (int cl = 0; cl < nCliques; cl++) {
			double sum = 0.0;
			double bestValue = 0.0;
			int bestVar = -1;
			bool bestPos;
			for (int lit: cliquetable.getClique(cl)) {
				const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
				if (isPos) {
					double xj = data.zeroobj_xlp[j];
					if (domain.lb(j) > 0.5) {
						// all variables in clique necessarily fixed
						bestVar = -1;
						break;
					}
					sum += xj;
					if ((xj > bestValue) && (domain.ub(j) > 0.5)) {
						bestValue = xj;
						bestVar = j;
						bestPos = true;
					}
				}
				else {
					double xj = 1.0 - data.zeroobj_xlp[j];
					if (domain.ub(j) < 0.5) {
						// all variables in clique necessarily fixed
						bestVar = -1;
						break;
					}
					sum += xj;
					if ((xj > bestValue) && (domain.lb(j) < 0.5)) {
						bestValue = xj;
						bestVar = j;
						bestPos = false;
					}
				}
			}

			if ((bestVar != -1) && sum >= (1.0 - 1e-5)) {
				// This clique is good for the rank: use it
				assert( domain.lb(bestVar) < domain.ub(bestVar) );
				if (!added[bestVar]) {
					sorted.push_back(bestVar);
					added[bestVar] = true;
					values[bestVar] = bestPos ? 1.0 : 0.0;
				}
				for (int lit: cliquetable.getClique(cl)) {
					const auto [j,isPos] = varFromLit(lit, data.mip.ncols);
					if (j == bestVar)  continue;
					if (domain.lb(j) == domain.ub(j))  continue;
					if (added[j])  continue;
					if (isPos) {
						sorted.push_back(j);
						added[j] = true;
						values[j] = 0.0;
					}
					else {
						sorted.push_back(j);
						added[j] = true;
						values[j] = 1.0;
					}
				}
			}
		}

		// now add remaining non-continuous vars
		// TODO: sort by increasing fractionality?
		for (int j = 0; j < data.mip.ncols; j++) {
			if ((domain.type(j) != 'C') && !added[j]) {
				sorted.push_back(j);
				// preferred value by nearest integer
				double frac = fractionalPart(data.zeroobj_xlp[j]);
				if (frac > 0.5)  values[j] = domain.ub(j);
				else             values[j] = domain.lb(j);
			}
		}

		return sorted;
	}
	virtual double choose(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		return values[var];
	}
protected:
	std::vector<double> values;
};


/************* Value Choosers *************/
class GoodObj : public ValueChooser
{
public:
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		const Domain& domain = engine.getDomain();
		return (data.mip.objSense*data.mip.obj[var] >= 0.0) ? domain.lb(var) : domain.ub(var);
	}
};


class BadObj : public ValueChooser
{
public:
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		const Domain& domain = engine.getDomain();
		return (data.mip.objSense*data.mip.obj[var] <= 0.0) ? domain.lb(var) : domain.ub(var);
	}
};


class Loose : public ValueChooser
{
public:
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		const Domain& domain = engine.getDomain();
		if (data.uplocks[var] > data.dnlocks[var])  return domain.lb(var);
		else if (data.uplocks[var] < data.dnlocks[var])  return domain.ub(var);
		// in case of ties consider objective, again in the loose direction
		return (data.mip.objSense*data.mip.obj[var] <  0.0) ? domain.lb(var) : domain.ub(var);
	}
};


class LooseDyn : public ValueChooser
{
public:
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		// compute dynamic locks based on current row activities
		int uplocks = 0;
		int dnlocks = 0;
		for (const auto& [i,coef]: data.mip.cols[var]) {
			if (engine.isRowRedundant(i))  continue;
			if (data.mip.sense[i] == 'E') {
				uplocks++;
				dnlocks++;
			}
			else {
				double mult = (data.mip.sense[i] == 'L') ? +1.0 : -1.0;
				if ((mult*coef) > 0.0)  uplocks++;
				else                    dnlocks++;
			}
		}

		// look at objective if variable has no locks
		if (!uplocks && !dnlocks) {
			double obj = data.mip.objSense * data.mip.obj[var];
			if (obj > 0.0)  uplocks++;
			else            dnlocks++;
		}

		const Domain& domain = engine.getDomain();
		if (uplocks >= dnlocks)  return domain.lb(var);
		return domain.ub(var);
	}
};


class RoundInt : public ValueChooser
{
public:
	RoundInt(const std::vector<double>& _xref) : xref(_xref) {}
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		const Domain& domain = engine.getDomain();
		assert( !xref.empty() );
		if (dominiqs::greaterEqualThan(xref[var], domain.ub(var), FEASTOL))  return domain.ub(var);
		if (dominiqs::lessEqualThan(xref[var], domain.lb(var), FEASTOL))  return domain.lb(var);
		double fp = dominiqs::fractionalPart(xref[var]);
		if (fp >= 0.5)  return dominiqs::ceilEps(xref[var]);
		else            return dominiqs::floorEps(xref[var]);
	}
protected:
	std::vector<double> xref;
};


class AlwaysUp : public ValueChooser
{
public:
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		return engine.getDomain().ub(var);
	}
};


class AlwaysDown : public ValueChooser
{
public:
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		return engine.getDomain().lb(var);
	}
};


class RandomValue : public ValueChooser
{
public:
	RandomValue(uint64_t seed) : ValueChooser{seed}, unif{0.0, 1.0}
	{
		rndgen.seed(seed);
	}
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		const Domain& domain = engine.getDomain();
		double lb = domain.lb(var);
		double ub = domain.ub(var);
		if (domain.type(var) == 'B') {
			if (unif(rndgen) < 0.5)  return lb;
			else                     return ub;
		}
		else {
			assert( domain.type(var) == 'I' );
			double value = dominiqs::floorEps(lb + unif(rndgen) * (ub - lb) + 0.5);
			value = std::min(value, ub);
			value = std::max(value, lb);
			return value;
		}
	}
protected:
	std::mt19937_64 rndgen;
	std::uniform_real_distribution<double> unif;
};


class RandomRelaxation : public ValueChooser
{
public:
	RandomRelaxation(uint64_t seed, const std::vector<double>& x)
		: ValueChooser{seed}, unif{0.0, 1.0}, xref{x}
	{
		rndgen.seed(seed);
	}
	virtual double operator()(const MIPData& data, const PropagationEngine& engine, int var) override
	{
		const Domain& domain = engine.getDomain();
		double lb = domain.lb(var);
		double ub = domain.ub(var);
		if (domain.type(var) == 'B') {
			if (unif(rndgen) >= xref[var])  return domain.lb(var);
			else                            return domain.ub(var);
		}
		else {
			assert( domain.type(var) == 'I' );
			double lb = dominiqs::floorEps(xref[var]);
			double ub = dominiqs::ceilEps(xref[var]);
			double fracPart = dominiqs::fractionalPart(xref[var]);
			double value = (unif(rndgen) >= fracPart) ? lb : ub;
			value = std::min(value, ub);
			value = std::max(value, lb);
			return value;
		}
	}
protected:
	std::mt19937_64 rndgen;
	std::uniform_real_distribution<double> unif;
	std::vector<double> xref;
};


class StaticStrategy : public DFSStrategy
{
public:
	StaticStrategy(const MIPData& _data) : data(_data) {}
	void setup(const PropagationEngine& engine, RankerPtr _ranker, ValuePtr _chooser)
	{
		ranker = _ranker;
		chooser = _chooser;
		assert( data.mip.ncols == (data.nBinaries + data.nIntegers + data.nContinuous) );
		sorted = (*ranker)(data, engine);
		assert( sorted.size() == (int)(data.mip.ncols - data.nContinuous) );
#if 0
		consoleLog("Variable order:");
		for (int j: sorted)  consoleLog("{}", data.mip.cNames[j]);
#endif
		// compute inverse map
		invmap.resize(data.mip.ncols);
		std::fill(invmap.begin(), invmap.end(), -1);
		for (int idx = 0; idx < (int)sorted.size(); idx++)  invmap[sorted[idx]] = idx;
	}
	std::vector<Branch> branch(const PropagationEngine& engine, bool nodeInfeas, const Branch& oldBranch) override
	{
		const Domain& domain = engine.getDomain();
		assert( oldBranch.type == BranchType::Variable );
		int startPos = (oldBranch.index >= 0) ? invmap[oldBranch.index] : 0;
		assert( startPos != -1 );

		// Branch according to order
		for (int idx = startPos; idx < (int)sorted.size(); idx++) {
			int var = sorted[idx];
			assert( data.mip.xtype[var] != 'C' );
			if (dominiqs::equal(domain.lb(var), domain.ub(var)))  continue;
			assert( (domain.ub(var) - domain.lb(var)) >= 1.0 );

			/* Choose a preferred value */
			double lb = domain.lb(var);
			double ub = domain.ub(var);
			double value = ranker->isChooser() ? ranker->choose(data, engine, var) : (*chooser)(data, engine, var);
			if (dominiqs::equal(value, lb)) {
				Branch preferred{BranchType::Variable, var, 'U', lb};
				Branch other{BranchType::Variable, var, 'L', lb + 1.0};
				return {preferred, other};
			}
			else {
				Branch preferred{BranchType::Variable, var, 'L', ub};
				Branch other{BranchType::Variable, var, 'U', ub - 1.0};
				return {preferred, other};
			}
		}
		return {};
	}
private:
	const MIPData& data;
	RankerPtr ranker;
	ValuePtr chooser;
	std::vector<int> sorted;
	std::vector<int> invmap;
};


inline RankerPtr makeRanker(const std::string& name, const Params& params, const MIPData& data)
{
	// select ranker
	if (name == "LR") {
		return RankerPtr{new LeftRight()};
	}
	else if (name == "type") {
		return RankerPtr{new ByType()};
	}
	else if (name == "typec") {
		return RankerPtr{new ByType(true)};
	}
	else if (name == "locks") {
		return RankerPtr{new Locks()};
	}
	else if (name == "random") {
		return RankerPtr{new RandomOrder(params.seed)};
	}
	else if (name == "ccover") {
		return RankerPtr{new ByCliqueStrategy(params.seed, "LR")};
	}
	else if (name == "ccover_badobj") {
		return RankerPtr{new ByCliqueStrategy(params.seed, "badobj")};
	}
	else if (name == "ccover_goodobj") {
		return RankerPtr{new ByCliqueStrategy(params.seed, "goodobj")};
	}
	else if (name == "cliques") {
		return RankerPtr{new ByCliqueCover(params.seed, true)};
	}
	else if (name == "cliques2") {
		return RankerPtr{new ByCliques()};
	}
	else if (name == "goodobj") {
		return RankerPtr{new ByObj(+1.0)};
	}
	else if (name == "badobj") {
		return RankerPtr{new ByObj(-1.0)};
	}
	else {
		throw std::runtime_error(fmt::format("Unknown ranker: {}", name));
	}
	return RankerPtr{};
}


inline ValuePtr makeValueChooser(const std::string& name, const Params& params, const MIPData& data)
{
	// select value chooser
	if (name == "good_obj") {
		return ValuePtr{new GoodObj()};
	}
	else if (name == "bad_obj") {
		return ValuePtr{new BadObj()};
	}
	else if (name == "random") {
		return ValuePtr{new RandomValue(params.seed)};
	}
	else if (name == "loose") {
		return ValuePtr{new Loose()};
	}
	else if (name == "loosedyn") {
		return ValuePtr{new LooseDyn()};
	}
	else if (name == "zeroobj_corepoint") {
		return ValuePtr{new RandomRelaxation(params.seed, data.zeroobj_corepoint)};
	}
	else if (name == "zeroobj_lp") {
		return ValuePtr{new RandomRelaxation(params.seed, data.zeroobj_xlp)};
	}
	else if (name == "corepoint") {
		return ValuePtr{new RandomRelaxation(params.seed, data.corepoint)};
	}
	else if (name == "lp") {
		return ValuePtr{new RandomRelaxation(params.seed, data.xlp)};
	}
	else if (name == "up") {
		return ValuePtr{new AlwaysUp()};
	}
	else if (name == "down") {
		return ValuePtr{new AlwaysDown()};
	}
	else {
		throw std::runtime_error(fmt::format("Unknown value chooser: {}", name));
	}
	return ValuePtr{};
}

#endif /* STRATEGIES_H */
