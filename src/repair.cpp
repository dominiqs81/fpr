//#define DEBUG_LOG
#ifdef DEBUG_LOG
static int DEBUG_LEVEL = 1;
#endif //< DEBUG_LOG


#include "fixproprep/repair.h"
#include "fixproprep/index_queue.h"
#include <utils/floats.h>
#include <utils/consolelog.h>
#include <utils/timer.h>
#include <utils/app.h>

using namespace dominiqs;


RepairMIP::RepairMIP(const MIPData& _data, const Params& _params)
	: data(_data), params(_params)
{
	// initialize random engine
	rndgen.seed(params.seed);
}


static void applyShift(PropagationEngine& solution, int var, double delta)
{
	const Domain& domain = solution.getDomain();
	const MIPInstance& mip = solution.getMIPData().mip;

	double oldLB = domain.lb(var);
	double oldUB = domain.ub(var);
	solution.shift(var, delta);
	if (oldLB == oldUB) {
		consoleDebug(3, "Apply shift: {} [{}] {} -> {}: viol = {}", mip.cNames[var], domain.type(var),
			oldLB, domain.lb(var), solution.violation());
	}
	else {
		consoleDebug(3, "Apply shift: {} [{}] [{},{}] -> [{},{}]: viol = {}", mip.cNames[var], domain.type(var),
			oldLB, oldUB, domain.lb(var), domain.ub(var), solution.violation());
	}
}


static void evalShift(const PropagationEngine& solution, int var, double delta, int pickedRow, bool& isCand, double& damage, uint64_t& work)
{
	const Domain& domain = solution.getDomain();
	const MIPInstance& mip = solution.getMIPData().mip;

	/* loop over affected rows and evaluate damage */
	isCand = true;
	damage = 0.0;
	for (const auto [i,a]: mip.cols[var]) {
		double minAct = solution.getMinAct(i);
		double maxAct = solution.getMaxAct(i);
		double deltaAct = delta*a;
		double oldViol = rowViol(minAct, maxAct, mip.sense[i], mip.rlb[i], mip.rub[i]);
		double newViol = rowViol(minAct+deltaAct, maxAct+deltaAct, mip.sense[i], mip.rlb[i], mip.rub[i]);
		if (newViol > oldViol)  damage += (newViol - oldViol);

		if ((i == pickedRow) && (newViol > oldViol)) {
			consoleDebug(4, "Not improving on picked {}: deltaAct = {} oldViol = {} newViol = {}",
				pickedRow, deltaAct, oldViol, newViol);
			isCand = false;
			break;
		}
	}
	work += mip.cols[var].size();
	consoleDebug(4, "Eval shift of {} for {} [cnt={}]: iscand = {} damage = {}",
			delta, mip.cNames[var], mip.cols[var].size(), isCand, damage);
}


void RepairMIP::walk(PropagationEngine& solution)
{
	int n = data.mip.ncols;
	int m = data.mip.nRows;
	const auto& sense = data.mip.sense;
	const Domain& domain = solution.getDomain();

	// init best
	const IndexSet<int>& violated = solution.violatedRows();
	double bestViol = solution.violation();
	int bestNviol = violated.size();
	PropagationEngine::iterator initMark = solution.mark();
	PropagationEngine::iterator bestMark = solution.mark();
	int nZeroDamage = 0;
	int nMinDamage = 0;
	int nRandom = 0;

	// main walk loop
	std::uniform_real_distribution<double> randomWalkDist(0.0, 1.0);
	int step = 0;
	score.resize(n);
	shifts.resize(n);
	int consecutiveNonBest = 0;
	int nSoftRestarts = 0;

	// tabu list to avoid short cycles
	IndexQueue<int> tabu(3, n); 

	for (;step < params.maxSteps; step++) {
		// Nothing to do if we got the violation to zero
		if (lessEqualThan(solution.violation(), FEASTOL))  break;

		// Pick a violated constraint at random
		assert( !violated.empty() );

		consoleDebug(3, "Repair step = {} #violated = {}", step, violated.size());

		std::uniform_int_distribution<int> uniformRow(0, violated.size()-1);
		int violrow = violated[uniformRow(rndgen)];
		const auto& row = data.mip.rows[violrow];
		const int* indices = row.idx();
		const double* coefs = row.coef();
		int cnt = row.size();

		bool violatedLessThanSense = (data.mip.sense[violrow] != 'G') &&
			greaterThan(solution.getMinAct(violrow), data.mip.rub[violrow], FEASTOL);

		consoleDebug(3, "Picked row {} violated {} act=[{},{}] rhs=[{},{}], violation={}", data.mip.rNames[violrow],
			violatedLessThanSense ? "<=" : ">=",
			solution.getMinAct(violrow), solution.getMaxAct(violrow), data.mip.rlb[violrow], data.mip.rub[violrow],
			rowViol(solution.getMinAct(violrow), solution.getMaxAct(violrow), data.mip.sense[violrow], data.mip.rlb[violrow], data.mip.rub[violrow]));

		// Retrieve list of candidates with corresponding score
		candidates.clear();
		double bestDamage = std::numeric_limits<double>::max();
		uint64_t work = 0;
		for (int pos = 0; pos < cnt; pos++) {
			int j = indices[pos];
			double coef = coefs[pos];

			/* Ignore tabu variables */
			if (tabu.has(j))  continue;

			if (domain.type(j) == 'B') {
				/* Ignore non-fixed binaries */
				if (domain.lb(j) < domain.ub(j))  continue;

				assert( (domain.ub(j) == 0.0) || (domain.lb(j) == 1.0) );
				/* Only consider variables that can reduce violation of current constraint */
				bool posDeltaX = (domain.ub(j) < 0.5);
				double mult = (posDeltaX != violatedLessThanSense) ? +1.0 : -1.0;
				if (mult*coef < 0.0)  continue;

				bool isCand;
				double damage;
				double delta = posDeltaX ? +1.0 : -1.0;
				evalShift(solution, j, delta, violrow, isCand, damage, work);

				if (isCand) {
					candidates.push_back(j);
					score[j] = damage;
					shifts[j] = delta;
					if (damage < bestDamage)  bestDamage = damage;
				}
			}
			else {
				assert( domain.type(j) != 'B' );
				assert( domain.lb(j) >= data.mip.lb[j] );
				assert( domain.ub(j) <= data.mip.ub[j] );

				/* Compute the (non necessarily integer) shift that would fix this constraint */
				double shift = 0.0;
				if (violatedLessThanSense)  shift = (data.mip.rub[violrow] - solution.getMinAct(violrow)) / coef;
				else                        shift = (data.mip.rlb[violrow] - solution.getMaxAct(violrow)) / coef;

				/* if variable is integer, round it */
				if (domain.type(j) == 'I') {
					if (shift > 0.0)  shift = ceilEps(shift, ZEROTOL);
					else              shift = floorEps(shift, ZEROTOL);
				}

				if (isNull(shift, ZEROTOL))  continue;

				/* Make sure we do not overshoot our own bounds */
				if (shift > 0.0) {
					shift = std::min(shift, data.mip.ub[j] - domain.ub(j));
				}
				else {
					shift = std::max(shift, data.mip.lb[j] - domain.lb(j));
				}
				if (isNull(shift, ZEROTOL))  continue;
				assert( domain.lb(j) + shift >= data.mip.lb[j] );
				assert( domain.ub(j) + shift <= data.mip.ub[j] );

				/* Now we can evaluate this shift */
				bool isCand;
				double damage;
				evalShift(solution, j, shift, violrow, isCand, damage, work);

				if (isCand) {
					candidates.push_back(j);
					score[j] = damage;
					shifts[j] = shift;
					if (damage < bestDamage)  bestDamage = damage;
				}
			}
		}
		work += cnt;

		consoleDebug(3, "candidates = {} / {}", candidates.size(), cnt);

		// It can happen that no flip can improve this constraint...skip it
		if (candidates.empty()) {
			if (violated.size() == 1) {
				if (tabu.empty()) {
					// there is only one violated constraint, the tabu list is empty and still we have no candidate: give up
					consoleDebug(3, "give up");
					break;
				}
				// try clearing the tabu list before giving up
				consoleDebug(3, "clear tabu list");
				tabu.clear();
			}
			continue;
		}
		assert( !candidates.empty() );

		int toFlip = data.mip.ncols;
		bool randomFlip = false;

		// If there is some damage, do a pure random walk flip with probability p
		if (isNotNull(bestDamage, FEASTOL) && randomWalkDist(rndgen) < params.p) {
			std::uniform_int_distribution<int> uniformBin(0, candidates.size()-1);
			toFlip = candidates[uniformBin(rndgen)];
			randomFlip = true;
			nRandom++;
		}
		else {
			// Otherwise, pick randomly from the candidates with the best damage
			// remove candidates with non lowest damage
			size_t first = 0;
			for (size_t k = 0; k < candidates.size(); k++) {
				if (lessEqualThan(score[candidates[k]], bestDamage, FEASTOL)) {
					candidates[first] = candidates[k];
					first++;
				}
			}
			assert( first > 0 );
			candidates.resize(first);
			std::uniform_int_distribution<int> uniformBin(0, candidates.size()-1);
			toFlip = candidates[uniformBin(rndgen)];
			if (isNull(bestDamage, FEASTOL))  nZeroDamage++;
			else                              nMinDamage++;
		}
		assert( toFlip != data.mip.ncols );

		consoleDebug(3, "Candidates after filter = {} [bestDamage = {}, randomflip = {}]", candidates.size(), bestDamage, randomFlip);

		// Apply flip/shift
		applyShift(solution, toFlip, shifts[toFlip]);
		tabu.push(toFlip);

		// update best if possible
		if (lessEqualThan(solution.violation(), bestViol, FEASTOL)) {
			consoleDebug(2, "Updated best: viol {} -> {}", bestViol, solution.violation());
			bestViol = solution.violation();
			bestNviol = violated.size();
			bestMark = solution.mark();
			consecutiveNonBest = 0;
			nSoftRestarts = 0;
		}
		else consecutiveNonBest++;

		// soft restart
		if (consecutiveNonBest >= params.maxNonImpr) {
			consoleDebug(2, "Soft restart");
			solution.undo(bestMark);
			consecutiveNonBest = 0;
			nSoftRestarts++;
		}

		solution.work += work;

		if (UserBreak)  break;
		if (nSoftRestarts >= 500)  break;
		if (gStopWatch().getElapsed() >= params.timeLimit)  break;
	}

	// Reset to best
	solution.undo(bestMark);
	consoleDebug(1, "Repair success {}. tot flips {} best flips {} [{}/{}/{}]", (int)violated.empty(), step,
		bestMark - initMark, nZeroDamage, nMinDamage, nRandom);
}


void RepairMIP::oneOpt(PropagationEngine& solution)
{
	int n = data.mip.ncols;
	int m = data.mip.nRows;
	const auto& sense = data.mip.sense;

	// nothing to do if there is some violation
	if (!solution.violatedRows().empty())  return;

	double objImpr = 0.0;
	int nShifted = 0;

	const Domain& domain = solution.getDomain();

	for (int j = 0; j < n; j++) {
		/* Cannot have unfixed variables */
		assert( equal(domain.lb(j), domain.ub(j), FEASTOL) );

		/* Cannot improve objective if variable does not appear in there */
		if (isNull(data.mip.obj[j]))  continue;

		double xj = domain.lb(j);
		bool posObj = (data.mip.objSense*data.mip.obj[j] > 0);
		bool canShift = true;
		double deltaUp = data.mip.ub[j] - xj;
		double deltaDown = xj - data.mip.lb[j];
		for (const auto [i,a]: data.mip.cols[j]) {
			/* Cannot shift variables that appear in equality constraints */
			if (data.mip.sense[i] == 'E') {
				canShift = false;
				break;
			}
			/* TODO: shift variables that appear in ranged constraints */
			if (data.mip.sense[i] == 'R') {
				canShift = false;
				break;
			}

			/* All vars are fixed, so activities should match.
			 * However, because of incremental updates, they might be a bit off.
			 */
			assert( relEqual(solution.getMinAct(i), solution.getMaxAct(i), domain.zeroTol) );

			/* Normalize to <= */
			double mult;
			double slack;
			if (data.mip.sense[i] == 'L') {
				mult = +1.0;
				slack = std::max(data.mip.rub[i] - solution.getMaxAct(i), 0.0);
			}
			else {
				mult = -1.0;
				slack = std::max(solution.getMinAct(i) - data.mip.rlb[i], 0.0);
			}

			double coef = mult*a;
			if (coef > 0) {
				deltaUp = std::min(deltaUp, slack / coef);
			}
			else {
				deltaDown = std::min(deltaDown, slack / -coef);
			}
		}

		if (canShift) {
			deltaUp = std::max(deltaUp, 0.0);
			deltaDown = std::max(deltaDown, 0.0);
			if (domain.type(j) != 'C') {
				deltaUp = floorEps(deltaUp);
				deltaDown = floorEps(deltaDown);
			}

			double deltaX = posObj ? -deltaDown : deltaUp;
			if (isNotNull(deltaX, ZEROTOL)) {
				solution.shift(j, deltaX);
				objImpr += deltaX*data.mip.obj[j];
				nShifted++;
				assert( solution.violatedRows().empty() );
			}
		}
	}

	// consoleLog("1-opt: {} flips improved objective by {}", nShifted, -objImpr);
}


struct RBranch
{
public:
	RBranch() = default;
	RBranch(int i, char s, double b) : index(i), sense(s), bound(b) {}
	int index = -1; /**< column/clique index */
	char sense; /**< 'L','U','B' for lower, upper, both respectively */
	double bound; /**< new bound */
};


/* Node data structure for repair DFS */
struct RNode
{
public:
	RNode(RBranch b, PropagationEngine::iterator _tp, PropagationEngine::iterator _stp, double _viol)
		: branch{b}, trailp(_tp), soltrailp(_stp), viol(_viol) {}
	RBranch branch;
	PropagationEngine::iterator trailp;  //< trail pointer for repair engine
	PropagationEngine::iterator soltrailp; //< trail pointer for solution engine
	double viol; //< current violation
};


/* Applies a branch to the engine */
static inline void applyBranch(const MIPData& data, const RBranch& branch, PropagationEngine& engine)
{
	assert( branch.index != -1 );

	if (branch.sense == 'U') {
		bool nodeInfeas = engine.changeUpperBound(branch.index, branch.bound);
		assert( !nodeInfeas );
	}
	else if (branch.sense == 'L') {
		bool nodeInfeas = engine.changeLowerBound(branch.index, branch.bound);
		assert( !nodeInfeas );
	}
	else {
		bool nodeInfeas = engine.changeLowerBound(branch.index, branch.bound);
		assert( !nodeInfeas );
		nodeInfeas = engine.changeUpperBound(branch.index, branch.bound);
		assert( !nodeInfeas );
	}
}


/* Apply all boundchanges from engine to solution */
static void applyChanges(const PropagationEngine& engine, PropagationEngine::iterator first, PropagationEngine::iterator last, PropagationEngine& solution)
{
	const Domain& domain = engine.getDomain();
	const Domain& solDomain = solution.getDomain();

	assert( 0 <= first );
	assert( first <= last );
	assert( last <= engine.changesEnd() );
	solution.work += (last - first);

	for (PropagationEngine::iterator itr = first; itr < last; itr++) {
		const BoundChange& bc = engine.change(itr);
		assert( bc.type != BoundChange::Type::SHIFT );
		int var = bc.var;
		assert( (0 <= var) && (var < domain.ncols()) );
		if (domain.type(var) == 'B') {
			assert( domain.lb(var) == domain.ub(var) );
			// binary case is easy. there are 3 cases:
			// 1) variable is fixed to the same value in the solution: nothing to do
			// 2) variable is fixed to opposite value in the solution: flip it
			// 3) variable is not fixed: fix it to the same value
			if (solDomain.lb(var) == solDomain.ub(var)) {
				// cases 1 and 2
				if (solDomain.lb(var) != domain.lb(var)) {
					// case 2: flip
					double delta = domain.lb(var) - solDomain.lb(var);
					assert( fabs(delta) == 1.0 );
					solution.shift(var, delta);
				}
				else {
					// case 1: nothing to do
				}
			}
			else {
				// case 3
				assert( solDomain.lb(var) == 0.0 );
				assert( solDomain.ub(var) == 1.0 );
				if (bc.type == BoundChange::Type::UPPER)  solution.changeUpperBound(var, 0.0);
				else                                      solution.changeLowerBound(var, 1.0);
			}
			assert( domain.lb(var) == solDomain.lb(var) );
			assert( domain.ub(var) == solDomain.ub(var) );
		}
		else {
			// nonbinary case is a bit more complicated
			if (lessThan(solDomain.ub(var), domain.lb(var), domain.zeroTol)) {
				// case 1: solDomain is to the left: shift to the right and fix
				double delta = domain.lb(var) - solDomain.ub(var);
				assert( delta > 0.0 );
				solution.shift(var, delta);
				assert( solDomain.ub(var) == domain.lb(var) );
				solution.changeLowerBound(var, domain.lb(var));
			}
			else if (greaterThan(solDomain.lb(var), domain.ub(var), domain.zeroTol)) {
				// case 2: solDomain is to the right: shift to the left and fix
				double delta = domain.ub(var) - solDomain.lb(var);
				assert( delta < 0.0 );
				solution.shift(var, delta);
				assert( solDomain.lb(var) == domain.ub(var) );
				solution.changeUpperBound(var, domain.ub(var));
			}
			else {
				// case3: we can take the intersection between the two intervals
				if (lessThan(domain.ub(var), solDomain.ub(var), domain.zeroTol))     solution.changeUpperBound(var, domain.ub(var));
				if (greaterThan(domain.lb(var), solDomain.lb(var), domain.zeroTol))  solution.changeLowerBound(var, domain.lb(var));
			}
			assert( lessEqualThan(domain.lb(var), solDomain.lb(var), domain.zeroTol) );
			assert( greaterEqualThan(domain.ub(var), solDomain.ub(var), domain.zeroTol) );
		}
	}
}


/* Find a repair disjunction */
std::pair<int, double> RepairMIP::findRepairMove(const PropagationEngine& engine, PropagationEngine& solution)
{
	const MIPInstance& mip = data.mip;
	const int n = mip.ncols;
	const Domain& repDomain = engine.getDomain();
	const Domain& solDomain = solution.getDomain();
	const IndexSet<int>& violated = solution.violatedRows();
	int nviolated = (int)violated.size();
	assert( nviolated > 0 );

	std::uniform_int_distribution<int> uniformRow(0, violated.size()-1);
	std::uniform_real_distribution<double> randomWalkDist(0.0, 1.0);

	int toFlip = n;
	bool randomFlip = false;
	score.resize(n);
	shifts.resize(n);

	int start = (nviolated > 1) ? uniformRow(rndgen) : 0;
	for (int i = 0; i < nviolated; i++) {
		int violrow = violated[(start + i) % nviolated];
		const auto& row = mip.rows[violrow];
		const int* indices = row.idx();
		const double* coefs = row.coef();
		int cnt = row.size();

		bool violatedLessThanSense = (mip.sense[violrow] != 'G') &&
			greaterThan(solution.getMinAct(violrow), mip.rub[violrow], FEASTOL);

		consoleDebug(3, "Picked row {} violated {} act=[{},{}] rhs=[{},{}], violation={}", mip.rNames[violrow],
				violatedLessThanSense ? "<=" : ">=",
				solution.getMinAct(violrow), solution.getMaxAct(violrow), mip.rlb[violrow], mip.rub[violrow],
				rowViol(solution.getMinAct(violrow), solution.getMaxAct(violrow), mip.sense[violrow], mip.rlb[violrow], mip.rub[violrow]));

		// Retrieve list of candidates with corresponding score
		candidates.clear();
		double bestDamage = std::numeric_limits<double>::max();
		for (int pos = 0; pos < cnt; pos++) {
			int j = indices[pos];
			double coef = coefs[pos];

			/* Domain in engine always contains the domain in the current solution */
			assert( lessEqualThan(repDomain.lb(j), solDomain.lb(j), solDomain.zeroTol) );
			assert( greaterEqualThan(repDomain.ub(j), solDomain.ub(j), solDomain.zeroTol) );

			/* Ignore variables fixed in repair engine */
			if (repDomain.lb(j) == repDomain.ub(j))  continue;

			if (solDomain.type(j) == 'B') {
				/* Ignore non-fixed binaries */
				if (solDomain.lb(j) < solDomain.ub(j))  continue;

				assert( (solDomain.ub(j) == 0.0) || (solDomain.lb(j) == 1.0) );
				/* Only consider variables that can reduce violation of current constraint */
				bool posDeltaX = (solDomain.ub(j) < 0.5);
				double mult = (posDeltaX != violatedLessThanSense) ? +1.0 : -1.0;
				if (mult*coef < 0.0)  continue;

				bool isCand;
				double damage;
				double delta = posDeltaX ? +1.0 : -1.0;
				evalShift(solution, j, delta, violrow, isCand, damage, solution.work);

				if (isCand) {
					candidates.push_back(j);
					score[j] = damage;
					shifts[j] = delta;
					if (damage < bestDamage)  bestDamage = damage;
				}
			}
			else {
				assert( solDomain.type(j) != 'B' );
				assert( repDomain.lb(j) >= mip.lb[j] );
				assert( repDomain.ub(j) <= mip.ub[j] );

				/* Compute the (non necessarily integer) shift that would fix this constraint */
				double shift = 0.0;
				if (violatedLessThanSense)  shift = (mip.rub[violrow] - solution.getMinAct(violrow)) / coef;
				else                        shift = (mip.rlb[violrow] - solution.getMaxAct(violrow)) / coef;

				/* if variable is integer, round it */
				if (solDomain.type(j) == 'I') {
					if (shift > 0.0)  shift = ceilEps(shift, ZEROTOL);
					else              shift = floorEps(shift, ZEROTOL);
				}

				if (isNull(shift, ZEROTOL))  continue;

				/* Make sure we do not overshoot our own bounds (as given by the repair engine!!!) */
				if (shift > 0.0) {
					shift = std::min(shift, repDomain.ub(j) - solDomain.ub(j));
				}
				else {
					shift = std::max(shift, repDomain.lb(j) - solDomain.lb(j));
				}
				if (isNull(shift, ZEROTOL))  continue;
				assert( solDomain.lb(j) + shift >= repDomain.lb(j) );
				assert( solDomain.ub(j) + shift <= repDomain.ub(j) );

				/* Now we can evaluate this shift */
				bool isCand;
				double damage;
				evalShift(solution, j, shift, violrow, isCand, damage, solution.work);

				if (isCand) {
					candidates.push_back(j);
					score[j] = damage;
					shifts[j] = shift;
					if (damage < bestDamage)  bestDamage = damage;
				}
			}
		}
		solution.work += cnt;

		consoleDebug(3, "candidates = {} / {}", candidates.size(), cnt);

		// It can happen that no flip can improve this constraint...skip it
		if (candidates.empty())  continue;

		// If there is some damage, do a pure random walk flip with probability p
		if (isNotNull(bestDamage, FEASTOL) && randomWalkDist(rndgen) < params.p) {
			std::uniform_int_distribution<int> uniformBin(0, candidates.size()-1);
			toFlip = candidates[uniformBin(rndgen)];
			randomFlip = true;
		}
		else {
			// Otherwise, pick randomly from the candidates with the best damage
			// remove candidates with non lowest damage
			size_t first = 0;
			for (size_t k = 0; k < candidates.size(); k++) {
				if (lessEqualThan(score[candidates[k]], bestDamage, FEASTOL)) {
					candidates[first] = candidates[k];
					first++;
				}
			}
			assert( first > 0 );
			candidates.resize(first);
			std::uniform_int_distribution<int> uniformBin(0, candidates.size()-1);
			toFlip = candidates[uniformBin(rndgen)];
		}
		consoleDebug(3, "Candidates after filter = {} [bestDamage = {}, randomflip = {}]", candidates.size(), bestDamage, randomFlip);

		// Found flip/shift to reduce violation
		assert( toFlip != n );
		assert( fabs(shifts[toFlip]) > 0.0 );

		return {toFlip, shifts[toFlip]};
	}

	// No repair move found!
	return {-1, 0.0};
}


/* Turn repair move into a repair disjunction */
static std::vector<RBranch> repairDisjunction(const PropagationEngine& engine, PropagationEngine& solution, int toFlip, double shift)
{
	assert( toFlip != -1 );
	assert( shift != 0.0 );
	const Domain& repDomain = engine.getDomain();
	const Domain& solDomain = solution.getDomain();

	// binary case: easy
	if (solDomain.type(toFlip) == 'B') {
		RBranch up{toFlip, 'L', 1.0};
		RBranch down{toFlip, 'U', 0.0};
		if (shift > 0.0) {
			return {up, down};
		}
		else {
			return {down, up};
		}
	}

	// We are left with the non-binary case
	// compute shifted domain
	double shiftedLB = solDomain.lb(toFlip) + shift;
	double shiftedUB = solDomain.ub(toFlip) + shift;
	double leftGap = shiftedLB - repDomain.lb(toFlip);
	double rightGap = repDomain.ub(toFlip) - shiftedUB;
	assert( (leftGap > 0.0) || (rightGap > 0.0) );

	// TODO: for integer variables we don't need overlap!
	if (leftGap <= rightGap) {
		RBranch left{toFlip, 'U', shiftedUB};
		RBranch right{toFlip, 'L', shiftedUB};
		return {left, right};
	}
	else {
		RBranch right{toFlip, 'L', shiftedLB};
		RBranch left{toFlip, 'U', shiftedLB};
		return {right, left};
	}
}


/** Perform repair DFS from a given infeasible partial solution */
void RepairMIP::search(PropagationEngine& solution, PropagationEngine& engine)
{
	const MIPInstance& mip = data.mip;
	const Domain& repDomain = engine.getDomain();
	const Domain& solDomain = solution.getDomain();

	//int n = mip.ncols;
	assert( engine.violatedRows().empty() );
	assert( !solution.violatedRows().empty() );

	std::vector<RNode> nodes;
	std::mt19937_64 rndgen;
	rndgen.seed(params.seed);

	// push root node
	PropagationEngine::iterator startMark = engine.mark();
	PropagationEngine::iterator solStartMark = solution.mark();
	double initViol = solution.violation();
	double bestViol = initViol;
	nodes.emplace_back(RBranch{}, startMark, solStartMark, initViol);
	int nodecnt = 0;
	int nonImpr = 0;

	auto branch2str = [&](const RBranch& br) {
		const char* sense;
		if (br.sense == 'B')       sense = "=";
		else if (br.sense == 'L')  sense = ">=";
		else                       sense = "<=";
		return fmt::format("{} {} {}", mip.cNames[br.index], sense, br.bound);
	};

	// DFS
	while (!nodes.empty()) {
		// pop node
		RNode node = nodes.back();
		nodes.pop_back();
		// backtrack
		engine.undo(node.trailp);
		solution.undo(node.soltrailp);
		const RBranch& branch = node.branch;
		// apply branch (if any)
		if (branch.index != -1) {
			consoleDebug(2, "Apply repair branching {}", branch2str(branch));
			applyBranch(data, branch, engine);
		}

		nodecnt++;

		bool nodeInfeas = (!engine.violatedRows().empty());

		// propagate
		if (params.repairPropagate) {
			nodeInfeas = engine.propagate(false);
#ifdef DEBUG_LOG
			if (DEBUG_LEVEL >= 3)  printChangesSinceMark(engine, node.trailp);
#endif
			consoleDebug(2, "Bound changes after repair propagate: {} [infeas={}]", engine.mark() - node.trailp, nodeInfeas);
		}

		if (!engine.violatedRows().empty()) {
			consoleDebug(3, "Repair backtrack");
			continue;
		}

		// apply branch (and its implications) to solution
		applyChanges(engine, node.trailp, engine.changesEnd(), solution);
		consoleDebug(3, "Repair node = {} #violated = {}", nodecnt, solution.violatedRows().size());

		// Did we repair the solution? If so we are done
		if (solution.violatedRows().empty()) {
			break;
		}

		// Check if we should backtrack for lack of recent progress
		double currViol = solution.violation();
		if (lessEqualThan(currViol, bestViol, FEASTOL)) {
			consoleDebug(2, "Updated best: viol {} -> {}", bestViol, currViol);
			bestViol = currViol;
			nonImpr = 0;
		}
		else nonImpr++;

		if ((nonImpr > params.maxNonImpr) && (!nodes.empty())) {
			// Backtrack to the best open node as far as violation is concerned
			// (in case of ties, go to the deepest)
			bestViol = std::numeric_limits<double>::max();
			size_t bestPos = nodes.size();
			for (size_t n = 0; n < nodes.size(); n++) {
				const RNode& node = nodes[n];
				if (node.viol <= bestViol) {
					bestViol = node.viol;
					bestPos = n;
				}
			}

			size_t oldOpenNodes = nodes.size();
			assert( (bestPos+1) <= oldOpenNodes );
			nodes.erase(nodes.begin() + bestPos + 1, nodes.end());
			consoleDebug(3, "Repair backtrack for lack of progress. Open nodes {} -> {}", oldOpenNodes, nodes.size());
			nonImpr = 0;
			continue;
		}

		// find repair move
		auto [toFlip, shift] = findRepairMove(engine, solution);

		if (toFlip != -1) {
			// turn move this into a repair disjunction
			std::vector<RBranch> branches = repairDisjunction(engine, solution, toFlip, shift);
			assert( !branches.empty() );
			consoleDebug(3, "{}-way repair branch", branches.size());
			// add them in reverse order (this is a stack after all)
			for (auto itr = branches.rbegin(); itr != branches.rend(); ++itr) {
				nodes.emplace_back(*itr, engine.mark(), solution.mark(), currViol);
			}
		}
		else {
			// Nothing to branch on: backtrack
			continue;
		}

		// termination criteria
		if ( (nodecnt >= params.maxSteps) ||
			 (gStopWatch().getElapsed() >= params.timeLimit) ||
			 dominiqs::UserBreak ) {
			consoleDebug(2, "Repair limits reached");
			break;
		}
	}

	// cleanup
	engine.undo(startMark);

	bool success = solution.violatedRows().empty();
	if (!success) {
		// there might still be a node with a reduced violation w.r.t. the initial one
		// find it and backtrack there
		bestViol = std::numeric_limits<double>::max();
		size_t bestPos = nodes.size();
		for (size_t n = 0; n < nodes.size(); n++) {
			const RNode& node = nodes[n];
			if (node.viol < bestViol) {
				bestViol = node.viol;
				bestPos = n;
			}
		}

		if (lessThan(bestViol, initViol, FEASTOL))  solution.undo(nodes[bestPos].soltrailp);
		else                                        solution.undo(solStartMark);
	}
	int nflips = solution.changesEnd() - solStartMark;

	consoleDebug(1, "Repair success {}. tot nodes {} best flips {}", (int)success, nodecnt, nflips);
}
