#ifndef REPAIR_H
#define REPAIR_H

#include <cstdint>
#include <random>
#include "fixproprep/mip.h"
#include "fixproprep/propagation.h"


// Solution Repair for MIP
class RepairMIP
{
public:
	RepairMIP(const MIPData& _data, const Params& _params);
	void walk(PropagationEngine& solution);
	void search(PropagationEngine& solution, PropagationEngine& engine);
	void oneOpt(PropagationEngine& solution);
private:
	// data
	const MIPData& data;
	const Params& params;
	// state
	std::mt19937_64 rndgen;
	// work
	std::vector<int> candidates;
	std::vector<double> score;
	std::vector<double> shifts;
	// helper
	std::pair<int, double> findRepairMove(const PropagationEngine& engine, PropagationEngine& solution);
};


#endif /* REPAIR_H */
