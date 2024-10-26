/**
 * @brief Worker Data structures and methods
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * Copyright 2021 Domenico Salvagnin
 */

#ifndef WORKER_H
#define WORKER_H

#include "mip.h"
#include "propagation.h"
#include <utils/mipmodel.h>
#include <mutex>


/* Worker (i.e., thread-local) data */
struct WorkerData
{
public:
	WorkerData(const MIPData& _mipdata, const PropagationEngine& _engine, bool& _forceStop)
		: mipdata{_mipdata}, engine{_engine}, repairEngine{_engine}, forceStop{_forceStop}
	{
		// set sense in pool
		solpool.setObjSense(mipdata.mip.objSense);
		// clone LP relaxation
		lp = mipdata.lp->clone();
	}
	// data
	const MIPData& mipdata;
	PropagationEngine engine;
	PropagationEngine repairEngine;
	uint64_t work = 0;
	bool& forceStop;
	SolutionPool solpool;
	MIPModelPtr lp; //< LP relaxation
};

using WorkerDataPtr = std::shared_ptr<WorkerData>;


/* Manager for worker data instances (basically a thread-safe arena) */
struct WorkerDataManager
{
public:
	WorkerDataManager(MIPData& _mipdata, const PropagationEngine& _engine, bool& _forceStop)
		: mipdata{_mipdata}, engine{_engine}, forceStop{_forceStop} {}
	WorkerDataPtr get()
	{
		std::lock_guard lock{m};
		if (unused.empty())  return std::make_shared<WorkerData>(mipdata, engine, forceStop);
		else {
			WorkerDataPtr w = unused.back();
			unused.pop_back();
			return w;
		}
	}
	void release(WorkerDataPtr w)
	{
		std::lock_guard lock{m};
		// merge local solution pool into global one
		mipdata.solpool.merge(w->solpool);
		// keep track of work
		maxWork = std::max(maxWork, w->work);
		w->work = 0;
		unused.push_back(w);
		// if we have a solution, we can set the force stop flag
		if (mipdata.solpool.hasFeas())  forceStop = true;
	}
	uint64_t maxWork = 0;
private:
	MIPData& mipdata;
	const PropagationEngine& engine;
	bool& forceStop;
	std::vector<WorkerDataPtr> unused;
	std::mutex m;
};

#endif /* WORKER_H */
