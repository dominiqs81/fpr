#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <cassert>
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>


class ThreadPool
{
public:
	using Task = std::function<void()>;
	/* Thread Pool with numThreads threads */
	ThreadPool(int numThreads) : stopped{false}, working{0}
	{
		for (int i = 0; i< numThreads; ++i)  workers.emplace_back(&ThreadPool::workerLoop, this);
	}
	/* Destructor (waiting for threads to finish) */
	~ThreadPool()
	{
		/* Set stop to false */
		{
			std::unique_lock<std::mutex> lock{m};
			stopped = true;
		}
		/* Notify all and join */
		cv_tasks.notify_all();
		for (auto& w: workers)  w.join();
	}
	/* Enqueue a task to be performed by the thread pool */
	template<class Func>
	void enqueue(Func&& f)
	{
		/* don't allow enqueueing after stopping the pool! */
		if (stopped) throw std::runtime_error("enqueue on stopped ThreadPool");

		/* push task to the queue */
		std::unique_lock<std::mutex> lock{m};
		tasks.push(Task{f});

		/* signal workers there is something to do */
		cv_tasks.notify_one();
	}
	/* Wait for all pending tasks to be done */
	void wait()
	{
		std::unique_lock<std::mutex> lock{m};
		cv_done.wait(lock, [this](){ return (tasks.empty() && (working == 0)); });
	}
private:
	/* Main loop executed by each thread in the pool */
	void workerLoop()
	{
		while (true) {
			/* Wait for a task to process or for the stop signal */
			std::unique_lock<std::mutex> lock{m};
			cv_tasks.wait(lock, [this](){ return (stopped || (!tasks.empty())); });

			/* Stop if signalled */
			if (stopped)  break;

			/* Pop task from queue */
			++working;
			assert( working <= (int)workers.size() );
			assert( !tasks.empty() );
			std::function<void()> task(tasks.front());
			tasks.pop();

			/* Release lock to do the work! */
			lock.unlock();

			/* Execute task */
			task();

			/* Grab lock again to decrease working counter */
			lock.lock();
			--working;
			assert( working >= 0);
			cv_done.notify_one();
		}
	}
	std::vector<std::thread> workers; //< worker threads
	std::queue<Task> tasks; //< task queue
	bool stopped;
	int working = 0;
	// synchronization
	std::mutex m;
	std::condition_variable cv_tasks;
	std::condition_variable cv_done;
};

#endif // THREAD_POOL_H
