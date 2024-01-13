#pragma once

#include <vector>
#include <queue>
#include <thread>
#include <functional>
#include <mutex>

using std::vector;
using std::thread;
using std::function;
using std::mutex;
using std::unique_lock;
using std::condition_variable;
using std::queue;

class ThreadPool {
public:
   ThreadPool(int num_threads);
   ~ThreadPool();
   void Start();
   void QueueJob(const function<void(int)>& job, int jid);
   bool busy();
   void Stop();

private:
   void ThreadLoop();

   bool m_should_terminate = false;
   int m_n_jobs = 0;
   int m_num_threads;
   vector<thread> m_threads;
   mutex m_queue_mutex;
   condition_variable m_mutex_condition;
   queue<function<void(int)>> m_jobs;
   queue<int> m_jids;
};
