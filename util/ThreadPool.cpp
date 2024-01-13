#include "ThreadPool.h"

ThreadPool::ThreadPool(int num_threads):
m_num_threads(num_threads) {
}

ThreadPool::~ThreadPool() {
}

void ThreadPool::Start() {
   for (uint32_t ii = 0; ii < m_num_threads; ++ii) {
      m_threads.emplace_back(thread(&ThreadPool::ThreadLoop, this));
   }
}

void ThreadPool::ThreadLoop() {
   while (true) {
      function<void(int)> job;
      int jid;
      {
         unique_lock<mutex> lock(m_queue_mutex);
         m_mutex_condition.wait(lock, [this] {
            return !m_jobs.empty() || m_should_terminate;
         });
         if (m_should_terminate) {
            return;
         }
         job = m_jobs.front();
         m_jobs.pop();
         jid = m_jids.front();
         m_jids.pop();
      }
      job(jid);
      {
         unique_lock<mutex> lock(m_queue_mutex);
         m_n_jobs--;
      }
   }
}

void ThreadPool::QueueJob(const function<void(int)>& job, int jid) {
   {
      unique_lock<mutex> lock(m_queue_mutex);
      m_jobs.push(job);
      m_jids.push(jid);
      m_n_jobs++;
   }
   m_mutex_condition.notify_one();
}

bool ThreadPool::busy() {
   bool poolbusy;
   {
      unique_lock<mutex> lock(m_queue_mutex);
      poolbusy = m_n_jobs != 0;
   }
   return poolbusy;
}

void ThreadPool::Stop() {
   {
      unique_lock<mutex> lock(m_queue_mutex);
      m_should_terminate = true;
   }
   m_mutex_condition.notify_all();
   for (thread& active_thread : m_threads) {
      active_thread.join();
   }
   m_threads.clear();
}
