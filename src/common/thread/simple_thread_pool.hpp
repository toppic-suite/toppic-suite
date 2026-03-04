//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_COMMON_THREAD_SIMPLE_THREAD_POOL_HPP_
#define TOPPIC_COMMON_THREAD_SIMPLE_THREAD_POOL_HPP_

#include <atomic>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <unordered_map>
#include <vector>

namespace toppic {

class SimpleThreadPool {
 public:
  SimpleThreadPool(int thread_num);
  ~SimpleThreadPool();

  SimpleThreadPool(const SimpleThreadPool&) = delete;
  SimpleThreadPool& operator=(const SimpleThreadPool&) = delete;
  SimpleThreadPool(SimpleThreadPool&&) = delete;
  SimpleThreadPool& operator=(SimpleThreadPool&&) = delete;

  // Adds task to a task queue.
  void enqueue(std::function<void()> f);

  // Shut down the pool.
  void shutDown();

  size_t getQueueSize() const;
  size_t getThreadNum() const;

  int getIdleThreadNum() const { return idle_thread_num_; }

  int getId(std::thread::id thread_id) const;

 private:
  // Thread pool storage.
  std::vector<std::thread> threads_;

  // Maps thread::id to index in threads_ for O(1) lookup in getId().
  std::unordered_map<std::thread::id, int> thread_id_map_;

  // Queue to keep track of incoming tasks.
  std::queue<std::function<void()> > tasks_;

  // Mutex protecting tasks_, threads_, and thread_id_map_.
  mutable std::mutex tasks_mutex_;

  // Condition variable.
  std::condition_variable condition_;

  // Indicates that pool needs to be shut down.
  std::atomic<bool> terminate_;

  // Idle thread number.
  std::atomic<int> idle_thread_num_;

  // Function that will be invoked by our threads.
  void invoke();
};

using SimpleThreadPoolPtr = std::shared_ptr<SimpleThreadPool>;

}  // namespace toppic

#endif
