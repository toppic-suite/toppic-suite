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
#include <vector>

namespace toppic {

using ThreadPtr = std::shared_ptr<std::thread>;

class ToppicThread {
 public:
  ToppicThread(int id, ThreadPtr thread_ptr):
      id_(id),
      thread_ptr_(thread_ptr) {}
  int getId() const { return id_; }
  const ThreadPtr& getThreadPtr() const { return thread_ptr_; }

 private:
  int id_;
  ThreadPtr thread_ptr_;
};

using ToppicThreadPtr = std::shared_ptr<ToppicThread>;
using ToppicThreadPtrVec = std::vector<ToppicThreadPtr>;

class SimpleThreadPool {
 public:
  // Constructor.
  SimpleThreadPool(int thread_num);

  // Destructor.
  ~SimpleThreadPool();

  SimpleThreadPool(const SimpleThreadPool&) = delete;
  SimpleThreadPool& operator=(const SimpleThreadPool&) = delete;
  SimpleThreadPool(SimpleThreadPool&&) = delete;
  SimpleThreadPool& operator=(SimpleThreadPool&&) = delete;

  // Adds task to a task queue.
  void enqueue(std::function<void()> f);

  // Shut down the pool.
  void shutDown();

  size_t getQueueSize() const {
    std::unique_lock<std::mutex> lock(tasks_mutex_);
    return tasks_.size();
  }

  size_t getThreadNum() const {
    std::unique_lock<std::mutex> lock(tasks_mutex_);
    return thread_ptr_vec_.size();
  }

  int getIdleThreadNum() const { return idle_thread_num_; }

  int getId(std::thread::id thread_id);

 private:
  // Thread pool storage.
  ToppicThreadPtrVec thread_ptr_vec_;

  // Queue to keep track of incoming tasks.
  std::queue<std::function<void()> > tasks_;

  // Task queue mutex.
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
