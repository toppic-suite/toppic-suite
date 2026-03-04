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

#include <stdexcept>

#include "common/util/logger.hpp"
#include "common/thread/simple_thread_pool.hpp"

namespace toppic {

SimpleThreadPool::SimpleThreadPool(int thread_num) :
    terminate_(false), idle_thread_num_(0) {
  threads_.reserve(thread_num);
  for (int i = 0; i < thread_num; i++) {
    threads_.emplace_back(&SimpleThreadPool::invoke, this);
    thread_id_map_[threads_.back().get_id()] = i;
  }
}

size_t SimpleThreadPool::getQueueSize() const {
  std::lock_guard<std::mutex> lock(tasks_mutex_);
  return tasks_.size();
}

size_t SimpleThreadPool::getThreadNum() const {
  std::lock_guard<std::mutex> lock(tasks_mutex_);
  return threads_.size();
}

void SimpleThreadPool::enqueue(std::function<void()> f) {
  // Scope based locking.
  {
    // Put lock on task mutex.
    std::lock_guard<std::mutex> lock(tasks_mutex_);

    // Reject tasks after shutdown to avoid silent data loss.
    if (terminate_) {
      LOG_ERROR("SimpleThreadPool::enqueue called after shutDown; task dropped.");
      return;
    }

    // Push task into queue.
    tasks_.push(std::move(f));
  }

  // Wake up one thread.
  condition_.notify_one();
}

void SimpleThreadPool::invoke() {
  std::function<void()> task;
  while (true) {
    // Scope based locking.
    {
      // Put unique lock on task mutex.
      std::unique_lock<std::mutex> lock(tasks_mutex_);
      ++idle_thread_num_;

      // Wait until queue is not empty or termination signal is sent.
      // The predicate is evaluated while holding the lock, so accessing
      // tasks_ and terminate_ here is safe.
      condition_.wait(lock, [this]{ return !tasks_.empty() || terminate_; });

      // If termination signal received and queue is empty then exit else continue clearing the queue.
      if (terminate_ && tasks_.empty()) {
        --idle_thread_num_;
        return;
      }

      // Get next task in the queue.
      task = tasks_.front();

      // Remove it from the queue.
      tasks_.pop();
      --idle_thread_num_;
    }

    // Execute the task.
    task();
  }
}

void SimpleThreadPool::shutDown() {
  // Scope based locking.
  {
    // Put lock on task mutex.
    std::lock_guard<std::mutex> lock(tasks_mutex_);

    // Set termination flag to true.
    terminate_ = true;
  }

  // Wake up all threads.
  condition_.notify_all();

  // Join all threads.
  for (std::thread& t : threads_) {
    if (t.joinable()) t.join();
  }

  // Empty workers vector.
  {
    std::lock_guard<std::mutex> lock(tasks_mutex_);
    threads_.clear();
    thread_id_map_.clear();
  }
}

int SimpleThreadPool::getId(std::thread::id thread_id) const {
  std::lock_guard<std::mutex> lock(tasks_mutex_);
  auto it = thread_id_map_.find(thread_id);
  if (it != thread_id_map_.end()) return it->second;
  LOG_ERROR("Thread Error: thread id " << thread_id << " cannot be found!");
  throw std::runtime_error("Thread id not found in pool");
}

SimpleThreadPool::~SimpleThreadPool() {
  if (!threads_.empty()) {
    shutDown();
  }
}

}  // namespace toppic
