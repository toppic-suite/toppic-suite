//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "common/thread/simple_thread_pool.hpp"

namespace toppic {

SimpleThreadPool::SimpleThreadPool(int thread_num) :
    terminate_(false), stopped_(false) {
      for (int i = 0; i < thread_num; i++) {
        ThreadPtr thread_ptr = std::make_shared<boost::thread>(&SimpleThreadPool::Invoke, this);
        ToppicThreadPtr toppic_thread_ptr = std::make_shared<ToppicThread>(i, thread_ptr);
        thread_ptr_vec_.emplace_back(toppic_thread_ptr);
      }
    }

void SimpleThreadPool::Enqueue(std::function<void()> f) {
  // Scope based locking.
  {
    // Put unique lock on task mutex.
    boost::unique_lock<boost::mutex> lock(tasks_mutex_);

    // Push task into queue.
    tasks_.push(f);
  }

  // Wake up one thread.
  condition_.notify_one();
}

void SimpleThreadPool::Invoke() {
  std::function<void()> task;
  while (true) {
    // Scope based locking.
    {
      // Put unique lock on task mutex.
      boost::unique_lock<boost::mutex> lock(tasks_mutex_);

      // Wait until queue is not empty or termination signal is sent.
      condition_.wait(lock, [this]{ return !tasks_.empty() || terminate_; });

      // If termination signal received and queue is empty then exit else continue clearing the queue.
      if (terminate_ && tasks_.empty()) {
        return;
      }

      // Get next task in the queue.
      task = tasks_.front();

      // Remove it from the queue.
      tasks_.pop();
    }

    // Execute the task.
    task();
  }
}

void SimpleThreadPool::ShutDown() {
  // Scope based locking.
  {
    // Put unique lock on task mutex.
    boost::unique_lock<boost::mutex> lock(tasks_mutex_);

    // Set termination flag to true.
    terminate_ = true;
  }

  // Wake up all threads.
  condition_.notify_all();

  // Join all threads.
  for (ToppicThreadPtr top_thread_ptr : thread_ptr_vec_) {
    ThreadPtr thread_ptr = top_thread_ptr->getThreadPtr();
    if (thread_ptr->joinable()) thread_ptr->join();
  }

  // Empty workers vector.
  thread_ptr_vec_.empty();

  // Indicate that the pool has been shut down.
  stopped_ = true;
}

int SimpleThreadPool::getId(boost::thread::id thread_id) {
  for (size_t i = 0; i < thread_ptr_vec_.size(); i++) {
    if (thread_ptr_vec_[i]->getThreadPtr()->get_id() == thread_id) {
      return thread_ptr_vec_[i]->getId();
    }
  }
  LOG_ERROR("Thread Error: no writer found");
  exit(EXIT_FAILURE);
}

// Destructor
SimpleThreadPool::~SimpleThreadPool() {
  if (!stopped_) {
    ShutDown();
  }
}

}  // namespace toppic
