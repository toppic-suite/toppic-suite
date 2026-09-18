// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_COMMON_THREAD_ORDERED_TASK_RUNNER_HPP_
#define TOPPIC_COMMON_THREAD_ORDERED_TASK_RUNNER_HPP_

#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <utility>

#include "common/thread/simple_thread_pool.hpp"

namespace toppic {

// Runs tasks on a SimpleThreadPool and delivers their results on the
// submitting thread in submission order, so a stage can compute in parallel
// while writing its output exactly as a sequential loop would. At most
// max_pending tasks are in flight: submit() blocks (consuming finished
// results) once that many are outstanding, which bounds memory. An exception
// thrown by a task is rethrown from submit() or finish().
template <typename Result>
class OrderedTaskRunner {
 public:
  OrderedTaskRunner(SimpleThreadPoolPtr pool_ptr, std::size_t max_pending)
      : pool_ptr_(std::move(pool_ptr)), max_pending_(max_pending) {}

  // compute runs on a pool thread; consume runs on the submitting thread,
  // after the consumers of every earlier task.
  void submit(std::function<Result()> compute,
              std::function<void(Result&)> consume) {
    while (pending_.size() >= max_pending_) {
      consumeFront();
    }
    auto promise_ptr = std::make_shared<std::promise<Result>>();
    pending_.emplace_back(promise_ptr->get_future(), std::move(consume));
    pool_ptr_->enqueue([promise_ptr, compute = std::move(compute)]() {
      try {
        promise_ptr->set_value(compute());
      } catch (...) {
        promise_ptr->set_exception(std::current_exception());
      }
    });
  }

  // Wait for the outstanding tasks and consume their results in order.
  void finish() {
    while (!pending_.empty()) {
      consumeFront();
    }
  }

 private:
  using Pending = std::pair<std::future<Result>, std::function<void(Result&)>>;

  void consumeFront() {
    Pending front = std::move(pending_.front());
    pending_.pop_front();
    Result result = front.first.get();
    front.second(result);
  }

  SimpleThreadPoolPtr pool_ptr_;
  std::size_t max_pending_;
  std::deque<Pending> pending_;
};

}  // namespace toppic

#endif
