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

#ifndef TOPPIC_FILTER_MASSMATCH_INDEX_LOAD_GATE_HPP_
#define TOPPIC_FILTER_MASSMATCH_INDEX_LOAD_GATE_HPP_

#include <condition_variable>
#include <mutex>

namespace toppic {

// Bounds how many filter worker threads deserialize their block index files
// at the same time (a C++17 stand-in for std::counting_semaphore). Unbounded
// concurrent loads saturate the disk just as well, but every load in flight
// holds its archive stream plus the partially built vectors, and letting all
// 12-16 workers load at once raised a run's peak resident memory by ~4 GB on
// a 14 GB machine (toppic_test experiment 66). Four concurrent loads keep an
// NVMe near its parallel-read throughput at a bounded transient cost. The
// gate is shared by all filter stages; within one stage it is the only
// load-concurrency limit once the stage's first block has warmed
// Boost.Serialization's singletons (see zero_ptm_filter.cpp).
class IndexLoadGate {
 public:
  static const int kMaxConcurrentLoads = 4;

  void acquire() {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this] { return in_flight_ < kMaxConcurrentLoads; });
    in_flight_++;
  }

  void release() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      in_flight_--;
    }
    cv_.notify_one();
  }

  // RAII holder for one load
  class Slot {
   public:
    explicit Slot(IndexLoadGate& gate) : gate_(gate) { gate_.acquire(); }
    ~Slot() { gate_.release(); }
    Slot(const Slot&) = delete;
    Slot& operator=(const Slot&) = delete;

   private:
    IndexLoadGate& gate_;
  };

 private:
  std::mutex mutex_;
  std::condition_variable cv_;
  int in_flight_ = 0;
};

inline IndexLoadGate index_load_gate;

}  // namespace toppic

#endif
