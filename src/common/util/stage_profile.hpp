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

#ifndef TOPPIC_COMMON_UTIL_STAGE_PROFILE_HPP_
#define TOPPIC_COMMON_UTIL_STAGE_PROFILE_HPP_

#include <chrono>
#include <string>

namespace toppic {

// Research timers (toppic_test experiment 150): named accumulators of
// thread time and counts inside the filter and search stages, summed over
// the threads and appended to the file named by the environment variable
// TOPPIC_PROFILE at the end of a stage (report). Without the variable the
// timers cost one branch per call. Not for release builds' output.
namespace stage_profile {

bool enabled();

// adds seconds (and one call) to the named accumulator
void add(const std::string& name, double seconds, long count = 1);

// one line per accumulator with the prefix (the stage's name), then they
// are cleared
void report(const std::string& stage);

class Timer {
 public:
  explicit Timer(const char* name) : name_(name) {
    if (enabled()) {
      start_ = std::chrono::steady_clock::now();
    }
  }
  ~Timer() {
    if (enabled()) {
      double s = std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start_).count();
      add(name_, s);
    }
  }
  Timer(const Timer&) = delete;
  Timer& operator=(const Timer&) = delete;

 private:
  const char* name_;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace stage_profile

}  // namespace toppic

#endif
