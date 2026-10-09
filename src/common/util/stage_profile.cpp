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

#include "common/util/stage_profile.hpp"

#include <cstdlib>
#include <fstream>
#include <map>
#include <mutex>
#include <string>

namespace toppic {

namespace stage_profile {

namespace {

struct Acc {
  double seconds = 0;
  long count = 0;
};

std::mutex profile_mutex;
std::map<std::string, Acc> accs;

const char* profileFile() {
  static const char* file = std::getenv("TOPPIC_PROFILE");
  return file;
}

}  // namespace

bool enabled() { return profileFile() != nullptr; }

void add(const std::string& name, double seconds, long count) {
  if (!enabled()) {
    return;
  }
  std::lock_guard<std::mutex> lock(profile_mutex);
  Acc& a = accs[name];
  a.seconds += seconds;
  a.count += count;
}

void report(const std::string& stage) {
  if (!enabled()) {
    return;
  }
  std::lock_guard<std::mutex> lock(profile_mutex);
  std::ofstream out(profileFile(), std::ios::app);
  for (const auto& kv : accs) {
    out << stage << "\t" << kv.first << "\t" << kv.second.seconds << "\t"
        << kv.second.count << "\n";
  }
  accs.clear();
}

}  // namespace stage_profile

}  // namespace toppic
