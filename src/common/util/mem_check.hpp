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

#ifndef TOPPIC_COMMON_UTIL_MEMCHECK_HPP_
#define TOPPIC_COMMON_UTIL_MEMCHECK_HPP_

#include <string>
#include <utility>
#include <vector>

namespace toppic {

namespace mem_check {

int getMaxThreads(const std::string& app_name);

// The thread number one protein sequence filtering stage of toppic or topmg
// will use: the stage's own memory budget applied to the memory available
// now, lowered to the requested thread number. Call it where the stage
// starts, and pass the stage's banner name: the returned number is the one
// the stage uses, and the function reports it when the available memory has
// moved it away from the number checkThreadNum reported for that stage at
// the start of the run.
int getStageThreadNum(const std::string& stage_name, int thread_num);

// The thread number of the graph alignment of topmg and of toppic's --topmg:
// a thread holds one alignment's DP table, consistent pairs and proteoform
// graph, whose size grows with the table's layer number, (unexpected shift
// number + 1) x (variable modification number + 1). The budget per thread
// is applied to the memory available now and lowered to thread_num; the
// function reports the number when it is lower.
int getGraphAlignThreadNum(int thread_num, int layer_num);

double getTotalMemInGb();

double getAvailMemInGb();

// the available memory as it was when the run started (sampled once, at the
// first call), the figure checkThreadNum prints and reports its per-stage
// thread numbers from
double getStartupAvailMemInGb();

// The protein sequence filtering steps of toppic or topmg whose thread number
// the memory available at the start of the run lowers below thread_num, each
// with the number it will use; empty when no step is lowered, and for a
// program that does not control filtering thread numbers itself. The console
// warning of checkThreadNum and the GUI dialogs are both built from it, so
// the two cannot describe the same run differently.
std::vector<std::pair<std::string, int>> getLimitedFilterStages(
    const std::string& prog, int thread_num);

bool checkThreadNum(int thread_num, const std::string& prog);

}  // namespace mem_check

}  // namespace toppic

#endif
