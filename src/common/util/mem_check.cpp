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

#include "common/util/mem_check.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <thread>
#include <utility>
#include <vector>

#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW32__) || \
    defined(__MINGW64__)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>

#include <cstdint>
#endif

#include "common/util/logger.hpp"

namespace toppic {

namespace mem_check {

const std::map<std::string, double> memory_per_thread_list{
    // topfd memory requirement per thread: about 0.15 gb
    {"topfd", 0.3},
    // topdia memory requirement per thread: about 0.15 gb
    {"topdia", 0.3},
    // topindex memory requirement per thread: about 0.4 gb
    {"topindex", 0.5},
    // toppic memory requirement per thread: about 0.75 gb
    {"toppic", 1.0},
    // topmg memory requirement per thread: about 0.75 gb
    {"topmg", 1.0},
    // zero or one shift filter memory requirement per thread: about 0.75 gb
    {"zero_one_shift_filter", 1.0},
    // the zero-shift filter's four block index files are about 0.95 gb per
    // thread, and with the block loads running in parallel every worker
    // holds its full block at once (serialized loading used to keep several
    // workers empty-handed at any moment, which is what the shared 1.0 gb
    // estimate was calibrated against): budgeting the true residency keeps
    // the stage's peak memory inside the machine (toppic_test experiment 66
    // measured 14.0 gb peak at 12 threads on a 14 gb machine without this).
    // The variable PTM filter reads the same four files per block and
    // budgets with this entry too (toppic_test report 0077_bfb05226: 13.8
    // gb peak and swapping at 13 threads on a 14 gb machine with the 1.0 gb
    // entry)
    {"zero_shift_filter", 1.25},
    // diag filter memory requirement per thread: about 0.4 gb
    {"diag_filter", 0.5},
    // graph alignment memory requirement per thread, for a DP table of
    // kGraphAlignLayers layers and scaled by the layer number
    // (getGraphAlignThreadNum). A thread holds one alignment: its DP table
    // (26 bytes for each of proteoform graph vertices x spectrum graph
    // vertices x layers), its consistent pairs and the proteoform graph's
    // distance lists. On the first 500 spectra of 480_F3_02 the largest was
    // 329 MiB at 8 layers (a 2,001-vertex proteoform graph against a
    // 592-vertex spectrum graph; 235 MiB of it the DP table), 0.32 GiB, and
    // 686 MiB at 18 layers (-S 5, --topmg-num-shift 2), 0.67 GiB; 0.5 gb at
    // 8 layers leaves room for spectra with more peaks (toppic_test
    // experiment 77)
    {"graph_align", 0.5}};

// the layer number the graph_align entry is for: one unexpected shift (TopMG's
// default) and three variable modifications (toppic's default -S), (1 + 1) x
// (3 + 1)
const int kGraphAlignLayers = 8;

// The protein sequence filtering stages of toppic and topmg budget their
// thread number per stage (see toppic_process.cpp and topmg_process.cpp):
// each one divides the memory available when it starts by its own entry of
// memory_per_thread_list, so the stages of one run do not all get the same
// number. These two lists are what checkThreadNum reports from and what
// getStageThreadNum computes from, so the message and the stages cannot name
// different budgets. A stage missing here is a stage whose reduced thread
// number the user is never told about, and a wrong memory_per_thread_list
// key is a reported number no stage uses. Keep them in step with the
// drivers; the names are the stage banners the drivers print.
const std::map<std::string, std::string> filter_stage_mem_entry_list{
    {"Zero unexpected shift filtering", "zero_shift_filter"},
    {"Variable PTM filtering", "zero_shift_filter"},
    {"One unexpected shift filtering", "zero_one_shift_filter"},
    {"Multiple unexpected shifts filtering", "diag_filter"},
    {"ASF-One PTM filtering", "zero_one_shift_filter"},
    {"ASF-Diagonal PTM filtering", "diag_filter"}};

const std::map<std::string, std::vector<std::string>> filter_stage_list{
    {"toppic",
     {"Zero unexpected shift filtering", "Variable PTM filtering",
      "One unexpected shift filtering",
      "Multiple unexpected shifts filtering"}},
    {"topmg", {"ASF-One PTM filtering", "ASF-Diagonal PTM filtering"}}};

namespace {

// the thread number a memory_per_thread_list entry allows with this much
// memory available, never zero
int threadNumForMem(const std::string& mem_entry, double avail_mem_in_gb) {
  double mem_per_thread = memory_per_thread_list.at(mem_entry);
  int thread_num = std::floor(avail_mem_in_gb / mem_per_thread);
  if (thread_num == 0) {
    thread_num = 1;
  }
  return thread_num;
}

// the memory_per_thread_list entry a filtering stage budgets with
const std::string& stageMemEntry(const std::string& stage_name) {
  auto it = filter_stage_mem_entry_list.find(stage_name);
  if (it == filter_stage_mem_entry_list.end()) {
    LOG_ERROR("invalid filtering stage name " << stage_name << "!");
    static const std::string default_entry = "zero_one_shift_filter";
    return default_entry;
  }
  return it->second;
}

// what checkThreadNum reports for a stage: its thread number computed from
// the memory available at the start of the run, which is the figure the same
// message prints
int getStageStartThreadNum(const std::string& stage_name, int thread_num) {
  int stage_thread_num =
      threadNumForMem(stageMemEntry(stage_name), getStartupAvailMemInGb());
  if (stage_thread_num > thread_num) {
    return thread_num;
  }
  return stage_thread_num;
}

}  // namespace

int getGraphAlignThreadNum(int thread_num, int layer_num) {
  double avail_mem_in_gb = getAvailMemInGb();
  double mem_per_thread =
      memory_per_thread_list.at("graph_align") *
      std::max(1.0, static_cast<double>(layer_num) / kGraphAlignLayers);
  int align_thread_num = std::floor(avail_mem_in_gb / mem_per_thread);
  if (align_thread_num < 1) {
    align_thread_num = 1;
  }
  if (align_thread_num >= thread_num) {
    return thread_num;
  }
  std::cout << "WARNING: Graph alignment will use " << align_thread_num
            << " threads, not the " << thread_num << " requested: "
            << std::setprecision(4) << avail_mem_in_gb
            << " GiB of memory is available and each thread may need "
            << mem_per_thread << " GiB." << std::endl;
  std::cout << std::setprecision(6);
  return align_thread_num;
}

double getTotalMemInGb() {
#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW32__) || \
    defined(__MINGW64__)
  MEMORYSTATUSEX mem_info;
  mem_info.dwLength = sizeof(MEMORYSTATUSEX);
  if (!GlobalMemoryStatusEx(&mem_info)) {
    return -1;
  }
  double bytes_per_gb = 1024.0 * 1024.0 * 1024.0;
  return mem_info.ullTotalPhys / bytes_per_gb;
#elif defined(__APPLE__)
  // hw.memsize is the physical memory size in bytes.
  int64_t total_mem = 0;
  size_t len = sizeof(total_mem);
  if (sysctlbyname("hw.memsize", &total_mem, &len, nullptr, 0) != 0) {
    return -1;
  }
  double bytes_per_gb = 1024.0 * 1024.0 * 1024.0;
  return total_mem / bytes_per_gb;
#else
  // Linux: MemTotal in /proc/meminfo is reported in kB.
  std::string token;
  std::ifstream file("/proc/meminfo");
  while (file >> token) {
    if (token == "MemTotal:") {
      double mem;
      if (file >> mem) {
        return mem / (1024.0 * 1024.0);
      } else {
        return -1;
      }
    }
    // Ignore the rest of the line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return 0;  // Nothing found
#endif
}

double getAvailMemInGb() {
  double avail_mem_in_gb = 0;
#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW32__) || \
    defined(__MINGW64__)
  // ullAvailPhys is the amount of physical memory currently available.
  MEMORYSTATUSEX mem_info;
  mem_info.dwLength = sizeof(MEMORYSTATUSEX);
  if (GlobalMemoryStatusEx(&mem_info)) {
    double bytes_per_gb = 1024.0 * 1024.0 * 1024.0;
    avail_mem_in_gb = mem_info.ullAvailPhys / bytes_per_gb;
  }
#elif defined(__APPLE__)
  // macOS has no MemAvailable; approximate it with the free and inactive
  // (reclaimable) pages, the closest analogue to Linux's MemAvailable.
  mach_port_t host_port = mach_host_self();
  vm_size_t page_size = 0;
  vm_statistics64_data_t vm_stat;
  mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
  if (host_page_size(host_port, &page_size) == KERN_SUCCESS &&
      host_statistics64(host_port, HOST_VM_INFO64,
                        reinterpret_cast<host_info64_t>(&vm_stat),
                        &count) == KERN_SUCCESS) {
    uint64_t avail_bytes =
        (static_cast<uint64_t>(vm_stat.free_count) + vm_stat.inactive_count) *
        page_size;
    double bytes_per_gb = 1024.0 * 1024.0 * 1024.0;
    avail_mem_in_gb = avail_bytes / bytes_per_gb;
  }
#else
  // Linux: MemAvailable in /proc/meminfo is reported in kB.
  std::string token;
  std::ifstream file("/proc/meminfo");
  while (file >> token) {
    if (token == "MemAvailable:") {
      double mem;
      if (file >> mem) {
        avail_mem_in_gb = mem / (1024.0 * 1024.0);
      }
      break;
    }
    // Ignore the rest of the line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
#endif
  if (avail_mem_in_gb < 0) {
    avail_mem_in_gb = 0;
  }
  return avail_mem_in_gb;
}

// return max thread number based on total memory size
int getMaxThreads(const std::string& app_name) {
  if (memory_per_thread_list.find(app_name) == memory_per_thread_list.end()) {
    LOG_ERROR("invalid application name!");
    return 0;
  }
  return threadNumForMem(app_name, getAvailMemInGb());
}

double getStartupAvailMemInGb() {
  // sampled at the first call, which is checkThreadNum's while it validates
  // the arguments, and returned unchanged afterwards: the figure the run
  // reports its thread numbers from stays the figure it printed
  static const double startup_avail_mem_in_gb = getAvailMemInGb();
  return startup_avail_mem_in_gb;
}

int getStageThreadNum(const std::string& stage_name, int thread_num) {
  // one sample of the available memory, so the number the stage uses and the
  // number the message below explains it with are the same measurement
  double avail_mem_in_gb = getAvailMemInGb();
  int stage_thread_num =
      threadNumForMem(stageMemEntry(stage_name), avail_mem_in_gb);
  if (stage_thread_num > thread_num) {
    stage_thread_num = thread_num;
  }
  // the drivers call this where the stage starts, so this is the number the
  // stage will really use; say so whenever the memory available now has moved
  // it away from the number checkThreadNum reported at the start of the run
  int start_thread_num = getStageStartThreadNum(stage_name, thread_num);
  if (stage_thread_num != start_thread_num) {
    std::cout << "WARNING: " << stage_name << " will use " << stage_thread_num
              << " threads, not the " << start_thread_num
              << " reported at the start of the run: "
              << (stage_thread_num < start_thread_num ? "less" : "more")
              << " memory is available now (" << std::setprecision(4)
              << avail_mem_in_gb << " GiB)." << std::endl;
    std::cout << std::setprecision(6);
  }
  return stage_thread_num;
}

std::vector<std::pair<std::string, int>> getLimitedFilterStages(
    const std::string& prog, int thread_num) {
  std::vector<std::pair<std::string, int>> limited_stages;
  if (filter_stage_list.find(prog) == filter_stage_list.end()) {
    return limited_stages;
  }
  for (const std::string& stage_name : filter_stage_list.at(prog)) {
    int stage_thread_num = getStageStartThreadNum(stage_name, thread_num);
    if (stage_thread_num < thread_num) {
      limited_stages.emplace_back(stage_name, stage_thread_num);
    }
  }
  return limited_stages;
}

bool checkThreadNum(int thread_number, const std::string& prog) {
  if (thread_number <= 0) {
    LOG_ERROR("Thread number " << thread_number
                               << " error! The value should be positive.");
    return false;
  }
  int total_thread_num = static_cast<int>(std::thread::hardware_concurrency());
  double total_mem_in_gb = mem_check::getTotalMemInGb();
  double avail_mem_in_gb = mem_check::getStartupAvailMemInGb();
  std::cout << "Total thread number: " << total_thread_num << std::endl;
  std::cout << "Total memory: " << std::setprecision(4) << total_mem_in_gb
            << " GiB" << std::endl;
  std::cout << "Available memory: " << avail_mem_in_gb << " GiB" << std::endl;
  std::cout << std::endl;
  // set precision to default
  std::cout << std::setprecision(6);

  if (thread_number > total_thread_num) {
    LOG_ERROR("Thread number "
              << thread_number << " error! The value is too large. At most "
              << total_thread_num << " threads are supported.");
    return false;
  }
  // in toppic and topmg, we automatically control thread numbers for
  // filtering, one budget per filtering stage: report the stages the
  // available memory limits and the number each one of them will use. Each
  // stage recomputes its number from the memory available when it starts,
  // through getStageThreadNum, which reports the difference if there is one.
  if (filter_stage_list.find(prog) != filter_stage_list.end()) {
    std::vector<std::pair<std::string, int>> limited_stages =
        getLimitedFilterStages(prog, thread_number);
    if (!limited_stages.empty()) {
      std::cout << "WARNING: With the memory available now, these protein"
                << " sequence filtering steps will use" << std::endl;
      std::cout << "WARNING: fewer than the " << thread_number
                << " threads requested:" << std::endl;
      for (const auto& stage : limited_stages) {
        std::cout << "WARNING:   " << stage.first << ": " << stage.second
                  << " threads" << std::endl;
      }
      std::cout << "WARNING: All other steps of proteoform identification will"
                << " use " << thread_number << " threads!" << std::endl;
      std::cout << std::endl;
    }
    return true;
  }
  int max_thread = mem_check::getMaxThreads(prog);
  if (max_thread < thread_number) {
    std::cout << "WARNING: Based on the available memory size, up to "
              << max_thread << " threads can be used!" << std::endl;
    std::cout << "WARNING: Please set the thread number to " << max_thread
              << " or the program may crash!" << std::endl;
    std::cout << std::endl;
  }
  return true;
}

}  // namespace mem_check

}  // namespace toppic
