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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <thread>

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
    // diag filter memory requirement per thread: about 0.4 gb
    {"diag_filter", 0.5}};

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
  double avail_mem_in_gb = getAvailMemInGb();
  double mem_per_thread = memory_per_thread_list.at(app_name);
  int max_thread_num = std::floor(avail_mem_in_gb / mem_per_thread);
  if (max_thread_num == 0) {
    max_thread_num = 1;
  }
  return max_thread_num;
}

bool checkThreadNum(int thread_number, const std::string& prog) {
  if (thread_number <= 0) {
    LOG_ERROR("Thread number " << thread_number
                               << " error! The value should be positive.");
    return false;
  }
  int total_thread_num = static_cast<int>(std::thread::hardware_concurrency());
  double total_mem_in_gb = mem_check::getTotalMemInGb();
  double avail_mem_in_gb = mem_check::getAvailMemInGb();
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
  int max_thread = mem_check::getMaxThreads(prog);
  if (max_thread < thread_number) {
    // in toppic, we automatically control thread numbers for filtering
    if (prog != "toppic" && prog != "topmg") {
      std::cout << "WARNING: Based on the available memory size, up to "
                << max_thread << " threads can be used!" << std::endl;
      std::cout << "WARNING: Please set the thread number to " << max_thread
                << " or the program may crash!" << std::endl;
      std::cout << std::endl;
    } else {
      std::cout << "WARNING: Based on the available memory size, " << max_thread
                << " threads will be used for protein sequence filtering and "
                << thread_number
                << " threads will be used for other steps in proteoform "
                   "identification!"
                << std::endl;
      std::cout << std::endl;
    }
  }
  return true;
}

}  // namespace mem_check

}  // namespace toppic
