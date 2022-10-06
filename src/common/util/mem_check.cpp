//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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

#include <math.h>
#include <iostream>
#include <fstream>
#include <limits>
#include "common/util/mem_check.hpp"
#include "common/util/logger.hpp"

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#include <windows.h> 
#else
#include <sys/sysinfo.h>
#endif

namespace toppic {

namespace mem_check {

std::map<std::string, double> memory_per_thread_list {
  {"topfd", 0.5}, 
    {"toppic", 2},
    {"toppic_filter", 2},
    {"topmg", 4}, 
    {"topmerge", 4}, 
    {"topdiff", 4},
    {"topindex", 1.5}  
};


double getTotalMemInGb () {
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  MEMORYSTATUSEX mem_info;
  mem_info.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&mem_info);
  DWORDLONG total_mem = mem_info.ullAvailPhys;
  free_mem_in_gb = total_mem / pow(2, 30);
  return free_mem_in_gb;
#else
  std::string token;
  std::ifstream file("/proc/meminfo");
  while(file >> token) {
    if(token == "MemTotal:") {
      double mem;
      if(file >> mem) {
        double total_mem_in_gb = mem/1000/1000;
        return total_mem_in_gb;
      } else {
        return -1;
      }
    }
    // Ignore the rest of the line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return 0; // Nothing found
#endif
}

int getMaxThreads(std::string app_name) {//return max thread number based on total memory size
  // total memory - 1
  double avail_mem_in_gb = getTotalMemInGb() - 1; 

  if (avail_mem_in_gb < 0) {
    LOG_ERROR("invalid memory size!");
    return 0;
  }
  if (memory_per_thread_list.find(app_name) == memory_per_thread_list.end()) {
    LOG_ERROR("invalid application name!");
    return 0;
  }
  double mem_per_thread = memory_per_thread_list[app_name];
  int max_thread_num =  static_cast<int>(avail_mem_in_gb / mem_per_thread);
  //std::cout << "Available memory " << avail_mem_in_gb << " memory per thread " << mem_per_thread << " max_thread_num " << max_thread_num << std::endl;
  if (max_thread_num == 0) {
    max_thread_num = 1;
  }
  return max_thread_num;
}

}
}
