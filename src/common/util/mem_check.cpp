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

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#else
double getAvailMemInGb () {
  std::string token;
  std::ifstream file("/proc/meminfo");
  while(file >> token) {
    if(token == "MemAvailable:") {
      double mem;
      if(file >> mem) {
        return mem/1024/1024;
      } else {
        return -1;
      }
    }
    // Ignore the rest of the line
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return 0; // Nothing found
}
#endif

int getMaxThreads(std::string app_name) {//return max thread number based on total memory size
  double freeMemInGb = -1;
  std::map<std::string, double> toppic_apps_memory_per_thread {
    {"topfd", 0.5}, 
      {"toppic", 2},
      {"toppic_filter", 2},
      {"topmg", 4}, 
      {"topmerge", 4}, 
      {"topdiff", 4},
      {"topindex", 1.5}  
  };
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  MEMORYSTATUSEX memInfo;
  memInfo.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&memInfo);
  DWORDLONG totalMem = memInfo.ullAvailPhys;
  freeMemInGb = totalMem / pow(2, 30);
#else
  //struct sysinfo si;
  //sysinfo(&si);
  //freeMemInGb = (int)(roundl(si.freeram / pow(10, 9)));
  freeMemInGb = getAvailMemInGb();
#endif
  if (freeMemInGb < 0) {
    LOG_ERROR("invalid memory size!");
    return 0;
  }
  if (toppic_apps_memory_per_thread.find(app_name) == toppic_apps_memory_per_thread.end()) {
    LOG_ERROR("invalid application name!");
    return 0;
  }
  LOG_DEBUG("Free ram " << freeMemInGb);
  int thread_num =  static_cast<int>(freeMemInGb / toppic_apps_memory_per_thread[app_name]);
  if (thread_num == 0) {
    thread_num = 1;
  }
  return thread_num;
}
}
}
