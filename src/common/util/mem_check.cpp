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
#include "common/util/mem_check.hpp"
#include "common/util/logger.hpp"

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#include <windows.h> 
#else
#include <sys/sysinfo.h>
#endif

namespace toppic {

namespace mem_check {

int getMaxThreads() {//return max thread number based on total memory size
int totalMemInGb = -1;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  MEMORYSTATUSEX memInfo;
  memInfo.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&memInfo);
  DWORDLONG totalMem = memInfo.ullTotalPhys;
  totalMemInGb = (int)(roundl(totalMem / pow(10, 9)));
#else
  struct sysinfo si;
  sysinfo(&si);
  totalMemInGb = (int)(roundl(si.totalram / pow(10, 9)));
#endif
  if (totalMemInGb < 0) {
    LOG_ERROR("invalid memory size!");
    return 0;
  }
  std::cout << "total memory: " << totalMemInGb << ", max thread: " << totalMemInGb / 4 << std::endl;
  return totalMemInGb / 4;
}
}
}