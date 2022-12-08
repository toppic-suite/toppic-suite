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
#include <iomanip>
#include <fstream>
#include <limits>

#include "boost/thread/thread.hpp"
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
  // topfd memory requirement per thread: about 0.15 gb
  {"topfd", 0.3}, 
    // topindex memory requirement per thread: about 0.4 gb
    {"topindex", 0.5},  
    // toppic memory requirement per thread: about 0.75 gb
    {"toppic", 1.0},
    // topmg memory requirement per thread: about 0.75 gb
    {"topmg", 1.0}, 
    // zero or one shift filter memory requirement per thread: about 0.75 gb
    {"zero_one_shift_filter", 1.0},
    // diag filter memory requirement per thread: about 0.4 gb
    {"diag_filter", 0.5}
};


#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
bool isWindows11() {
  DWORD dwVersion = 0; 
  DWORD dwMajorVersion = 0;
  DWORD dwBuild = 0;

  dwVersion = GetVersion();

  // Get the Windows version.
  dwMajorVersion = (DWORD)(LOBYTE(LOWORD(dwVersion)));

  // Get the build number.
  if (dwVersion < 0x80000000)              
    dwBuild = (DWORD)(HIWORD(dwVersion));

  //std::cout << "Version is " << dwMajorVersion <<  " " << dwBuild << std::endl;
  if (dwMajorVersion >= 10 && dwBuild >= 22000) {
    return true;
  }
  else {
    return false;
  }  
}

#endif


double getTotalMemInGb () {
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  MEMORYSTATUSEX mem_info;
  mem_info.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&mem_info);
  DWORDLONG total_mem = mem_info.ullTotalPhys;
  double total_mem_in_gb = total_mem / pow(2, 30);
  return total_mem_in_gb;
#else
  std::string token;
  std::ifstream file("/proc/meminfo");
  while(file >> token) {
    if(token == "MemTotal:") {
      double mem;
      if(file >> mem) {
        double total_mem_in_gb = mem/1024/1024;
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

double getAvailMemInGb () {
  double avail_mem_in_gb = 0;
#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
  avail_mem_in_gb = getTotalMemInGb();
  if (isWindows11()) {
    // minus 3 for windows 11
    avail_mem_in_gb = avail_mem_in_gb - 3;
  }
  else {
    // minus 2 for Windows 10
    avail_mem_in_gb = avail_mem_in_gb - 2;
  }
#else
  std::string token;
  std::ifstream file("/proc/meminfo");
  while(file >> token) {
    if(token == "MemAvailable:") {
      double mem;
      if(file >> mem) {
        avail_mem_in_gb = mem/1024/1024;
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

int getMaxThreads(std::string app_name) {//return max thread number based on total memory size
  double avail_mem_in_gb = getAvailMemInGb(); 
  if (memory_per_thread_list.find(app_name) == memory_per_thread_list.end()) {
    LOG_ERROR("invalid application name!");
    return 0;
  }
  double mem_per_thread = memory_per_thread_list[app_name];
  int max_thread_num = std::floor(avail_mem_in_gb / mem_per_thread);
  //std::cout << "Available memory " << avail_mem_in_gb << " memory per thread " << mem_per_thread << " max_thread_num " << max_thread_num << std::endl;
  if (max_thread_num == 0) {
    max_thread_num = 1;
  }
  return max_thread_num;
}

bool checkThreadNum(int thread_number, std::string prog) {
  if (thread_number <= 0) {
    LOG_ERROR("Thread number " << thread_number << " error! The value should be positive.");
    return false;
  }
  int total_thread_num = static_cast<int>(boost::thread::hardware_concurrency());
  float total_mem_in_gb = mem_check::getTotalMemInGb(); 
  float avail_mem_in_gb = mem_check::getAvailMemInGb();
  std::cout << "Total thread number: " << total_thread_num << std::endl;
  std::cout << "Total memory: " << std::setprecision(4) << total_mem_in_gb << " GiB" << std::endl;
  std::cout << "Available memory: " << avail_mem_in_gb << " GiB" << std::endl;
  std::cout << std::endl;
  // set precision to default
  std::cout << std::setprecision(6);

  if(thread_number > total_thread_num){
    LOG_ERROR("Thread number " << thread_number << " error! The value is too large. At most " << total_thread_num << " threads are supported.");
    return false;
  }
  int max_thread = mem_check::getMaxThreads(prog);
  // std::cout << "Maximum thread number: " << max_thread << std::endl;
  if (max_thread < thread_number) {
    // in toppic, we automatically control thread numbers for filtering
    if (prog != "toppic" && prog != "topmg") {
      std::cout << "WARNING: Based on the available memory size, up to " << max_thread << " threads can be used!" << std::endl;
      std::cout << "WARNING: Please set the thread number to " << max_thread << " or the program may crash!" << std::endl;
      std::cout << std::endl;
    }
    else {
      std::cout << "WARNING: Based on the available memory size, " << max_thread << " threads will be used for protein sequence filtering and " << thread_number << " threads will be used for other steps in proteoform identification!"  << std::endl;
      std::cout << std::endl;
    }
  }
  return true;
}

}
}