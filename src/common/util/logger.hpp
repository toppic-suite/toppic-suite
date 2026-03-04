//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_UTIL_LOGGER_HPP_
#define TOPPIC_COMMON_UTIL_LOGGER_HPP_

#include <iostream>
#include <cstring>  // for strstr, used by SRC_FILENAME__ macro

namespace toppic {

namespace logger {

extern int log_level;
void setLogLevel(int level);

}  // namespace logger

inline constexpr int LOG_LEVEL_TRACE = 1;
inline constexpr int LOG_LEVEL_DEBUG = 2;
inline constexpr int LOG_LEVEL_INFO  = 3;
inline constexpr int LOG_LEVEL_WARN  = 4;
inline constexpr int LOG_LEVEL_ERROR = 5;

#if defined (_WIN32) || defined (_WIN64) || defined (__MINGW32__) || defined (__MINGW64__)
#define SRC_FILENAME__ (strstr(__FILE__, "\\src") ? strstr(__FILE__, "\\src") : __FILE__)
#else
#define SRC_FILENAME__ (strstr(__FILE__, "/src") ? strstr(__FILE__, "/src") + 1 : __FILE__)
#endif

#define LOG_TRACE(X)                    \
  {                                     \
    if (logger::log_level <= LOG_LEVEL_TRACE) { \
      std::cout << "LOG TRACE: " << SRC_FILENAME__ << "[" << __LINE__ << "]: " << X << std::endl;     \
    }                                   \
  }

#define LOG_DEBUG(X)                    \
  {                                     \
    if (logger::log_level <= LOG_LEVEL_DEBUG) { \
      std::cout << "LOG DEBUG: " << SRC_FILENAME__ << "[" << __LINE__ << "]: " << X << std::endl;     \
    }                                   \
  }

#define LOG_INFO(X)                     \
  {                                     \
    if (logger::log_level <= LOG_LEVEL_INFO) {  \
      std::cout << "LOG INFO: " << SRC_FILENAME__ << "[" << __LINE__ << "]: " << X << std::endl;      \
    }                                   \
  }

#define LOG_WARN(X)                     \
  {                                     \
    if (logger::log_level <= LOG_LEVEL_WARN) {  \
      std::cerr << "LOG WARN: " << SRC_FILENAME__ << "[" << __LINE__ << "]: " << X << std::endl;      \
    }                                   \
  }

#define LOG_ERROR(X)                    \
  {                                     \
    if (logger::log_level <= LOG_LEVEL_ERROR) { \
      std::cerr << "LOG ERROR: " << SRC_FILENAME__ << "[" << __LINE__ << "]: " << X << std::endl;     \
    }                                   \
  }
}  // namespace toppic

#endif
