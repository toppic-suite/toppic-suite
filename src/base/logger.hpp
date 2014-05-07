/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#ifndef PROT_LOGGER_HPP_
#define PROT_LOGGER_HPP_

#include <string>
#include <iostream>


namespace prot {

#define LOG_LEVEL 4

#define LOG_LEVEL_TRACE 1
#define LOG_LEVEL_DEBUG 2
#define LOG_LEVEL_INFO  3
#define LOG_LEVEL_WARN  4
#define LOG_LEVEL_ERROR 5

#define LOG_TRACE(X)                    \
  {                                     \
    if (LOG_LEVEL <= LOG_LEVEL_TRACE) { \
      std::cout << "LOG TRACE: " << X << std::endl;      \
    }                                   \
  }

#define LOG_DEBUG(X)                    \
  {                                     \
    if (LOG_LEVEL <= LOG_LEVEL_DEBUG) { \
      std::cout << "LOG DEBUG: " << X << std::endl;      \
    }                                   \
  }

#define LOG_INFO(X)                     \
  {                                     \
    if (LOG_LEVEL <= LOG_LEVEL_INFO) {  \
      std::cout << "LOG INFO: " << X << std::endl;      \
    }                                   \
  }

#define LOG_WARN(X)                     \
  {                                     \
    if (LOG_LEVEL <= LOG_LEVEL_WARN) {  \
      std::cout << "LOG WARN: " << X << std::endl;      \
    }                                   \
  }

#define LOG_ERROR(X)                    \
  {                                     \
    if (LOG_LEVEL <= LOG_LEVEL_ERROR) { \
      std::cout << "LOG ERROR: " << X << std::endl;      \
    }                                   \
  }

}
#endif
