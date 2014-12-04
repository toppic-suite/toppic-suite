/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#ifndef PROT_LOGGER_HPP_
#define PROT_LOGGER_HPP_

#include <string>
#include <iostream>
#include <fstream>

namespace prot {

extern int log_level;

#define LOG_LEVEL_TRACE 1
#define LOG_LEVEL_DEBUG 2
#define LOG_LEVEL_INFO  3
#define LOG_LEVEL_WARN  4
#define LOG_LEVEL_ERROR 5

#define LOG_TRACE(X)                    \
  {                                     \
    if (log_level <= LOG_LEVEL_TRACE) { \
      std::cout << "LOG TRACE: " << __FILE__ << ": "  << X << std::endl;      \
    }                                   \
  }

#define LOG_DEBUG(X)                    \
  {                                     \
    if (log_level <= LOG_LEVEL_DEBUG) { \
      std::cout << "LOG DEBUG: " << __FILE__ << ": " << X << std::endl;      \
    }                                   \
  }

#define LOG_INFO(X)                     \
  {                                     \
    if (log_level <= LOG_LEVEL_INFO) {  \
      std::cout << "LOG INFO: " << __FILE__ << ": " << X << std::endl;      \
    }                                   \
  }

#define LOG_WARN(X)                     \
  {                                     \
    if (log_level <= LOG_LEVEL_WARN) {  \
      std::cout << "LOG WARN: " << __FILE__ << ": " << X << std::endl;      \
    }                                   \
  }

#define LOG_ERROR(X)                    \
  {                                     \
    if (log_level <= LOG_LEVEL_ERROR) { \
      std::cout << "LOG ERROR: " << __FILE__ << ": " << X << std::endl;      \
    }                                   \
  }

class WebLog {
 public:
  static void init(std::string log_file_name);
  static void percent_log(double p);
  static void close();

 private:
  static std::ofstream log_;
  static std::string log_file_name_;
};

}
#endif
