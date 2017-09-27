// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_BASE_LOGGER_HPP_
#define PROT_BASE_LOGGER_HPP_

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
      std::cout << "LOG TRACE: " << __FILE__ << "[" << __LINE__ << "]: " << X << std::endl;     \
    }                                   \
  }

#define LOG_DEBUG(X)                    \
  {                                     \
    if (log_level <= LOG_LEVEL_DEBUG) { \
      std::cout << "LOG DEBUG: " << __FILE__ << "[" << __LINE__ << "]: " << X << std::endl;     \
    }                                   \
  }

#define LOG_INFO(X)                     \
  {                                     \
    if (log_level <= LOG_LEVEL_INFO) {  \
      std::cout << "LOG INFO: " << __FILE__ << "[" << __LINE__ << "]: " << X << std::endl;      \
    }                                   \
  }

#define LOG_WARN(X)                     \
  {                                     \
    if (log_level <= LOG_LEVEL_WARN) {  \
      std::cout << "LOG WARN: " << __FILE__ << "[" << __LINE__ << "]: " << X << std::endl;      \
    }                                   \
  }

#define LOG_ERROR(X)                    \
  {                                     \
    if (log_level <= LOG_LEVEL_ERROR) { \
      std::cout << "LOG ERROR: " << __FILE__ << "[" << __LINE__ << "]: " << X << std::endl;     \
    }                                   \
  }
}

#endif
