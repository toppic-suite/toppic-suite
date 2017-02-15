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


#include "base/web_logger.hpp"

namespace prot {

std::ofstream WebLog::log_;
std::string WebLog::log_file_name_;

double WebLog::total_time_;
double WebLog::function_start_time_;

void WebLog::init(std::string log_file_name, bool use_gf, bool localization, int ptm_num) {
  log_file_name_ = log_file_name;
  if (log_file_name_.length() > 0) {
    log_.open(log_file_name, std::ios::out | std::ios::app);
  }
  function_start_time_ = 0;
  total_time_ = ZeroPtmFilterTime() + ZeroPtmSearchTime();
  if (ptm_num >= 1) {
    total_time_ += OnePtmFilterTime();
    total_time_ += OnePtmSearchTime();
  }
  if (ptm_num >= 2) {
    total_time_ += DiagFilterTime();
    total_time_ += TwoPtmSearchTime();
  }
  if (use_gf) {
    total_time_ += GfEvalueTime();  
  } else {
    total_time_ += TableEvalueTime();
  }
  if (localization) {
    total_time_ += LocalizationTime();
  }
  total_time_ += OutPutTime();
}

void WebLog::close() {
  if (log_.is_open()) {
    log_ << 1 << std::endl;
    log_.close();
  }
}

void WebLog::percentLog(int spectrum_index, int spectrum_num, int block_index, 
                        int block_num, double func_time) {
  double func_percentage = (double)block_index/block_num 
      + (double)spectrum_index/spectrum_num/block_num;
  double time = function_start_time_ + func_time * func_percentage;
  double total_percentage = time / total_time_;
  if (log_.is_open()) {
    log_ << total_percentage << std::endl;
  }
}

void WebLog::percentLog(int spectrum_index, int spectrum_num, double func_time) {
  return percentLog(spectrum_index, spectrum_num, 0, 1, func_time);
}

void WebLog::completeFunction(double func_time) {
  function_start_time_ += func_time;
}

}
