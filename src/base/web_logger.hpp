// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_BASE_WEB_LOGGER_HPP_
#define PROT_BASE_WEB_LOGGER_HPP_

#include <string>
#include <iostream>
#include <fstream>

namespace prot {

class WebLog {
 public:
  static double ZeroPtmFilterTime() {return 1.0;}
  static double ZeroPtmSearchTime() {return 1.0;}
  static double OnePtmFilterTime() {return 1.0;}
  static double OnePtmSearchTime() {return 2.0;}
  static double DiagFilterTime() {return 6.0;}
  static double TwoPtmSearchTime() {return 6.0;}
  static double TableEvalueTime() {return 2.0;}
  static double GfEvalueTime() {return 76.0;}
  static double SelectingTime() {return 1.0;}
  static double LocalizationTime() {return 3.0;}
  static double OutPutTime() {return 1.0;}

  static void init(std::string log_file_name, bool use_gf, bool localization, int ptm_num);
  static void close();

  static void percentLog(int spectrum_index, int spectrum_num,
                         int block_index, int block_num, double func_time);

  static void percentLog(int spectrum_index, int spectrum_num, double func_time);

  static void completeFunction(double func_time);

 private:
  static std::string log_file_name_;
  static std::ofstream log_;

  static double total_time_;
  static double function_start_time_;
};

}
#endif
