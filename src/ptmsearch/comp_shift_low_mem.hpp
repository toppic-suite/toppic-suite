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


#ifndef PROT_COMP_SHIFT_LOW_MEM_HPP_
#define PROT_COMP_SHIFT_LOW_MEM_HPP_

#include <memory>
#include <vector>

namespace prot {

class CompShiftLowMem {
 public:

  std::vector<std::vector<int>> findBestShift(const std::vector<int> &a,
                                              const std::vector<int> &b);

  std::vector<double> findBestShift(const std::vector<std::pair<int,int>> &a_e,
                                    const std::vector<int> &b,
                                    int total,int minimum_gap, double scale);

  /*
  std::vector<double> findBestShift(const std::vector<int> &a,
                                    const std::vector<int> &errors,
                                    const std::vector<int> &b,
                                    int total,int minimum_gap,
                                    double scale);
                                    */

  std::vector<std::vector<int>> findBestShift(const std::vector<int> &a,
                                              const std::vector<int> &errors,
                                              const std::vector<int> &b,
                                              int total,int minimum_gap) ;

 private:
  std::vector<short> num_;
  std::vector<std::vector<int>> findBestShift(const std::vector<int> &a,
                                              const std::vector<int> &b,
                                              int total,int minimum_gap) ;
  int checkD(std::vector<std::vector<int>> &ans,int d,int cur_min,
             int total,int min_gap);

  void resetNumbers(const std::vector<int> &a, const std::vector<int> &errors,
                    const std::vector<int> &b);
  void resetNumbers(const std::vector<int> &a, const std::vector<int> &b);
  void resetNumbers();
};

typedef std::shared_ptr<CompShiftLowMem> CompShiftLowMemPtr;
} /* namespace prot */

#endif /* COMP_SHIFT_LOW_MEM_HPP_ */
