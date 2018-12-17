//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_COMP_SHIFT_LOW_MEM_HPP_
#define PROT_COMP_SHIFT_LOW_MEM_HPP_

#include <memory>
#include <vector>

namespace toppic {

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
} /* namespace toppic */

#endif /* COMP_SHIFT_LOW_MEM_HPP_ */
