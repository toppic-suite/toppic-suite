//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SET_XIC_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SET_XIC_HPP

#include <memory>
#include <vector>

namespace toppic {

class Xic {
 public:
  Xic(std::vector<double> inte_ratio_list,
      std::vector<double> top_three_inte_list, 
      std::vector<double> all_peak_inte_list);

  std::vector<double> getInteRatioList() {return inte_ratio_list_; }

  double getInteRatio(int idx) {return inte_ratio_list_[idx]; }

  std::vector<double> getTopThreeInteList() { return top_three_inte_list_; }

  double getTopThreeInteSum() {return top_three_inte_sum_; }

  std::vector<double> getSmoothedInteList() { return smoothed_inte_list_; }

  std::vector<double> getAllPeakInteList() {return all_peak_inte_list_; }

  double getAllPeakInte(int idx) {return all_peak_inte_list_[idx]; }

  double getAllPeakInteSum() {return all_peak_inte_sum_;}

 private:
  std::vector<double> inte_ratio_list_;
  // the intensities are based on scaled theoretical envelopes
  std::vector<double> top_three_inte_list_;
  std::vector<double> all_peak_inte_list_;
  // smoothed top three intensity list
  std::vector<double> smoothed_inte_list_;
  double top_three_inte_sum_;
  double all_peak_inte_sum_;

  void moving_avg(int n);

  int moving_avg_win_ = 2;

};

typedef std::shared_ptr<Xic> XicPtr;
typedef std::vector<XicPtr> XicPtrVec;

}

#endif //TOPPIC_XIC_HPP
