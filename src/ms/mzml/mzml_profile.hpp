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

#ifndef TOPPIC_MS_MZML_MZML_PROFILE_HPP_
#define TOPPIC_MS_MZML_MZML_PROFILE_HPP_

#include <memory>
#include <vector>
#include <map>

namespace toppic {

class MzmlProfile {
 public:
  MzmlProfile(int ms1_cnt, int ms2_cnt, 
              std::map<double, std::pair<int,int>> volt_map);
  int getMs1Cnt() {return total_ms1_cnt_;}
  int getMs2Cnt() {return total_ms2_cnt_;}
  bool isFaims() {return volt_map_.size() > 0;}

  std::map<double, std::pair<int, int>> getVoltageMap() {return volt_map_;}

 private:
  int total_ms1_cnt_ = 0;
  int total_ms2_cnt_ = 0;
  std::map<double, std::pair<int,int>> volt_map_;
};

typedef std::shared_ptr<MzmlProfile> MzmlProfilePtr;

}

#endif
