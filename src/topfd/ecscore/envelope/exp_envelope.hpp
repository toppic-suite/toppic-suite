//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#ifndef TOPPIC_ECSCORE_ENVELOPE_EXP_ENVELOPE_HPP
#define TOPPIC_ECSCORE_ENVELOPE_EXP_ENVELOPE_HPP

#include <vector>

#include "ms/msmap/ms_map_peak.hpp"

namespace toppic {

class ExpEnvelope {
 public:
  ExpEnvelope(int spec_id, MsMapPeakPtrVec peak_list);

  int getMatchPeakNum(int base_idx);

  std::vector<double> getInteList();

  std::vector<double> getPosList();

  std::vector<double> getNonEmptyPosList();

  int getPeakNum() { return peak_list_.size(); }

  MsMapPeakPtr getPeakPtr(int idx) { return peak_list_[idx]; }

  int getSpecId() const { return spec_id_; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  MsMapPeakPtrVec getExpPeakList() { return peak_list_; }

  void setExpPeakList(const MsMapPeakPtrVec &peak_list) { peak_list_ = peak_list; }

  void setPeakPtr(int idx, MsMapPeakPtr peak_ptr) { peak_list_[idx] = peak_ptr; }

 private:
  int spec_id_;
  MsMapPeakPtrVec peak_list_;
};

typedef std::shared_ptr<ExpEnvelope> ExpEnvelopePtr;
typedef std::vector<ExpEnvelopePtr> ExpEnvelopePtrVec;

}

#endif //TOPPIC_EXP_ENVELOPE_HPP
