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

#ifndef TOPPIC_TOPFD_ECSCORE_MS_MAP_MS_MAP_PEAK_HPP
#define TOPPIC_TOPFD_ECSCORE_MS_MAP_MS_MAP_PEAK_HPP

#include "ms/spec/peak.hpp"

namespace toppic {

class MsMapPeak : public Peak {
 public:
  MsMapPeak(int peak_id, int spec_id, PeakPtr peak);

  std::string getString();

  int getPeakId() const { return peak_id_; }

  void setPeakId(int peak_id) { peak_id_ = peak_id; }

  int getSpecId() const { return spec_id_; }

  void setSpecId(int spec_id) { spec_id_ = spec_id; }

  double getOriInte() const { return ori_inte_; }

  void setOriInte(double ori_inte) { ori_inte_ = ori_inte; }

  int getStartIdx() const { return start_idx_; }

  void setStartIdx(int start_idx) { start_idx_ = start_idx; }

  int getEndIdx() const { return end_idx_; }

  void setEndIdx(int end_idx) { end_idx_ = end_idx; }

  bool getNeighbor() const { return neighbor_; }

  void setNeighbor(bool neighbor) { neighbor_ = neighbor; }

 private:
  int peak_id_;
  int spec_id_;
  double ori_inte_;
  int start_idx_;
  int end_idx_;
  bool neighbor_;
};

typedef std::shared_ptr<MsMapPeak> MsMapPeakPtr;
typedef std::vector<MsMapPeakPtr> MsMapPeakPtrVec;
typedef std::vector<MsMapPeakPtrVec> MsMapPeakPtr2D;

}
#endif
