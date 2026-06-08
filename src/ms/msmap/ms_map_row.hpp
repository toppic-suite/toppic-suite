//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_MS_MSMAP_MS_MAP_ROW_HPP_
#define TOPPIC_MS_MSMAP_MS_MAP_ROW_HPP_

#include <memory>
#include <vector>

#include "ms/msmap/ms_map_peak.hpp"
#include "ms/msmap/ms_map_row_header.hpp"

namespace toppic {

class MsMapRow {
 public:
  MsMapRow(const MsMapRowHeaderPtr &spec_ptr, int bin_num);

  MsMapRowHeaderPtr getHeaderPtr() const { return header_ptr_; }

  void setHeaderPtr(const MsMapRowHeaderPtr &header_ptr) { header_ptr_ = header_ptr; }

  int getSpecID() const { return header_ptr_->getSpecId(); }

  int getScanNum() const { return header_ptr_->getScanNum(); }

  double getRT() const { return header_ptr_->getRt(); }

  void addPeak(int idx, const MsMapPeakPtr &peak_ptr) { peak_ptr_2d_[idx].push_back(peak_ptr);}

  const MsMapPeakPtrVec& getPeakPtrVec(int i) const {return peak_ptr_2d_[i];}

  void setPeakPtrVec(int i, const MsMapPeakPtrVec &peak_ptr_vec) { peak_ptr_2d_[i] = peak_ptr_vec;}

  void clearPeaks();

  void print() const;

 private:
  MsMapRowHeaderPtr header_ptr_;
  MsMapPeakPtr2D peak_ptr_2d_;
};

using MsMapRowPtr = std::shared_ptr<MsMapRow>;
using MsMapRowPtrVec = std::vector<MsMapRowPtr>;

}  // namespace toppic

#endif
