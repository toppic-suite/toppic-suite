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

#ifndef TOPPIC_TOPFD_ECSCORE_MS_MAP_MS_MAP_HPP
#define TOPPIC_TOPFD_ECSCORE_MS_MAP_MS_MAP_HPP

#include "ms/spec/deconv_ms.hpp"
#include "topfd/ecscore/spectrum/ms_map_row.hpp"

namespace toppic {

class MsMap {
 public:

  MsMap(PeakPtrVec2D &raw_peaks, DeconvMsPtrVec &ms1_ptr_vec,
        double bin_size, double sn_ratio);

  int getColNum() { return col_num_; }

  int getRowNum() { return row_ptr_list_.size(); }

  double getBinSize() { return bin_size_;}

  double getMinMz() { return min_mz_; }

  double getMaxMz() { return max_mz_; }

  double getBaseInte() { return base_inte_; }

  MsMapRowPtr getRowPtr(int i) {return row_ptr_list_[i];}

  MsMapRowHeaderPtrVec getHeaderPtrList();

  void removeNonNeighbors(double mass_tol);

  int getBinIndex(double mz);

  MsMapPeakPtrVec getBinPeakList(int row_idx, int bin_idx) {
    return row_ptr_list_[row_idx]->getPeakPtrVec(bin_idx);}

  void setBinPeakList(int row_idx, int bin_idx, MsMapPeakPtrVec &peaks) {
    return row_ptr_list_[row_idx]->setPeakPtrVec(bin_idx, peaks);}

 private:
  void initMap(DeconvMsPtrVec &ms1_ptr_vec, double sn_ratio);

  void findNeighbors(int spec_id, int search_bin_num, double mass_tol);

  MsMapRowPtrVec row_ptr_list_;

  double bin_size_;
  int col_num_;
  double min_mz_;
  double max_mz_;
  double base_inte_;

  MsMapPeakPtrVec peaks_;
};

typedef std::shared_ptr<MsMap> MsMapPtr;

}

#endif
