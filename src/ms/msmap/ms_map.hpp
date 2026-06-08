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

#ifndef TOPPIC_MS_MSMAP_MS_MAP_HPP_
#define TOPPIC_MS_MSMAP_MS_MAP_HPP_

#include <memory>
#include <vector>

#include "ms/msmap/ms_map_row.hpp"
#include "ms/spec/deconv_ms.hpp"

namespace toppic {

class MsMap {
 public:
  MsMap(const PeakPtrVec2D &raw_peak_2d, const MsHeaderPtrVec &ms1_header_ptr_vec,
        double bin_size, double sn_ratio, bool single_scan_noise);

  int getColNum() const { return col_num_; }

  int getRowNum() const { return row_ptr_list_.size(); }

  double getBinSize() const { return bin_size_;}

  double getMinMz() const { return min_mz_; }

  double getMaxMz() const { return max_mz_; }

  double getBaseInte() const { return base_inte_; }

  const MsMapPeakPtr2D& get2DPeaks() const {return peaks_;}

  MsMapRowPtr getRowPtr(int i) const {return row_ptr_list_[i];}

  MsMapRowHeaderPtrVec getHeaderPtrList() const;

  std::vector<int> getScanListBySpecId(const std::vector<int> &spec_id_list) const;

  std::vector<double> getRtListBySpecId(const std::vector<int> &spec_id_list) const;

  void removeNonNeighbors(double mass_tol);

  int getColIndex(double mz) const;

  const MsMapPeakPtrVec& getBinPeakList(int row_idx, int bin_idx) const {
    return row_ptr_list_[row_idx]->getPeakPtrVec(bin_idx);}

  void setBinPeakList(int row_idx, int bin_idx, const MsMapPeakPtrVec &peaks) {
    return row_ptr_list_[row_idx]->setPeakPtrVec(bin_idx, peaks);}

  void reconstruct(double sn_ratio, bool single_scan_noise); 

 private:
  void initMap(const PeakPtrVec2D &raw_peak_2d, const MsHeaderPtrVec &ms1_header_ptr_vec,
               double sn_ratio, bool single_scan_noise);

  void findNeighbors(int spec_id, int search_bin_num, double mass_tol);

  MsMapRowPtrVec row_ptr_list_;

  double bin_size_;
  int col_num_;
  double min_mz_;
  double max_mz_;
  double base_inte_;

  MsMapPeakPtr2D peaks_;
};

using MsMapPtr = std::shared_ptr<MsMap>;

}  // namespace toppic

#endif
