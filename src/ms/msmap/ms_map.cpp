//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include <algorithm>
#include <numeric>

#include "common/util/logger.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/msmap/ms_map.hpp"

namespace toppic {

MsMap::MsMap(PeakPtrVec2D &raw_peak_2d, DeconvMsPtrVec &ms1_ptr_vec,
             double bin_size, double sn_ratio, bool single_scan_noise) {
  bin_size_ = bin_size;
  /// get min mz value
  std::vector<double> intes;
  std::vector<double> mz;
  for (auto &row_peaks: raw_peak_2d) {
    for (auto &p: row_peaks) {
      mz.push_back(p->getPosition());
      intes.push_back(p->getIntensity());
    }
  }
  /// set values!
  min_mz_ = *std::min_element(mz.begin(), mz.end());
  max_mz_ = *std::max_element(mz.begin(), mz.end());
  base_inte_ = baseline_util::getBaseLine(intes);
  col_num_ = int((max_mz_ - min_mz_) / bin_size_) + 1;
  LOG_DEBUG("Data Level Noise Intensity Level: " << base_inte_);
  LOG_DEBUG("Min and Max m/z values: " << min_mz_ << " , " << max_mz_);

  initMap(raw_peak_2d, ms1_ptr_vec, sn_ratio, single_scan_noise);
}

void MsMap::initMap(PeakPtrVec2D &raw_peak_2d, DeconvMsPtrVec &ms1_ptr_vec, 
                    double sn_ratio, bool single_scan_noise) {
  double base_inte = base_inte_;
  for (size_t row_id = 0; row_id < ms1_ptr_vec.size(); row_id++) {
    MsHeaderPtr ms_header_ptr = ms1_ptr_vec[row_id]->getMsHeaderPtr();
    MsMapRowHeaderPtr row_header_ptr = std::make_shared<MsMapRowHeader>(ms_header_ptr->getSpecId(),
                                                                        ms_header_ptr->getFirstScanNum(),
                                                                        ms_header_ptr->getRetentionTime());
    PeakPtrVec row_peaks = raw_peak_2d[row_id]; 
    // get base peak intensity for each row
    std::vector<double> intes;
    for (auto &peak: row_peaks) { intes.push_back(peak->getIntensity()); }
    double row_base_inte = 0.0;
    if (std::accumulate(intes.begin(), intes.end(), 0.0) > 0) {
      row_base_inte = baseline_util::getBaseLine(intes);
    }
    row_header_ptr->setBaseInte(row_base_inte);
    if (single_scan_noise) {
      base_inte = row_base_inte;
    }
    // add a new row
    MsMapRowPtr row_ptr = std::make_shared<MsMapRow>(row_header_ptr, col_num_);
    row_ptr_list_.push_back(row_ptr);
    // init indexes
    MsMapPeakPtrVec ms_map_row_peaks;
    for (size_t j = 0; j < row_peaks.size(); j++) {
      PeakPtr p_ptr = row_peaks[j];
      MsMapPeakPtr new_peak_ptr = std::make_shared<MsMapPeak>(p_ptr);
      ms_map_row_peaks.push_back(new_peak_ptr);
      // filter low intensity peak
      if (new_peak_ptr->getIntensity() > sn_ratio * base_inte) {
        int bin_idx = getColIndex(new_peak_ptr->getPosition());
        row_ptr_list_[row_id]->addPeak(bin_idx, new_peak_ptr);
      }
    }
    peaks_.push_back(ms_map_row_peaks);
  }
}

void MsMap::reconstruct(double sn_ratio, bool single_scan_noise) {
  double base_inte = base_inte_;
  for (size_t row_id = 0; row_id < row_ptr_list_.size(); row_id++) {
    MsMapRowPtr row_ptr = row_ptr_list_[row_id];
    MsMapRowHeaderPtr row_header_ptr = row_ptr->getHeaderPtr(); 
    if (single_scan_noise) {
      base_inte = row_header_ptr->getBaseInte();
    }
    MsMapPeakPtrVec row_peaks = peaks_[row_id];
    row_ptr->clearPeaks();
    // init indexes
    for (size_t j = 0; j < row_peaks.size(); j++) {
      MsMapPeakPtr new_peak_ptr = row_peaks[j]; 
      // filter low intensity peak
      if (new_peak_ptr->getIntensity() > sn_ratio * base_inte) {
        int bin_idx = getColIndex(new_peak_ptr->getPosition());
        row_ptr->addPeak(bin_idx, new_peak_ptr);
      }
    }
  }
}

MsMapRowHeaderPtrVec MsMap::getHeaderPtrList() {
  MsMapRowHeaderPtrVec header_list;
  for (size_t i = 0; i < row_ptr_list_.size(); i++) {
    header_list.push_back(row_ptr_list_[i]->getHeaderPtr());
  }
  return header_list;
}

int MsMap::getColIndex(double mz) {
  double mz_diff = mz - min_mz_;
  int bin_idx = int(mz_diff / bin_size_);
  if (bin_idx < 0) bin_idx = 0;
  if (bin_idx >= col_num_) bin_idx = col_num_ - 1;
  return bin_idx;
}

void MsMap::findNeighbors(int row_id, int search_bin_num, double mass_tol) {
  MsMapRowPtr first_row = row_ptr_list_[row_id];
  MsMapRowPtr second_row = row_ptr_list_[row_id + 1];
  for (int bin_idx = 0; bin_idx < col_num_; bin_idx++) {
    int start = std::max(0, bin_idx - search_bin_num);
    int end = std::min(bin_idx + search_bin_num, col_num_ - 1);
    MsMapPeakPtrVec first_row_peaks = first_row->getPeakPtrVec(bin_idx);
    for (auto &first_peak: first_row_peaks) {
      for (int second_idx = start; second_idx < end + 1; second_idx++) {
        MsMapPeakPtrVec second_row_peaks = first_row->getPeakPtrVec(second_idx);
        for (auto &second_peak: second_row_peaks) {
          double mass_diff = std::abs(first_peak->getPosition() - second_peak->getPosition());
          if (mass_diff <= mass_tol) {
            first_peak->setNeighbor(true);
            second_peak->setNeighbor(true);
          }
        }
      }
    }
  }
}

void MsMap::removeNonNeighbors(double mass_tol) {
  int search_bin_num = int(mass_tol / bin_size_) + 1;
  for (size_t i = 0; i < peaks_.size(); i++) {
    MsMapPeakPtrVec ms_map_row_peaks = peaks_[i];
    for (auto peak: ms_map_row_peaks) {
      peak->setNeighbor(false);
    }
  }
  size_t row_num = getRowNum();
  for (size_t row_id = 0; row_id < row_num - 1; row_id++) {
    findNeighbors(row_id, search_bin_num, mass_tol);
  }
  for (size_t row_id = 0; row_id < row_num; row_id++) {
    MsMapRowPtr peak_row_ptr = row_ptr_list_[row_id];
    for (int bin_idx = 0; bin_idx < col_num_; bin_idx++) {
      MsMapPeakPtrVec peak_ptrs = peak_row_ptr->getPeakPtrVec(bin_idx);
      MsMapPeakPtrVec new_peak_ptrs;
      for (auto peak: peak_ptrs) {
        if (peak->getNeighbor()) {
          new_peak_ptrs.push_back(peak);
        }
      }
      peak_row_ptr->setPeakPtrVec(bin_idx, new_peak_ptrs);
    }
  }
}

}
