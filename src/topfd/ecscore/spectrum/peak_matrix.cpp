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

#include <algorithm>
#include <numeric>

#include "common/util/logger.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/ecscore/spectrum/peak_matrix.hpp"

namespace toppic {

PeakMatrix::PeakMatrix(PeakPtrVec2D &raw_peaks, DeconvMsPtrVec &ms1_ptr_vec, 
                       double bin_size, double sn_ratio) {
  bin_size_ = bin_size;
  /// get min mz value
  std::vector<double> intes;
  std::vector<double> mz;
  for (auto &raw_peak: raw_peaks) {
    for (auto &p: raw_peak) {
      mz.push_back(p->getPosition());
      intes.push_back(p->getIntensity());
    }
  }
  /// set values!
  min_mz_ = *std::min_element(mz.begin(), mz.end());
  max_mz_ = *std::max_element(mz.begin(), mz.end());
  base_inte_ = baseline_util::getBaseLine(intes);
  initSpectrumNoiseIntensities(raw_peaks);
  initSpecList(ms1_ptr_vec);
  bin_num_ = int((max_mz_ - min_mz_) / bin_size_) + 1;
  LOG_DEBUG("Data Level Noise Intensity Level: " << base_inte_);
  LOG_DEBUG("Min and Max m/z values: " << min_mz_ << " , " << max_mz_);
  initMatrix(raw_peaks, sn_ratio);
}

void PeakMatrix::initSpectrumNoiseIntensities(PeakPtrVec2D &raw_peaks) {
  for (auto &peaks: raw_peaks) {
    std::vector<double> intes;
    for (auto &peak: peaks) {
      intes.push_back(peak->getIntensity());
    }
    double noise = 0.0;
    if (std::accumulate(intes.begin(), intes.end(), 0.0) > 0)
      noise = baseline_util::getBaseLine(intes);
    spec_noise_intes_.push_back(noise);
  }
}

void PeakMatrix::initSpecList(DeconvMsPtrVec &ms1_ptr_vec) {
  for (auto &ms1: ms1_ptr_vec) {
    MsHeaderPtr header_ptr = ms1->getMsHeaderPtr();
    MsMapRowHeaderPtr spec_ptr = std::make_shared<MsMapRowHeader>(header_ptr->getSpecId(),
                                                                  header_ptr->getFirstScanNum(),
                                                                  header_ptr->getRetentionTime());
    spec_list_.push_back(spec_ptr);
  }
}

int PeakMatrix::getBinIndex(double mz) {
  double mz_diff = mz - min_mz_;
  int bin_idx = int(mz_diff / bin_size_);
  if (bin_idx < 0) bin_idx = 0;
  if (bin_idx >= bin_num_) bin_idx = bin_num_ - 1;
  return bin_idx;
}


void PeakMatrix::initMatrix(PeakPtrVec2D &raw_peaks, double sn_ratio) {
  for (size_t spec_id = 0; spec_id < raw_peaks.size(); spec_id++) {
    MsMapRowPtr row_ptr = std::make_shared<MsMapRow>(spec_list_[spec_id], bin_num_);
    matrix_.push_back(row_ptr);
  }
  int peak_id = 0;
  for (size_t i = 0; i < raw_peaks.size(); i++) {
    for (size_t j = 0; j < raw_peaks[i].size(); j++) {
      PeakPtr p_ptr = raw_peaks[i][j];
      MsMapPeakPtr new_peak_ptr = std::make_shared<MsMapPeak>(peak_id, i, p_ptr);
      all_peaks_.push_back(new_peak_ptr);
      // filter low intensity peak
      if (new_peak_ptr->getIntensity() > sn_ratio * base_inte_) {
        int bin_idx = getBinIndex(new_peak_ptr->getPosition());
        matrix_[i]->addPeak(bin_idx, new_peak_ptr);
      }
      peak_id++;
    }
  }
}

void PeakMatrix::findNeighbors(int spec_id, int search_bin_num, double mass_tol) {
  MsMapRowPtr first_row = matrix_[spec_id];
  MsMapRowPtr second_row = matrix_[spec_id + 1];
  for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++) {
    int start = std::max(0, bin_idx - search_bin_num);
    int end = std::min(bin_idx + search_bin_num, bin_num_ - 1);
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

void PeakMatrix::removeNonNeighbors(double mass_tol) {
  int search_bin_num = int(mass_tol / bin_size_) + 1;
  for (auto peak: all_peaks_) {
    peak->setNeighbor(false);
  }
  size_t spec_num = getSpecNum();
  for (size_t spec_id = 0; spec_id < spec_num - 1; spec_id++) {
    findNeighbors(spec_id, search_bin_num, mass_tol);
  }
  for (size_t spec_id = 0; spec_id < spec_num; spec_id++) {
    MsMapRowPtr peak_row_ptr = matrix_[spec_id];
    for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++) {
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
