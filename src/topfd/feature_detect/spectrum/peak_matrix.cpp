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

#include "peak_matrix.hpp"
#include "topfd/feature_detect/test_output_functions/write_out_files.hpp"

namespace toppic {
  PeakMatrix::PeakMatrix() {
    bin_num_ = -1;
    spec_num_ = -1;
    min_mz_ = -1;
    max_mz_ = -1;
    min_inte_ = -1;
    bin_size_ = -1;
  };

  PeakMatrix::PeakMatrix(PeakPtrVec2D &raw_peaks, DeconvMsPtrVec &ms1_ptr_vec, double bin_size, double snr) {
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
    min_inte_ = getDataLevelNoiseIntensities(intes);
    spec_noise_inte_ = getSpectrumNoiseIntensities(raw_peaks);
    specs_ = get_spec_list(ms1_ptr_vec);
    bin_size_ = bin_size;
    spec_num_ = raw_peaks.size();
    bin_num_ = int((max_mz_ - min_mz_) / bin_size_) + 1;
    std::cout << "Data Level Noise Intensity Level: " << min_inte_ << std::endl;
    std::cout << "Min and Max m/z values: " << min_mz_ << " , " << max_mz_ << std::endl;
    init_matrix(raw_peaks, snr);
  }

  spec_list PeakMatrix::get_spec_list(DeconvMsPtrVec &ms1_ptr_vec) {
    std::vector<Spectrum> spec_list;
    for (auto &i: ms1_ptr_vec)
      spec_list.emplace_back(i->getMsHeaderPtr()->getId(), i->getMsHeaderPtr()->getFirstScanNum(),
                             i->getMsHeaderPtr()->getRetentionTime());
    return spec_list;
  }

  std::vector<double> PeakMatrix::getSpectrumNoiseIntensities(PeakPtrVec2D &raw_peaks) {
    std::vector<double> spectrum_noise_levels;
    for (auto &peaks: raw_peaks) {
      std::vector<double> intes;
      for (auto &peak: peaks) {
        intes.push_back(peak->getIntensity());
      }
      double noise = baseline_util::getBaseLine(intes);
      spectrum_noise_levels.push_back(noise);
    }
    return spectrum_noise_levels;
  }

  void PeakMatrix::init_matrix(PeakPtrVec2D &raw_peaks, double snr) {
    for (int spec_id = 0; spec_id < spec_num_; spec_id++) {
      matrix_[spec_id] = PeakRow(specs_[spec_id], bin_num_);
      std::vector<PeakPtr> spec_peak_list = raw_peaks[spec_id];
      std::vector<ExpPeak> exp_peak_list;
      for (auto &p: spec_peak_list) {
        int peak_id = peaks_.size();
        ExpPeak new_peak = ExpPeak(peak_id, spec_id, p);
        exp_peak_list.push_back(new_peak);
        peaks_.push_back(new_peak);
      }
      for (auto &cur_peak: exp_peak_list) {
        // Only keep peaks above data level noise intensity * SNR
        if (cur_peak.getInte() > snr * min_inte_) {
          int bin_idx = get_index(cur_peak.getPos());
          matrix_[spec_id].addPeak(bin_idx, cur_peak);
        }
      }
    }
  }

  int PeakMatrix::get_index(double mz) {
    double mz_diff = mz - min_mz_;
    int bin_idx = int(mz_diff / bin_size_);
    if (bin_idx < 0) bin_idx = 0;
    if (bin_idx >= bin_num_) bin_idx = bin_num_ - 1;
    return bin_idx;
  }

  void PeakMatrix::find_pair_neighbors(int spec_id, int search_bin_num, double mass_tol) {
    std::vector<std::vector<ExpPeak>> first_bin_list = matrix_[spec_id].getRow();
    std::vector<std::vector<ExpPeak>> second_bin_list = matrix_[spec_id + 1].getRow();
    for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++) {
      int start = std::max(0, bin_idx - search_bin_num);
      int end = std::min(bin_idx + search_bin_num, bin_num_ - 1);
      for (auto &first_peak: first_bin_list[bin_idx]) {
        for (int second_idx = start; second_idx < end + 1; second_idx++) {
          for (auto &second_peak: second_bin_list[second_idx]) {
            double mass_diff = std::abs(first_peak.getPos() - second_peak.getPos());
            if (mass_diff <= mass_tol) {
              first_peak.setNeighbor(true);
              second_peak.setNeighbor(true);
            }
          }
        }
      }
    }
    matrix_[spec_id].setRow(first_bin_list);
    matrix_[spec_id + 1].setRow(second_bin_list);
  }

  void PeakMatrix::find_remove_non_neighbors(double mass_tol) {
    int search_bin_num = int(mass_tol / bin_size_) + 1;
    for (auto &peak: peaks_)
      peak.setNeighbor(false);
    for (int spec_id = 0; spec_id < spec_num_ - 1; spec_id++)
      find_pair_neighbors(spec_id, search_bin_num, mass_tol);
    for (int spec_id = 0; spec_id < spec_num_ - 1; spec_id++) {
      std::vector<std::vector<ExpPeak>> bin_list = matrix_[spec_id].getRow();
      for (int bin_idx = 0; bin_idx < bin_num_; bin_idx++)
        for (auto &peak: bin_list[bin_idx])
          if (!peak.getNeighbor())
            peak = ExpPeak();
      matrix_[spec_id].setRow(bin_list);
    }
  }

  std::vector<toppic::ExpPeak> PeakMatrix::get_bin_peak(int sp_id, int peak_id) {
    return matrix_[sp_id].getRowPeak(peak_id);
  }

  void PeakMatrix::set_bin_peak(int sp_id, int peak_id, std::vector<ExpPeak> &bin_peaks) {
    matrix_[sp_id].setRowPeak(peak_id, bin_peaks);
  }
}