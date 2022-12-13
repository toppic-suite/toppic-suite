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

#ifndef TOPPIC_PEAK_MATRIX_HPP
#define TOPPIC_PEAK_MATRIX_HPP

#include <map>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "spectrum.hpp"
#include "exp_peak.hpp"
#include "peak_row.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/peak.hpp"

namespace toppic {
  class PeakMatrix {
  public:
    PeakMatrix();

    PeakMatrix(PeakPtrVec2D &raw_peaks, DeconvMsPtrVec &ms1_ptr_vec, double bin_size, double snr);

    std::vector<double> getSpectrumNoiseIntensities(PeakPtrVec2D &raw_peaks);

    spec_list get_spec_list(DeconvMsPtrVec &ms1_ptr_vec);

    void init_matrix(PeakPtrVec2D &raw_peaks, double snr);

    int get_index(double mz);

    void find_remove_non_neighbors(double mass_tol);

    void find_pair_neighbors(int spec_id, int search_bin_num, double mass_tol);

    PeakRow get_row(int idx) { return matrix_[idx]; }

    std::vector<ExpPeak> get_bin_peak(int sp_id, int peak_id);

    void set_bin_peak(int sp_id, int peak_id, std::vector<ExpPeak> &bin_peaks);

    std::vector<std::vector<ExpPeak>> getRow(int idx) { return matrix_[idx].getRow(); }

    double getDataLevelNoiseIntensities(std::vector<double> intes) { return baseline_util::getBaseLine(intes); }

    int get_bin_num() { return bin_num_; }

    int get_spec_num() { return spec_num_; }

    double get_min_mz() { return min_mz_; }

    double get_max_mz() { return max_mz_; }

    double get_min_inte() { return min_inte_; }

    spec_list get_spectra_list() { return specs_; }

    std::vector<double> get_spec_noise_inte() { return spec_noise_inte_; }

  private:
    std::map<int, PeakRow> matrix_;
    std::vector<ExpPeak> peaks_;
    std::vector<double> spec_noise_inte_;
    spec_list specs_;
    int bin_num_;
    int spec_num_;
    double min_mz_;
    double max_mz_;
    double min_inte_;
    double bin_size_;
  };
}


#endif //TOPPIC_PEAK_MATRIX_HPP
