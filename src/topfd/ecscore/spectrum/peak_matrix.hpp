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

#ifndef TOPPIC_TOPFD_ECSCORE_SPECTRUM_PEAK_MATRIX_HPP
#define TOPPIC_TOPFD_ECSCORE_SPECTRUM_PEAK_MATRIX_HPP


#include "ms/spec/peak.hpp"
#include "ms/spec/deconv_ms.hpp"

#include "topfd/ecscore/spectrum/matrix_spectrum.hpp"
#include "topfd/ecscore/spectrum/peak_row.hpp"

namespace toppic {

class PeakMatrix {
 public:

  PeakMatrix(PeakPtrVec2D &raw_peaks, DeconvMsPtrVec &ms1_ptr_vec, 
             double bin_size, double sn_ratio);

  int getBinNum() { return bin_num_; }

  int getSpecNum() { return spec_list_.size(); }

  double getMinMz() { return min_mz_; }

  double getMaxMz() { return max_mz_; }

  double getBaseInte() { return base_inte_; }

  MatrixSpectrumPtrVec getSpecList() { return spec_list_; }

  std::vector<double> getSpecNoiseInte() { return spec_noise_intes_; }

  void removeNonNeighbors(double mass_tol);

 private:
  void initSpectrumNoiseIntensities(PeakPtrVec2D &raw_peaks);

  void initSpecList(DeconvMsPtrVec &ms1_ptr_vec);

  void initMatrix(PeakPtrVec2D &raw_peaks, double sn_ratio);

  int getBinIndex(double mz);

  void findNeighbors(int spec_id, int search_bin_num, double mass_tol);

  PeakRowPtrVec matrix_;
  MatrixSpectrumPtrVec spec_list_;
  std::vector<double> spec_noise_intes_;
  double bin_size_;
  int bin_num_;
  double min_mz_;
  double max_mz_;
  double base_inte_;

  MatrixPeakPtrVec all_peaks_;
};
}


#endif //TOPPIC_PEAK_MATRIX_HPP
