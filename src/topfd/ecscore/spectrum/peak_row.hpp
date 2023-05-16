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

#ifndef TOPPIC_TOPFD_ECSCORE_SPECTRUM_PEAK_ROW_HPP
#define TOPPIC_TOPFD_ECSCORE_SPECTRUM_PEAK_ROW_HPP

#include "ms/spec/baseline_util.hpp"
#include "topfd/ecscore/spectrum/matrix_peak.hpp"
#include "topfd/ecscore/spectrum/matrix_spectrum.hpp"

namespace toppic {

class PeakRow {
 public:
  PeakRow(MatrixSpectrumPtr spec_ptr, int bin_num);

  MatrixSpectrumPtr getSpectrumPtr() { return spectrum_ptr_; }

  void setSpectrumPtr(MatrixSpectrumPtr spectrum_ptr) { spectrum_ptr = spectrum_ptr; }

  int getSpecID() const { return spectrum_ptr_->getSpecId(); }

  int getScanNum() const { return spectrum_ptr_->getScanNum(); }

  double getRT() const { return spectrum_ptr_->getRt(); }

  void addPeak(int idx, MatrixPeakPtr peak_ptr) { peak_ptrs_[idx].push_back(peak_ptr);}

  MatrixPeakPtrVec getPeakPtrVec(int i) {return peak_ptrs_[i];}

  void setPeakPtrVec(int i, MatrixPeakPtrVec &peak_ptr_vec) {peak_ptrs_[i] = peak_ptr_vec;}

 private:
  MatrixSpectrumPtr spectrum_ptr_;
  MatrixPeakPtr2D peak_ptrs_;
};

typedef std::shared_ptr<PeakRow> PeakRowPtr;
typedef std::vector<PeakRowPtr> PeakRowPtrVec;

}

#endif //TOPPIC_PEAK_ROW_HPP
