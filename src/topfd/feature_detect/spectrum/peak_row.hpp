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


#ifndef TOPPIC_PEAK_ROW_HPP
#define TOPPIC_PEAK_ROW_HPP

#include "spectrum.hpp"
#include "exp_peak.hpp"

namespace toppic {
  class PeakRow {
  public:
    PeakRow() = default;

    PeakRow(Spectrum spectrum, int bin_num) {
      spectrum_ = spectrum;
      for (int i = 0; i < bin_num; i++) {
        std::vector<ExpPeak> vec;
        row_.push_back(vec);
      }
    }

    Spectrum getSpectrum() const { return spectrum_; }

    void setSpectrum(Spectrum spectrum) { spectrum_ = spectrum; }

    std::vector<std::vector<ExpPeak>> getRow() const { return row_; }

    void setRow(std::vector<std::vector<ExpPeak>> row) {
      int num_rows = row.size();
      for (int i = 0; i < num_rows; i++)
        row_[i] = row[i];
    }

    std::vector<ExpPeak> getRowPeak(int peak_id) const { return row_[peak_id]; }

    void setRowPeak(int peak_id, std::vector<ExpPeak> &bin_peaks) { row_[peak_id] = bin_peaks; }

    int getSpecID() const { return spectrum_.getSpecId(); }

    int getScanNum() const { return spectrum_.getScanNum(); }

    double getRT() const { return spectrum_.getRt(); }

    void addPeak(int idx, ExpPeak &peak) { row_[idx].push_back(peak); }

  private:
    Spectrum spectrum_;
    std::vector<std::vector<ExpPeak>> row_;
  };

  typedef std::vector<Spectrum> spec_list;
}

#endif //TOPPIC_PEAK_ROW_HPP
