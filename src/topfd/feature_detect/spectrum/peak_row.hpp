//
// Created by abbash on 8/22/22.
//

#ifndef TOPPIC_PEAK_ROW_HPP
#define TOPPIC_PEAK_ROW_HPP

#include "spectrum.hpp"
#include "exp_peak.hpp"

namespace toppic {
    class PeakRow {
    public:
        PeakRow(){}
        PeakRow(Spectrum spectrum, int bin_num){
          spectrum_ = spectrum;
          for (int i = 0; i < bin_num; i++) {
            std::vector<ExpPeak> vec;
            row_.push_back(vec);
          }
        }

        Spectrum getSpectrum() const { return spectrum_; }
        void setSpectrum(Spectrum spectrum) { spectrum_ = spectrum; }

        std::vector<std::vector<ExpPeak>> getRow() const { return row_; }
        void setRow(std::vector<std::vector<ExpPeak>> row) { row_ = row; }

        int getSpecID() const { return spectrum_.getSpecId(); }
        int getScanNum() const { return spectrum_.getScanNum(); }
        double getRT() const { return spectrum_.getRt(); }

        void addPeak(int idx, ExpPeak peak) { row_[idx].push_back(peak); }

    private:
        Spectrum spectrum_;
        std::vector<std::vector<ExpPeak>> row_;
    };
    typedef std::vector<Spectrum> spec_list;
}

#endif //TOPPIC_PEAK_ROW_HPP
