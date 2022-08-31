//
// Created by abbash on 8/23/22.
//

#ifndef TOPPIC_PEAK_MATRIX_HPP
#define TOPPIC_PEAK_MATRIX_HPP

#include <map>
#include "spectrum.hpp"
#include "exp_peak.hpp"
#include "peak_row.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/spec/deconv_ms.hpp"

namespace toppic {
    class PeakMatrix {
    public:
        PeakMatrix(const PeakPtrVec2D& raw_peaks, DeconvMsPtrVec ms1_ptr_vec, double bin_size, double snr);
        spec_list get_spec_list(DeconvMsPtrVec ms1_ptr_vec);
        void init_matrix(PeakPtrVec2D raw_peaks, double snr);
        int get_index(double mz);
        void find_all_neighbors(double mass_tol);
        void find_pair_neighbors(PeakRow first_row, PeakRow second_row, int search_bin_num, double mass_tol);
        void remove_peak(ExpPeak peak);
        void remove_peak_in_range(int spec_id, double min_pos, double max_pos);
        PeakRow get_row(int idx) { return matrix_[idx]; }
        double getDataLevelNoiseIntensities(std::vector<double> intes){ return baseline_util::getBaseLine(intes); }

        int get_bin_num(){ return bin_num_; }
        int get_spec_num(){ return spec_num_; }
        double get_min_mz(){ return min_mz_; }
        double get_max_mz(){ return max_mz_; }
        double get_min_inte(){ return min_inte_; }
        spec_list get_spectra_list(){ return specs_; }

    private:
        std::map<int, PeakRow> matrix_;
        std::vector<ExpPeak> peaks_;
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
