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

#ifndef TOPPIC_GENERATE_PSEUDO_SPECTRUM_HPP
#define TOPPIC_GENERATE_PSEUDO_SPECTRUM_HPP

#include <cmath>
#include <set>
#include <memory>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include "ms/spec/msalign_reader_util.hpp"
#include "ms/env/env_para.hpp"
#include "topdia/common/topdia_para.hpp"
#include "topdia/pseudo_spec/mzrt_feature.hpp"
#include "topdia/pseudo_spec/pseudo_spectrum.hpp"

namespace toppic {
    class GeneratePseudoSpectrum {
    public:
        explicit GeneratePseudoSpectrum(const TopdiaParaPtr& para_ptr);
        void process(const TopdiaParaPtr& para_ptr);

        ///////////////////////
        static bool compareFeaturesInte(const mzrtFeaturePtr &a, const mzrtFeaturePtr &b) { return a->getIntensity() > b->getIntensity(); }

    private:
        double get_max_rt();
        static std::vector<double> interp(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &fp);
        mzrtFeaturePtrVec get_iso_win_ms1_features(int isolation_window_base_index, double window_size);
        static int get_apex_cycle_distance(const mzrtFeaturePtr &ms1_feature, const mzrtFeaturePtr &ms2_feature);
        static std::vector<double> moving_avg(std::vector<double> xic, int size);
        static double computeSharedArea(const std::vector<double> &xic1, const std::vector<double> &xic2);
        static void score_pseudo_peaks(std::vector<PseudoPeaks> &pseudo_peak_list, const mzrtFeaturePtr &ms1_feature);
        static double get_pred(double intensity_ratio, double shared_area, double length_ratio);
        static std::vector<PseudoPeaks> filterPseudoPeaks(const EnvParaPtr &env_para_ptr, const mzrtFeaturePtr &ms1_feature,
                                                  mzrtFeaturePtrVec &ms2_features_window, std::vector<PseudoPeaks> &pseudo_peak_list,
                                                  double cutoff, int min_peak_num);
        static void write_pseudo_spectrum(const TopdiaParaPtr &topdia_para_ptr, int ms1_feature_idx, const mzrtFeaturePtr &ms1_feature, std::vector<PseudoPeaks> &assigned_ms2_features);

        static bool comparePseudoPeaksInte(const PseudoPeaks &a, const PseudoPeaks &b) { return a.getIntensity() > b.getIntensity(); }
        static bool comparePseudoPeaksScore(const PseudoPeaks &a, const PseudoPeaks &b) { return a.getScore() > b.getScore(); }


        std::vector<double> isolation_window_base_mz_;
        std::vector<double> rt_ms1_;
        std::vector<std::vector<double>> rt_ms2_;
        mzrtFeaturePtrVec ms1_features_;
        std::vector<mzrtFeaturePtrVec> ms2_features_;
    };

    typedef std::shared_ptr <GeneratePseudoSpectrum> GeneratePseudoSpectrumPtr;
}
#endif //TOPPIC_GENERATE_PSEUDO_SPECTRUM_HPP
