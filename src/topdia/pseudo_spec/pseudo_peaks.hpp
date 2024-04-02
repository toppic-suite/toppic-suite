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

#ifndef TOPPIC_PSEUDO_PEAKS_HPP
#define TOPPIC_PSEUDO_PEAKS_HPP

namespace toppic {

    class PseudoPeaks {
    public:

        PseudoPeaks(double mass, double monoMz, int charge, double intensity, double score, double corr, double shared_inte,
                    int ms2_cycle_span, double apexDiffScan, double rtLow, double rtHigh, int ms2_apex_cycle);

        PseudoPeaks(const PseudoPeaks &peaks);

        double getMass() const { return mass_; }

        double getMonoMz() const { return mono_mz_; }

        int getCharge() const { return charge_; }

        double getIntensity() const { return intensity_; }

        double getScore() const { return score_; }

        double getCorr() const { return corr_; }

        double getSharedInte() const { return shared_inte_; }

        int getMS2CycleSpan() const { return ms2_cycle_span_; }

        double getApexDiffScan() const { return apex_diff_scan_; }

        double getRtLow() const { return rt_low_; }

        double getRtHigh() const { return rt_high_; }

        int getMs2FeatureIdx() const { return ms2_feature_idx_; }

        void setMs2FeatureIdx(int ms2FeatureIdx) { ms2_feature_idx_ = ms2FeatureIdx; }

        double getRank() const { return rank_; }

        void setRank(double rank) { rank_ = rank; }

        void setScore(double score) { score_ = score; }

        int getMs2ApexCycle() const { return ms2_apex_cycle_; }

        void setMs2ApexCycle(int ms2ApexSpec) { ms2_apex_cycle_ = ms2ApexSpec; }

    private:
        double mass_;
        double mono_mz_;
        int charge_;
        double intensity_;
        double score_;
        double corr_;
        double rank_;
        double shared_inte_;
        int ms2_cycle_span_;
        double apex_diff_scan_;
        double rt_low_;
        double rt_high_;
        int ms2_feature_idx_;
        int ms2_apex_cycle_;
    };
}

#endif //TOPPIC_PSEUDO_PEAKS_HPP
