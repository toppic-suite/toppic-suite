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

#ifndef TOPPIC_MZRT_FEATURE_HPP
#define TOPPIC_MZRT_FEATURE_HPP

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <memory>

namespace toppic {

    class mzrtFeature;
    typedef std::shared_ptr<mzrtFeature> mzrtFeaturePtr;
    typedef std::vector<mzrtFeaturePtr> mzrtFeaturePtrVec;

    class mzrtFeature {
    public:
        mzrtFeature(int id, int fraction_id, int env_num, double mass, double mono_mz, int charge, double intensity,
                    int mz_begin, int mz_end, double time_begin, double time_end, int spec_id_begin, int spec_id_end,
                    double time_apex, double ec_score, std::vector<double> xic, std::vector<double> normalized_xic,
                    std::vector<double> envelope_mz, std::vector<double> envelope_inte, int apex_cycle);
        static mzrtFeaturePtrVec read_record(std::string filename);

        int getId() const {return id_;}

        int getFractionId() const { return fraction_id_; }

        int getEnvNum() const { return env_num_; }

        double getMass() const { return mass_; }

        double getMonoMz() const { return mono_mz_; }

        int getCharge() const { return charge_; }

        double getIntensity() const { return intensity_; }

        int getMzBegin() const { return mz_begin_; }

        int getMzEnd() const { return mz_end_; }

        double getTimeBegin() const { return time_begin_; }

        double getTimeEnd() const { return time_end_; }

        int getSpecIDBegin() const { return spec_id_begin_; }

        int getSpecIDEnd() const { return spec_id_end_; }

        double getTimeApex() const { return time_apex_; }

        double getEcScore() const { return ec_score_; }

        const std::vector<double> &getXic() const { return xic_; }

        const std::vector<double> &getNormalizedXic() const { return normalized_xic_; }

        const std::vector<double> &getEnvelopeMz() const { return envelope_mz_; }

        const std::vector<double> &getEnvelopeInte() const { return envelope_inte_; }

        double getBaseMz() const { return base_mz_; }

        void setBaseMz(double base_mz) { base_mz_ = base_mz; }

        bool getUsedStatus() const { return used_; }

        void setUsedStatus(bool used) { used_ = used; }

        static std::vector<double> normalizeXIC(const std::vector<double> &xic);

        int getApexCycle() const { return apex_cycle_; }

        void setApexCycle(int apexCycle) { apex_cycle_ = apexCycle; }

        int getCycleSpan() const { return cycle_span_; }

        void setCycleSpan(int cycleSpan) { cycle_span_ = cycleSpan; }

        std::vector<double> getInterpolatedXic() { return interpolated_xic_; }

        void setInterpolatedXic(std::vector<double> interpolatedXic) { interpolated_xic_ = std::vector<double> (interpolatedXic.begin(), interpolatedXic.end()); }

        int getPseudoPeakNum() const { return pseudo_peak_num_; }

        void setPseudoPeakNum(int pseudoPeakNum) { pseudo_peak_num_ = pseudoPeakNum; }

    private:
        static std::vector<double> parseXIC(const std::string& line);
        static void parseEnvelope(const std::string& input, std::vector<double>& array1, std::vector<double>& array2);
        static int countNonZero(const std::vector<double> &xic);
        static std::vector<double> interp(const std::vector<double> &x, const std::vector<double> &xp, const std::vector<double> &fp);

        int id_;
        int fraction_id_;
        int env_num_;
        double mass_;
        double mono_mz_;
        int charge_;
        double intensity_;
        int mz_begin_;
        int mz_end_;
        double time_begin_;
        double time_end_;
        int spec_id_begin_;
        int spec_id_end_;
        double time_apex_;
        double ec_score_;
        std::vector<double> xic_;
        std::vector<double> interpolated_xic_;
        std::vector<double> normalized_xic_;
        std::vector<double> envelope_mz_;
        std::vector<double> envelope_inte_;

        int apex_cycle_;
        int cycle_span_;
        double base_mz_;
        bool used_;
        int pseudo_peak_num_ = 0;
    };
}

#endif //TOPPIC_MZRT_FEATURE_HPP
