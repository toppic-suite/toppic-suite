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


#ifndef TOPPIC_PSEUDO_SPECTRUM_HPP
#define TOPPIC_PSEUDO_SPECTRUM_HPP

#include <vector>
#include "pseudo_peaks.hpp"
#include "mzrt_feature.hpp"

namespace toppic {

    class PseudoSpectrum {
    public:
        PseudoSpectrum(mzrtFeaturePtr ms1Feature, const std::vector<PseudoPeaks> &fragmentFeatures) : ms1_feature_(
                ms1Feature), fragment_features_(fragmentFeatures), assigned_fragment_count_(fragmentFeatures.size()) {}

        void add_fragments(const std::vector<PseudoPeaks> &fragmentFeatures);
        void add_fragment(PseudoPeaks &fragmentFeature);
        mzrtFeaturePtr getMs1Feature() const { return ms1_feature_;  }
        std::vector<PseudoPeaks>& getFragmentFeatures() { return fragment_features_; }
        int getAssignedFragmentCount() const { return assigned_fragment_count_;  }

    private:
        mzrtFeaturePtr ms1_feature_;
        std::vector<PseudoPeaks> fragment_features_;
        int assigned_fragment_count_;

    };
    typedef std::shared_ptr<PseudoSpectrum> PseudoSpectrumPtr;
    typedef std::vector<PseudoSpectrumPtr> PseudoSpectrumPtrVec;

} // toppic

#endif //TOPPIC_PSEUDO_SPECTRUM_HPP
