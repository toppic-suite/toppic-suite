//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_FEATURE_DETECT_OLD_HPP
#define TOPPIC_FEATURE_DETECT_OLD_HPP


#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/time_util.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/mod_util.hpp"
#include "seq/fasta_util.hpp"
#include "seq/fasta_index_reader.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/match_env.hpp"
#include "ms/feature/frac_feature.hpp"
#include "ms/feature/single_charge_feature.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/peak_cluster.hpp"
#include "ms/feature/sample_feature.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "topfd/feature_detect/feature_para.hpp"

namespace toppic {

namespace feature_detect_old {
    void getSampleFeatures(SampleFeaturePtrVec &sample_features, FracFeaturePtrVec &frac_features, SpecFeaturePtrVec &spec_features);
    void readHeaders(const std::string & file_name, MsHeaderPtrVec &header_ptr_vec);
    void getMs2Features(DeconvMsPtrVec &ms1_ptr_vec, MsHeaderPtrVec &header_ptr_vec, FracFeaturePtrVec &features, FeatureParaPtr para_ptr, SpecFeaturePtrVec &ms2_features);
    void process(int frac_id, const std::string &sp_file_name,
                 bool miss_level_one, const std::string &resource_dir, const std::string &activation, bool isFaims, const std::vector<std::pair<double, int>> voltage_vec, double score_thr);
};

}

#endif //TOPPIC_FEATURE_DETECT_OLD_HPP
