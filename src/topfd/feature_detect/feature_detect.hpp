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

#ifndef TOPPIC_TOPFD_FEATURE_DETECT_FEATURE_DETECT_HPP_
#define TOPPIC_TOPFD_FEATURE_DETECT_FEATURE_DETECT_HPP_

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include "topfd/feature_detect/feature_detect_old.hpp"
#include "common/util/file_util.hpp"
#include "seq/fasta_util.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/env/env_para.hpp"
#include "ms/feature/frac_feature.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "topfd/feature_detect/feature_para.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/env_collection/env_coll_util.hpp"
#include "topfd/feature_detect/feature/feature.hpp"
#include "topfd/feature_detect/util/write_feature.hpp"
#include "topfd/feature_detect/envelope/evaluate_envelope.hpp"
#include "topfd/feature_detect/test_output_functions/write_out_files.hpp"
#include "ms/feature/spec_feature.hpp"

namespace toppic {

namespace feature_detect {

void process(int frac_id, const std::string &sp_file_name,
             bool miss_level_one, const std::string &resource_dir, const std::string &activation, bool isFaims, const std::vector<std::pair<double, int>> voltage_vec);
};

}

#endif
