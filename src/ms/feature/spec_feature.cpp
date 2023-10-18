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

#include <fstream>

#include "common/util/str_util.hpp"
#include "ms/feature/spec_feature.hpp"

namespace toppic {

SpecFeature::SpecFeature(std::string line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  frac_id_ = std::stoi(strs[0]);
  file_name_ = strs[1];
  spec_id_ = std::stoi(strs[2]);
  scans_ = strs[3];
  ms_one_id_ = std::stoi(strs[4]);
  ms_one_scan_ = std::stoi(strs[5]);
  frac_feature_id_ = std::stoi(strs[6]);
  frac_feature_inte_ = std::stod(strs[7]);
  frac_feature_score_ = std::stod(strs[8]);
  frac_feature_time_apex_ = std::stod(strs[9]);
  sample_feature_id_ = std::stoi(strs[10]);
  sample_feature_inte_ = std::stod(strs[11]);
  prec_mono_mz_ = std::stod(strs[12]);
  prec_charge_ = std::stoi(strs[13]);
  prec_inte_ = std::stod(strs[14]);
}

SpecFeature::SpecFeature(MsHeaderPtr header, FracFeaturePtr feature,
                         double prec_mono_mz, int prec_charge, double prec_inte) {
  frac_id_ = header->getFractionId();
  file_name_ = header->getFileName();
  spec_id_ = header->getSpecId();
  scans_ = header->getScansString();
  ms_one_id_ = header->getMsOneId();
  ms_one_scan_ = header->getMsOneScan();
  frac_feature_id_ = feature->getId();
  frac_feature_inte_ = feature->getIntensity();
  frac_feature_score_ = feature->getEcscore();
  frac_feature_time_apex_ = feature->getApexTime();
  sample_feature_id_ = feature->getSampleFeatureId();
  sample_feature_inte_ = feature->getSampleFeatureInte();
  prec_mono_mz_ = prec_mono_mz;
  prec_charge_ = prec_charge; 
  prec_inte_ = prec_inte; 
}

}
