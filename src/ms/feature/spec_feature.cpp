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

#include <fstream>

#include "common/util/str_util.hpp"
#include "ms/feature/spec_feature.hpp"

namespace toppic {

SpecFeature::SpecFeature(std::string line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  file_name_ = strs[0];
  frac_id_ = std::stoi(strs[1]);
  spec_id_ = std::stoi(strs[2]);
  scans_ = strs[3];
  ms_one_id_ = std::stoi(strs[4]);
  ms_one_scan_ = std::stoi(strs[5]);
  frac_feature_id_ = std::stoi(strs[6]);
  frac_feature_inte_ = std::stod(strs[7]);
  frac_feature_score_ = std::stod(strs[8]);
  frac_feature_min_time_ = std::stod(strs[9]) * 60;
  frac_feature_max_time_ = std::stod(strs[10]) * 60;
  frac_feature_apex_time_ = std::stod(strs[11]) * 60;
  prec_mono_mz_ = std::stod(strs[13]);
  prec_avg_mz_ = std::stod(strs[14]);
  prec_charge_ = std::stoi(strs[15]);
  prec_inte_ = std::stod(strs[16]);
}

SpecFeature::SpecFeature(MsHeaderPtr header, FracFeaturePtr feature,
                         double prec_mono_mz, double prec_avg_mz, 
                         int prec_charge, double prec_inte) {
  file_name_ = header->getFileName();
  frac_id_ = feature->getFracId();
  spec_id_ = header->getSpecId();
  scans_ = header->getScansString();
  ms_one_id_ = header->getMsOneId();
  ms_one_scan_ = header->getMsOneScan();
  frac_feature_id_ = feature->getFeatId();
  frac_feature_inte_ = feature->getIntensity();
  frac_feature_score_ = feature->getEcScore();
  frac_feature_min_time_ = feature->getTimeBegin();
  frac_feature_max_time_ = feature->getTimeEnd();
  frac_feature_apex_time_ = feature->getApexTime();
  prec_mono_mz_ = prec_mono_mz;
  prec_avg_mz_ = prec_avg_mz;
  prec_charge_ = prec_charge; 
  prec_inte_ = prec_inte; 
}

SpecFeature::SpecFeature(double prec_mono_mz, 
                         double prec_charge) {
  frac_id_ = 0;
  frac_feature_id_ = -1;
  prec_mono_mz_ = prec_mono_mz;
  prec_charge_ = prec_charge; 
  prec_inte_ = 0; 
}

}
