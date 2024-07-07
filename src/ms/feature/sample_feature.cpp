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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/feature/sample_feature.hpp"

namespace toppic {

SampleFeature::SampleFeature(const std::string &line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  sample_id_ = std::stoi(strs[0]);
  feat_id_ = std::stoi(strs[1]);
  mono_mass_ = std::stod(strs[2]);
  intensity_ = std::stod(strs[3]);
  time_begin_ = std::stod(strs[4]) * 60;
  time_end_ = std::stod(strs[5]) * 60;
  min_scan_ = std::stoi(strs[6]);
  max_scan_ = std::stoi(strs[7]);
  min_charge_ = std::stoi(strs[8]);
  max_charge_ = std::stoi(strs[9]);
  min_frac_id_ = std::stoi(strs[10]);
  max_frac_id_ = std::stoi(strs[11]);
  rep_frac_id_ = std::stoi(strs[12]);
  rep_apex_time_ = std::stod(strs[13]) * 60;
  rep_apex_scan_ = std::stoi(strs[14]);
  rep_apex_inte_ = std::stod(strs[15]);
  rep_charge_ = std::stoi(strs[16]);
  rep_avg_mz_ = std::stoi(strs[17]);
  env_num_ = std::stoi(strs[18]);
  rep_ec_score_ = std::stod(strs[19]);
}

SampleFeature::SampleFeature(FracFeaturePtr frac_feature, int feat_id) {
  feat_id_ = feat_id;
  init(frac_feature);
}

void SampleFeature::init(FracFeaturePtr frac_feature) {
  mono_mass_ = frac_feature->getMonoMass();
  intensity_ = frac_feature->getIntensity();
  time_begin_ = frac_feature->getTimeBegin();
  time_end_ = frac_feature->getTimeEnd();
  min_scan_ = frac_feature->getScanBegin();
  max_scan_ = frac_feature->getScanEnd();
  min_charge_ = frac_feature->getMinCharge();
  max_charge_ = frac_feature->getMaxCharge();
  min_frac_id_ = frac_feature->getFracId();
  max_frac_id_ = frac_feature->getFracId();
  rep_frac_id_ = frac_feature->getFracId();
  rep_apex_time_ = frac_feature->getApexTime();
  rep_apex_scan_ = frac_feature->getApexScan();
  rep_apex_inte_ = frac_feature->getApexInte();
  rep_charge_ = frac_feature->getRepCharge(); 
  rep_avg_mz_ = frac_feature->getRepAvgMz();
  env_num_ = frac_feature->getEnvNum(); 
  rep_ec_score_ = frac_feature->getEcScore(); 
}

SampleFeature::SampleFeature(FracFeaturePtrVec &frac_features, 
                             int feat_id) {
  if (frac_features.size() == 0) {
    LOG_ERROR("Fraction feature size is 0!");
    exit(EXIT_FAILURE);  
  }
  feat_id_ = feat_id;
  init(frac_features[0]);

  for (size_t i = 1; i < frac_features.size(); i++) {
    FracFeaturePtr cur_ft = frac_features[i];
    intensity_ = intensity_+ cur_ft->getIntensity();
    if (cur_ft->getTimeBegin() < time_begin_) {
      time_begin_ = cur_ft->getTimeBegin();
    }
    if (cur_ft->getTimeEnd() > time_end_) {
      time_end_ = cur_ft->getTimeEnd();
    }
    if (cur_ft->getScanBegin() < min_scan_) {
      min_scan_ = cur_ft->getScanBegin(); 
    }
    if (cur_ft->getScanEnd() > max_scan_) {
      max_scan_ = cur_ft->getScanEnd();
    }
    if (cur_ft->getMinCharge() < min_charge_) {
      min_charge_ = cur_ft->getMinCharge();
    }
    if (cur_ft->getMaxCharge() > max_charge_) {
      max_charge_ = cur_ft->getMaxCharge();
    }
    if (cur_ft->getFracId() < min_frac_id_) {
      min_frac_id_ = cur_ft->getFracId();
    }
    if (cur_ft->getFracId() > max_frac_id_) {
      max_frac_id_ = cur_ft->getFracId();
    }
    env_num_ = env_num_ + cur_ft->getEnvNum();
    if (cur_ft->getEcScore() > rep_ec_score_) {
      rep_ec_score_ = cur_ft->getEcScore();
    }
  }
}

}
