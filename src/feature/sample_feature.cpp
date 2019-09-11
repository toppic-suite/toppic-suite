//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
#include "feature/sample_feature.hpp"

namespace toppic {

SampleFeature::SampleFeature(const std::string &line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  sample_id_ = std::stoi(strs[0]);
  id_ = std::stoi(strs[1]);
  mono_mass_ = std::stod(strs[2]);
  intensity_ = std::stod(strs[3]);
  time_begin_ = std::stod(strs[4]);
  time_end_ = std::stod(strs[5]);
  min_charge_ = std::stoi(strs[6]);
  max_charge_ = std::stoi(strs[7]);
  min_frac_id_ = std::stoi(strs[8]);
  max_frac_id_ = std::stoi(strs[9]);
}

SampleFeature::SampleFeature(FracFeaturePtr frac_feature, int id) {
  id_ = id;
  FracFeaturePtr first_ft = frac_feature;
  mono_mass_ = first_ft->getMonoMass();
  intensity_ = first_ft->getIntensity();
  time_begin_ = first_ft->getTimeBegin();
  time_end_ = first_ft->getTimeEnd();
  min_charge_ = first_ft->getMinCharge();
  max_charge_ = first_ft->getMaxCharge();
  min_frac_id_ = first_ft->getFracId();
  max_frac_id_ = first_ft->getFracId();
}

SampleFeature::SampleFeature(FracFeaturePtrVec &frac_features, int id) {
  if (frac_features.size() == 0) {
    LOG_ERROR("Fraction feature size is 0!");
    exit(EXIT_FAILURE);  
  }
  id_ = id;
  FracFeaturePtr first_ft = frac_features[0];
  mono_mass_ = first_ft->getMonoMass();
  intensity_ = first_ft->getIntensity();
  time_begin_ = first_ft->getTimeBegin();
  time_end_ = first_ft->getTimeEnd();
  min_charge_ = first_ft->getMinCharge();
  max_charge_ = first_ft->getMaxCharge();
  min_frac_id_ = first_ft->getFracId();
  max_frac_id_ = first_ft->getFracId();
  for (size_t i = 1; i < frac_features.size(); i++) {
    FracFeaturePtr cur_ft = frac_features[i];
    intensity_ = intensity_+ cur_ft->getIntensity();
    if (cur_ft->getTimeBegin() < time_begin_) {
      time_begin_ = cur_ft->getTimeBegin();
    }
    if (cur_ft->getTimeEnd() > time_end_) {
      time_end_ = cur_ft->getTimeEnd();
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
  }
}

}
