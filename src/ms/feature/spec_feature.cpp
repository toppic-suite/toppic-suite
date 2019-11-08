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

#include <fstream>

#include "common/util/str_util.hpp"
#include "ms/feature/spec_feature.hpp"

namespace toppic {

/*
SpecFeature::SpecFeature(int spec_id, int frac_id, 
                         const std::string &file_name,
                         std::string &scans,
                         int ms_one_id, int ms_one_scan, 
                         double prec_mass, double prec_inte,
                         int frac_feature_id, double frac_feature_inte,
                         int sample_feature_id, double sample_feature_inte): 
    spec_id_(spec_id),
    frac_id_(frac_id),
    file_name_(file_name),
    scans_(scans),
    ms_one_id_(ms_one_id),
    ms_one_scan_(ms_one_scan),
    prec_mass_(prec_mass),
    prec_inte_(prec_inte),
    frac_feature_id_(frac_feature_id),
    frac_feature_inte_(frac_feature_inte),
    sample_feature_id_(sample_feature_id),
    sample_feature_inte_(sample_feature_inte) {
    }
    */

SpecFeature::SpecFeature(std::string line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  spec_id_ = std::stoi(strs[0]);
  frac_id_ = std::stoi(strs[1]);
  file_name_ = strs[2];
  scans_ = strs[3];
  ms_one_id_ = std::stoi(strs[4]);
  ms_one_scan_ = std::stoi(strs[5]);
  prec_mass_ = std::stod(strs[6]);
  prec_inte_ = std::stod(strs[7]);
  frac_feature_id_ = std::stoi(strs[8]);
  frac_feature_inte_ = std::stod(strs[9]);
  frac_feature_score_ = std::stod(strs[10]);
  sample_feature_id_ = std::stoi(strs[11]);
  sample_feature_inte_ = std::stod(strs[12]);
}

SpecFeature::SpecFeature(MsHeaderPtr header, FracFeaturePtr feature) {
  spec_id_ = header->getId();
  frac_id_ = header->getFractionId();
  file_name_ = header->getFileName();
  scans_ = header->getScansString();
  ms_one_id_ = header->getMsOneId();
  ms_one_scan_ = header->getMsOneScan();
  prec_mass_ = header->getPrecMonoMass();
  prec_inte_ = header->getPrecInte();
  frac_feature_id_ = feature->getId();
  frac_feature_inte_ = feature->getIntensity();
  frac_feature_score_ = feature->getPromexScore();
  sample_feature_id_ = feature->getSampleFeatureId();
  sample_feature_inte_ = feature->getSampleFeatureInte();
}

/*
bool FracFeature::cmpFracIncInteDec(const FracFeaturePtr &a, 
                                    const FracFeaturePtr &b) { 
  if (a->getFracId() < b->getFracId()) {
    return true;
  }
  else if (a->getFracId() > b->getFracId()) {
    return false;
  }
  else if (a->getIntensity() > b->getIntensity()) {
    return true;
  }
  else {
    return false;
  }
}
*/

}
