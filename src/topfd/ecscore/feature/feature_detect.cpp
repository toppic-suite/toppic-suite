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

#include <map>

#include "topfd/ecscore/feature/feature_detect.hpp"

namespace toppic {

namespace feature_detect {

void getSampleFeatures(SampleFeaturePtrVec &sample_features, FracFeaturePtrVec &frac_features,
                       SpecFeaturePtrVec &spec_features) {
  //sample features;
  for (size_t i = 0; i < frac_features.size(); i++) {
    SampleFeaturePtr sample_feature = std::make_shared<SampleFeature>(frac_features[i], frac_features[i]->getId());
    sample_features.push_back(sample_feature);
    frac_features[i]->setSampleFeatureId(frac_features[i]->getId());
    frac_features[i]->setSampleFeatureInte(frac_features[i]->getIntensity());
  }

  //spec features
  std::map<int,FracFeaturePtr> feature_map;
  for (size_t i = 0; i < frac_features.size(); i++) {
    feature_map[frac_features[i]->getId()] =  frac_features[i];
  }

  for (size_t j = 0; j < spec_features.size(); j++) {
    SpecFeaturePtr spec_feature = spec_features[j];
    int frac_feature_id = spec_feature->getFracFeatureId(); 
    FracFeaturePtr frac_feature = feature_map.find(frac_feature_id)->second;
    spec_feature->setSampleFeatureId(frac_feature->getSampleFeatureId());
    spec_feature->setSampleFeatureInte(frac_feature->getSampleFeatureInte());
  }
}

}  // namespace 

}  // namespace toppic 
