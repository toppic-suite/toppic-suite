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

#include <map>

#include "ms/feature/feature_util.hpp"

namespace toppic {

namespace feature_util {

void getSampleFeatures(SampleFeaturePtrVec &sample_features, FracFeaturePtrVec &frac_features,
                       SpecFeaturePtrVec &spec_features) {
  //sample features;
  for (size_t i = 0; i < frac_features.size(); i++) {
    SampleFeaturePtr sample_feature = std::make_shared<SampleFeature>(frac_features[i], frac_features[i]->getFeatId());
    sample_features.push_back(sample_feature);
  }

  //spec features
  std::map<int,FracFeaturePtr> feature_map;
  for (size_t i = 0; i < frac_features.size(); i++) {
    feature_map[frac_features[i]->getFeatId()] =  frac_features[i];
  }

  for (size_t j = 0; j < spec_features.size(); j++) {
    SpecFeaturePtr spec_feature = spec_features[j];
    int frac_feature_id = spec_feature->getFracFeatureId(); 
    FracFeaturePtr frac_feature = feature_map.find(frac_feature_id)->second;
  }
}

}  // namespace 

}  // namespace toppic 
