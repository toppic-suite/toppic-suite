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

#include <algorithm>

#include "common/util/file_util.hpp"
#include "ms/feature/frac_feature.hpp"
#include "ms/feature/frac_feature_cluster.hpp"

namespace toppic {

namespace frac_feature_cluster {

bool matchFeature(FracFeaturePtr a, FracFeaturePtr b, 
                  double mass_tolerance, double time_tolerance) {
  double mass_diff = a->getMonoMass() - b->getMonoMass();
  if (std::abs(mass_diff) > mass_tolerance) {
    return false;
  }
  double retent_diff = a->getTimeMiddle() - b->getTimeMiddle();
  if (std::abs(retent_diff) > time_tolerance) {
    return false;
  }
  if (a->getFracId() == b->getFracId()) {
    return false;
  }
  return true;
}

void getClusters(FracFeaturePtrVec& features, 
                 FracFeaturePtrVec2D& clusters,
                 double mass_tolerance, double time_tolerance) {
  for (size_t i = 0; i < features.size(); i++) {
    bool is_found = false;
    FracFeaturePtr cur_ptr = features[i];
    for (size_t j = 0; j < clusters.size(); j++) {
      for (size_t k = 0; k < clusters[j].size(); k++) {
        FracFeaturePtr ref_ptr = clusters[j][k];
        if (matchFeature(cur_ptr, ref_ptr, mass_tolerance, time_tolerance)) {
          is_found = true;
          break;
        }
      }
      if (is_found) {
        clusters[j].push_back(cur_ptr);
        break;
      }
    }
    if (!is_found) {
      FracFeaturePtrVec new_cluster;
      new_cluster.push_back(cur_ptr); 
      clusters.push_back(new_cluster);
    }
  }
}

void cluster(FracFeaturePtrVec &features, FracFeaturePtrVec2D &clusters,
             double mass_tolerance, double time_tolerance) {
  std::sort(features.begin(), features.end(), FracFeature::cmpInteDec);
  getClusters(features, clusters, mass_tolerance, time_tolerance);
  std::sort(features.begin(), features.end(), FracFeature::cmpFracIncInteDec);
}

void simpleCluster(FracFeaturePtrVec &features, FracFeaturePtrVec2D &clusters) {
  std::sort(features.begin(), features.end(), FracFeature::cmpInteDec);
  for (size_t i = 0; i < features.size(); i++) {
    FracFeaturePtrVec new_cluster;
    new_cluster.push_back(features[i]);
    clusters.push_back(new_cluster);
  }
  std::sort(features.begin(), features.end(), FracFeature::cmpFracIncInteDec);
}

}

}  // namespace toppic
