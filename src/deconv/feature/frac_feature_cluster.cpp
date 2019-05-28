//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include "deconv/feature/frac_feature.hpp"
#include "deconv/feature/frac_feature_cluster.hpp"

namespace toppic {

namespace frac_feature_cluster {

bool matchFeature(FracFeaturePtr a, FracFeaturePtr b, 
                  double mass_tolerance, double time_tolerance) {
  double mass_diff = a->getMonoMass() - b->getMonoMass();
  if (std::abs(mass_diff) > mass_tolerance) {
    return false;
  }
  double retent_diff = a->getRetentMiddle() - b->getRetentMiddle();
  if (std::abs(retent_diff) > time_tolerance) {
    return false;
  }
  return true;
}

void setSampleFeatureId(FracFeaturePtrVec& features, 
                        double mass_tolerance, double time_tolerance) {
  std::vector<FracFeaturePtrVec> clusters; 
  std::vector<double> cluster_intes;
  for (size_t i = 0; i < features.size(); i++) {
    bool is_found = false;
    FracFeaturePtr cur_ptr = features[i];
    for (size_t j = 0; j < clusters.size(); j++) {
      FracFeaturePtr ref_ptr = clusters[j][0];
      if (matchFeature(cur_ptr, ref_ptr, mass_tolerance, time_tolerance)) {
          is_found = true;
          clusters[j].push_back(cur_ptr);
          cluster_intes[j] += cur_ptr->getIntensity();
          break;
      }
    }
    if (!is_found) {
      FracFeaturePtrVec new_clusters;
      new_clusters.push_back(cur_ptr); 
      clusters.push_back(new_clusters);
      cluster_intes.push_back(cur_ptr->getIntensity());
    }
  }
  for (size_t i = 0; i < clusters.size(); i++) {
    for (size_t j = 0; j < clusters[i].size(); j++) {
      clusters[i][j]->setSampleFeatureId(i);
      clusters[i][j]->setSampleFeatureInte(cluster_intes[i]);
    }
  }
}

void cluster(FracFeaturePtrVec &features, double mass_tolerance, 
             double time_tolerance) {
  std::sort(features.begin(), features.end(), FracFeature::cmpInteDec);
  setSampleFeatureId(features, mass_tolerance, time_tolerance);
  std::sort(features.begin(), features.end(), FracFeature::cmpFracIncInteDec);
}

}

}  // namespace toppic
