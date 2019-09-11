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

#ifndef TOPPIC_FEATURE_FRAC_FEATURE_CLUSTER_HPP_
#define TOPPIC_FEATURE_FRAC_FEATURE_CLUSTER_HPP_

#include "feature/frac_feature.hpp"

namespace toppic {

namespace frac_feature_cluster {

void cluster(FracFeaturePtrVec &features, FracFeaturePtrVec2D &clusters,
             double mass_tolerance, double time_tolerance);

}

}  // namespace toppic

#endif
