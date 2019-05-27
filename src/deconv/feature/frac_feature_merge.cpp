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

#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "deconv/feature/frac_feature_reader.hpp"
#include "deconv/feature/frac_feature_writer.hpp"
#include "deconv/feature/frac_feature_cluster.hpp"
#include "deconv/feature/frac_feature_merge.hpp"

namespace toppic {

namespace frac_feature_merge {

void mergeFiles(const std::vector<std::string> &spec_file_lst,
                const std::string &output_file_name, 
                int max_num_per_file,
                const std::string &para_str) {

  FracFeaturePtrVec all_features;
  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    FracFeatureReader ft_reader(spec_file_lst[i]); 
    FracFeaturePtrVec features = ft_reader.readAllFeatures();
    ft_reader.close();
    for (size_t j = 0; j < features.size(); j++) {
      int feature_id = features[j]->getId() + i * max_num_per_file;
      features[j]->setId(feature_id);
      //frac_feature_writer::writeOneFeature(outfile, features[j]);
    }
    all_features.insert(features.begin(), features.end(), all_features.end());
  }

  double mass_tolerance = 0.2;
  double time_tolerance = 600;
  frac_feature_cluster::cluster(all_features, mass_tolerance, time_tolerance);

  frac_feature_writer::writeFeatures(output_file_name, all_features);
}

}

} /* namespace toppic */
