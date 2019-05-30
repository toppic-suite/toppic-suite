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
#include <map>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "deconv/feature/frac_feature_reader.hpp"
#include "deconv/feature/frac_feature_writer.hpp"
#include "deconv/feature/frac_feature_cluster.hpp"
#include "deconv/feature/frac_ms2_feature.hpp"
#include "deconv/feature/frac_ms2_feature_reader.hpp"
#include "deconv/feature/frac_ms2_feature_writer.hpp"
#include "deconv/feature/frac_feature_merge.hpp"

namespace toppic {

namespace frac_feature_merge {

void mergeFiles(const std::vector<std::string> &feature_file_lst,
                const std::string &feature_output_file_name, 
                const std::vector<std::string> &ms2_feature_file_lst,
                const std::string &ms2_feature_output_file_name,
                int max_num_per_file, 
                const std::string &para_str) {
  FracFeaturePtrVec all_features;
  for (size_t i = 0; i < feature_file_lst.size(); i++) {
    FracFeatureReader ft_reader(feature_file_lst[i]); 
    FracFeaturePtrVec features = ft_reader.readAllFeatures();
    ft_reader.close();
    for (size_t j = 0; j < features.size(); j++) {
      int feature_id = features[j]->getId() + i * max_num_per_file;
      features[j]->setId(feature_id);
      //frac_feature_writer::writeOneFeature(outfile, features[j]);
    }
    all_features.insert(all_features.end(), features.begin(), features.end());
  }

  double mass_tolerance = 0.2;
  double time_tolerance = 600;
  frac_feature_cluster::cluster(all_features, mass_tolerance, time_tolerance);
  frac_feature_writer::writeFeatures(feature_output_file_name, all_features);

  std::map<int,FracFeaturePtr> feature_map;
  for (size_t i = 0; i < all_features.size(); i++) {
    feature_map[all_features[i]->getId()] =  all_features[i];
  }

  FracMs2FeaturePtrVec all_ms2_features;
  for (size_t i = 0; i < ms2_feature_file_lst.size(); i++) {
    FracMs2FeatureReader ft_reader(ms2_feature_file_lst[i]);
    FracMs2FeaturePtrVec ms2_features = ft_reader.readAllFeatures();
    for (size_t j = 0; j < ms2_features.size(); j++) {
      FracMs2FeaturePtr ms2_feature = ms2_features[j];
      int id = ms2_feature->getId() + i * max_num_per_file;
      ms2_feature->setId(id);
      int ms_one_id = ms2_feature->getMsOneId() + i * max_num_per_file;
      ms2_feature->setMsOneId(ms_one_id);
      int frac_feature_id = ms2_feature->getFracFeatureId()+ i * max_num_per_file;
      ms2_feature->setFracFeatureId(frac_feature_id);
      FracFeaturePtr ms1_feature = feature_map.find(frac_feature_id)->second;
      ms2_feature->setSampleFeatureId(ms1_feature->getSampleFeatureId());
      ms2_feature->setSampleFeatureInte(ms1_feature->getSampleFeatureInte());
    }
    all_ms2_features.insert(all_ms2_features.end(), ms2_features.begin(), ms2_features.end());
  }

  frac_ms2_feature_writer::writeFeatures(ms2_feature_output_file_name, all_ms2_features);
}

}

} /* namespace toppic */
