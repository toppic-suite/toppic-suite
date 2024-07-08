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

#include <set>
#include <map>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "para/sp_para.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/feature/spec_feature_reader.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/frac_feature_reader.hpp"
#include "ms/feature/frac_xml_feature_reader.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/frac_feature_cluster.hpp"
#include "ms/feature/sample_feature.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "ms/feature/feature_merge.hpp"

namespace toppic {

namespace feature_merge {

void mergeFiles(const std::vector<std::string> &frac_feature_xml_file_list,
                const std::string &frac_feature_xml_output_file_name, 
                const std::string &frac_feature_tsv_output_file_name, 
                const std::vector<std::string> &spec_feature_file_lst,
                const std::string &spec_feature_output_file_name,
                const std::string &sample_feature_output_file_name,
                int max_spec_num_per_file, 
                int max_feature_num_per_file, 
                const std::string &para_str) {
  FracFeaturePtrVec all_frac_features;
  for (size_t i = 0; i < frac_feature_xml_file_list.size(); i++) {
    FracXmlFeatureReader ft_reader(frac_feature_xml_file_list[i]); 
    FracFeaturePtrVec features = ft_reader.readAllFeatures();
    ft_reader.close();
    for (size_t j = 0; j < features.size(); j++) {
      int feature_id = features[j]->getFeatId() + i * max_feature_num_per_file;
      features[j]->setFracId(i);
      features[j]->setFeatId(feature_id);
    }
    all_frac_features.insert(all_frac_features.end(), features.begin(), features.end());
  }

  FracFeaturePtrVec2D clusters;
  //In simpleCluster, we do not merge fraction features with similar 
  //masses and charge states. 
  frac_feature_cluster::simpleCluster(all_frac_features, clusters);
  //double mass_tolerance = 0.2;
  //double time_tolerance = 600;
  //frac_feature_cluster::cluster(all_frac_features, clusters, mass_tolerance, time_tolerance);
  
  //sample features;
  SampleFeaturePtrVec sample_features;
  for (size_t i = 0; i < clusters.size(); i++) {
    SampleFeaturePtr sample_feature = std::make_shared<SampleFeature>(clusters[i][0], clusters[i][0]->getFeatId());
    sample_features.push_back(sample_feature);
  }
  sample_feature_writer::writeFeatures(sample_feature_output_file_name, sample_features);

  frac_feature_writer::writeXmlFeatures(frac_feature_xml_output_file_name, all_frac_features);
  frac_feature_writer::writeFeatures(frac_feature_tsv_output_file_name, all_frac_features);

  //spec features
  std::map<int,FracFeaturePtr> feature_map;
  for (size_t i = 0; i < all_frac_features.size(); i++) {
    feature_map[all_frac_features[i]->getFeatId()] =  all_frac_features[i];
  }
  SpecFeaturePtrVec all_spec_features;
  for (size_t i = 0; i < spec_feature_file_lst.size(); i++) {
    SpecFeatureReader ft_reader(spec_feature_file_lst[i]);
    SpecFeaturePtrVec spec_features = ft_reader.readAllFeatures();
    for (size_t j = 0; j < spec_features.size(); j++) {
      SpecFeaturePtr spec_feature = spec_features[j];
      int spec_id = spec_feature->getSpecId() + i * max_spec_num_per_file;
      spec_feature->setSpecId(spec_id);
      int ms_one_id = spec_feature->getMsOneId() + i * max_spec_num_per_file;
      spec_feature->setMsOneId(ms_one_id);
      int frac_feature_id = spec_feature->getFracFeatureId()+ i * max_feature_num_per_file;
      spec_feature->setFracFeatureId(frac_feature_id);
      FracFeaturePtr frac_feature = feature_map.find(frac_feature_id)->second;
    }
    all_spec_features.insert(all_spec_features.end(), spec_features.begin(), spec_features.end());
  }

  spec_feature_writer::writeFeatures(spec_feature_output_file_name, all_spec_features);
}

void process(const std::vector<std::string> &raw_file_names,
             const std::string &output_file_name, 
             std::string &para_str) {
  std::vector<std::string> frac_feature_xml_names;
  std::vector<std::string> frac_feature_names;
  std::vector<std::string> spec_feature_names;
  for (size_t i = 0; i < raw_file_names.size(); i++) { 
    std::string base_name = file_util::basename(raw_file_names[i]);
    std::string frac_feature_xml = base_name + "_frac_feature.xml";
    frac_feature_xml_names.push_back(frac_feature_xml);
    std::string spec_feature = base_name + "_ms2.feature";
    spec_feature_names.push_back(spec_feature);
  }
  
  std::string frac_feature_xml_output_name = output_file_name + "_frac_feature.xml";
  std::string frac_feature_tsv_output_name = output_file_name + "_ms1.frac_feature";
  std::string spec_feature_output_name = output_file_name + "_ms2.feature";
  std::string sample_feature_output_name = output_file_name + "_ms1.feature";

  mergeFiles(frac_feature_xml_names, frac_feature_xml_output_name, 
             frac_feature_tsv_output_name,
             spec_feature_names, spec_feature_output_name,
             sample_feature_output_name,
             SpPara::getMaxSpecNumPerFile(), 
             SpPara::getMaxFeatureNumPerFile(), 
             para_str); 
}


}

} /* namespace toppic */
