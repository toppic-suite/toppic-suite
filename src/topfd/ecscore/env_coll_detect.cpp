//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#include <ctime>
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/msalign_reader_util.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "topfd/feature_detect/feature_detect.hpp"
#include "topfd/ecscore/spectrum/peak_matrix.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"
#include "topfd/ecscore/envelope/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/feature/feature.hpp"
#include "topfd/ecscore/feature/ecscore_write_feature.hpp"
#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/env_coll_detect.hpp"

namespace toppic {

namespace env_coll_detect {

void process_single_file(std::string &base_file_name,
                         std::string &ms1_file_name, 
                         std::string &ms2_file_name, 
                         const std::string &mzml_file_name, 
                         TopfdParaPtr topfd_para_ptr,
                         EcscoreParaPtr score_para_ptr) {

  //logger::setLogLevel(LOG_LEVEL_DEBUG);  
  /// Read msalign file and get the seed envelopes.
  DeconvMsPtrVec ms1_ptr_vec;
  msalign_reader_util::readAllMsOneSpectra(ms1_file_name, ms1_ptr_vec);
  DeconvMsPtrVec ms2_ptr_vec;
  msalign_reader_util::readAllMsTwoSpectra(ms2_file_name, ms2_ptr_vec);
  LOG_DEBUG("Processed msalign file."); 

  /// Read mzml data
  RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(mzml_file_name, 
                                                                topfd_para_ptr->getActivation(),
                                                                topfd_para_ptr->getPrecWindowWidth());
  PeakPtrVec2D ms1_raw_peaks;
  std::vector<double> ms2_prec_mzs;
  raw_reader_ptr->getMs1Map(ms1_ptr_vec, ms2_ptr_vec, ms1_raw_peaks, ms2_prec_mzs); 
  raw_reader_ptr = nullptr;
  LOG_DEBUG("Processed mzML file."); 

  /// Prepare data -- seed envelopes
  SeedEnvelopePtrVec seed_ptrs;
  for (auto &ms1_data: ms1_ptr_vec) {
    DeconvPeakPtrVec peaks = ms1_data->getPeakPtrVec();
    for (auto &peak: peaks) {
      SeedEnvelopePtr seed_ptr = std::make_shared<SeedEnvelope>(peak);
      seed_ptrs.push_back(seed_ptr);
    }
  }
  std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnvelope::cmpInteDec);
  // write_out_files::write_seed_envelopes(seed_envs, "envs.csv");

  /// Prepare data -- Peak Matrix
  PeakMatrixPtr matrix_ptr = std::make_shared<PeakMatrix>(ms1_raw_peaks, ms1_ptr_vec, 
                                                          score_para_ptr->bin_size_,
                                                          topfd_para_ptr->getMsOneSnRatio());

  if (score_para_ptr->filter_neighboring_peaks_) {
    matrix_ptr->removeNonNeighbors(score_para_ptr->neighbor_mass_tole_);
  }

  /// Extract Fetures
  LOG_DEBUG("Number of seed envelopes: " << seed_ptrs.size());
  int seed_num = seed_ptrs.size();
  EnvCollPtrVec env_coll_list;

  int feat_id = 0;
  double perc = 0;
  FeaturePtrVec features;
  FracFeaturePtrVec frac_features;
  for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
    int count = seed_env_idx + 1;
    if (count % 100 == 0 || count == seed_num) {
      perc = static_cast<int>(count * 100 / seed_num);
      std::cout << "\r" << "Processing seed " << count << " ...       " << perc << "\% finished." << std::flush;
    }
    SeedEnvelopePtr seed_ptr = seed_ptrs[seed_env_idx];
    bool valid = seed_env_util::preprocessEnv(matrix_ptr, seed_ptr, 
                                              score_para_ptr, topfd_para_ptr->getMsOneSnRatio());
    if (!valid) continue;
    EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(matrix_ptr, seed_ptr,
                                                         score_para_ptr, 
                                                         topfd_para_ptr->getMsOneSnRatio());
    if (env_coll_ptr != nullptr) {
      if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                               env_coll_list, score_para_ptr)) {
        continue;
      }
      env_coll_ptr->refineMonoMass();

      FeaturePtr feat_ptr = std::make_shared<Feature>(env_coll_ptr, matrix_ptr,
                                                      feat_id, topfd_para_ptr->getMsOneSnRatio());
      feat_id++;
      features.push_back(feat_ptr);
      env_coll_ptr->removePeakData(matrix_ptr);
      if (feat_ptr->getScore() < topfd_para_ptr->getEcscoreCutoff()) {
        continue;
      }
      env_coll_ptr->setEcscore(feat_ptr->getScore());
      env_coll_list.push_back(env_coll_ptr);
      
      FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(feat_id, ms1_ptr_vec, score_para_ptr->frac_id_,
                                                                   score_para_ptr->file_name_,
                                                                   env_coll_ptr, matrix_ptr,
                                                                   topfd_para_ptr->getMsOneSnRatio());
      frac_feat_ptr->setEcscore(feat_ptr->getScore());
      frac_features.push_back(frac_feat_ptr);
    }
  }
  

  // map MS2 features
  double zero_sn_ratio = 0;
  PeakMatrixPtr raw_matrix_ptr = std::make_shared<PeakMatrix>(ms1_raw_peaks, ms1_ptr_vec, 
                                                              score_para_ptr->bin_size_,
                                                              zero_sn_ratio); 
  SpecFeaturePtrVec ms2_features;
  Feature::assignFeatures(ms2_file_name, frac_features, env_coll_list, ms2_features,
                          ms2_prec_mzs, topfd_para_ptr, score_para_ptr, 
                          raw_matrix_ptr,features, ms1_ptr_vec); 

  std::cout << std::endl << "Number of proteoform features: " << features.size() << std::endl;
  /// output files
  std::string feat_file_name = base_file_name + "_ms1.csv";
  ecscore_write_feature::writeFeatures(feat_file_name, features);

  std::string output_file_name = base_file_name + "_" + "feature.xml";
  frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);
  std::string batmass_file_name = base_file_name + "_" + "frac.mzrt.csv";
  frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);

  SampleFeaturePtrVec sample_features;
  feature_detect::getSampleFeatures(sample_features, frac_features, ms2_features);
  std::string sample_feature_file_name = base_file_name + "_"  + "ms1.feature";
  sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);
  output_file_name = base_file_name + "_"  + "ms2.feature";
  spec_feature_writer::writeFeatures(output_file_name, ms2_features);
}

void process(int frac_id, const std::string &sp_file_name, TopfdParaPtr para_ptr,
             bool is_faims, const std::vector<std::pair<double, int>> voltage_vec) {
  if (para_ptr->isMissingLevelOne()) {
    return;
  }

  std::clock_t start;
  start = std::clock();

  EcscoreParaPtr ecscore_para_ptr 
    = std::make_shared<EcscorePara>(frac_id, sp_file_name, para_ptr->getResourceDir(), para_ptr);
  std::string base_name = file_util::basename(sp_file_name);
  // if FAIME Data
  if (is_faims) { 
    for (size_t i = 0; i < voltage_vec.size(); i++) {
      std::string file_num = str_util::toString(i) + "_"; 
      std::string cur_base_name = base_name + "_" + file_num;
      std::string ms1_file_name = cur_base_name + "_ms1.msalign";
      std::string ms2_file_name = cur_base_name + "_ms2.msalign";
      process_single_file(cur_base_name, ms1_file_name, ms2_file_name, sp_file_name, para_ptr, ecscore_para_ptr);
    }
  }
  else {
    std::string ms1_file_name = base_name + "_ms1.msalign";
    std::string ms2_file_name = base_name + "_ms2.msalign";
    process_single_file(base_name, ms1_file_name, ms2_file_name, sp_file_name, para_ptr, ecscore_para_ptr);
  }

  double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
  std::cout << "Feature Extraction duration (minutes): " << duration / 60.0 << '\n';
}

}  // namespace

}  // namespace toppic
