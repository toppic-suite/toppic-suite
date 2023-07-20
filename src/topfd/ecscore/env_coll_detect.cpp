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
#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_reader_util.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/ecscore/spectrum/peak_matrix.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"
#include "topfd/ecscore/envelope/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/feature/feature.hpp"
#include "topfd/ecscore/feature/feature_util.hpp"
#include "topfd/ecscore/feature/ecscore_write_feature.hpp"
#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/env_coll_detect.hpp"

namespace toppic {

namespace env_coll_detect {

void process(TopfdParaPtr topfd_para_ptr) {
  if (topfd_para_ptr->isMissingLevelOne()) {
    return;
  }
  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(topfd_para_ptr->getFracId(), 
                                                                topfd_para_ptr->getMzmlFileName(),
                                                                topfd_para_ptr->getResourceDir(), 
                                                                topfd_para_ptr);
  // read deconvoluted MS1 peaks
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();
  std::string ms1_file_name = output_base_name + "_ms1.msalign";
  DeconvMsPtrVec deconv_ms1_ptr_vec;
  msalign_reader_util::readAllSpectra(ms1_file_name, deconv_ms1_ptr_vec);

  // read ms1 raw peaks and ms2_headers
  PeakPtrVec2D ms1_mzml_peaks;
  MsHeaderPtr2D ms2_header_ptr_2d;
  MzmlMsGroupReaderPtr mzml_reader_ptr = 
    std::make_shared<MzmlMsGroupReader>(topfd_para_ptr->getMzmlFileName(), 
                                        topfd_para_ptr->getPrecWindowWidth(),
                                        topfd_para_ptr->getActivation(),
                                        topfd_para_ptr->getFracId(),
                                        topfd_para_ptr->isFaims(), 
                                        topfd_para_ptr->getFaimsVoltage(), 
                                        topfd_para_ptr->isMissingLevelOne());
  mzml_reader_ptr->getMs1Map(ms1_mzml_peaks, ms2_header_ptr_2d); 
  mzml_reader_ptr = nullptr;

  //Prepare seed envelopes
  SeedEnvelopePtrVec seed_ptrs;
  SeedEnvelopePtr2D seed_ptr_2d;
  for (auto &ms1_data: deconv_ms1_ptr_vec) {
    DeconvPeakPtrVec peaks = ms1_data->getPeakPtrVec();
    SeedEnvelopePtrVec one_spec_seed_ptrs;
    for (auto &peak: peaks) {
      SeedEnvelopePtr seed_ptr_1 = std::make_shared<SeedEnvelope>(peak);
      seed_ptrs.push_back(seed_ptr_1);
      SeedEnvelopePtr seed_ptr_2 = std::make_shared<SeedEnvelope>(peak);
      one_spec_seed_ptrs.push_back(seed_ptr_2);
    }
    seed_ptr_2d.push_back(one_spec_seed_ptrs);
  }
  std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnvelope::cmpInteDec);
  // write_out_files::write_seed_envelopes(seed_envs, "envs.csv");

  double sn_ratio = topfd_para_ptr->getMsOneSnRatio();
  /// Prepare data -- Peak Matrix
  PeakMatrixPtr matrix_ptr = std::make_shared<PeakMatrix>(ms1_mzml_peaks, deconv_ms1_ptr_vec, 
                                                          score_para_ptr->bin_size_, sn_ratio); 

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
                                              score_para_ptr, sn_ratio); 
    if (!valid) continue;
    EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(matrix_ptr, seed_ptr,
                                                         score_para_ptr, sn_ratio); 
    if (env_coll_ptr != nullptr) {
      if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                               env_coll_list, score_para_ptr)) {
        continue;
      }
      env_coll_ptr->refineMonoMass();
      FeaturePtr feat_ptr = std::make_shared<Feature>(env_coll_ptr, matrix_ptr,
                                                      feat_id, sn_ratio); 
      feat_id++;
      features.push_back(feat_ptr);
      env_coll_ptr->removePeakData(matrix_ptr);
      if (feat_ptr->getScore() < topfd_para_ptr->getEcscoreCutoff()) {
        continue;
      }
      env_coll_ptr->setEcscore(feat_ptr->getScore());
      env_coll_list.push_back(env_coll_ptr);
      FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(feat_id, deconv_ms1_ptr_vec, 
                                                                   score_para_ptr->frac_id_,
                                                                   score_para_ptr->file_name_,
                                                                   env_coll_ptr, matrix_ptr, sn_ratio);
      frac_feat_ptr->setEcscore(feat_ptr->getScore());
      frac_features.push_back(frac_feat_ptr);
    }
  }
  
  /*
  for (size_t z = 0; z < env_coll_list.size(); z++) {
    double mass = env_coll_list[z]->getMass();
    if (mass > 10059.3 && mass < 10059.4) {
      LOG_ERROR("Step 3 Mass " << mass << " env set list length " <<
                env_coll_list[z]->getEnvSetList().size()); 
    }
  }
  */

  // map MS2 features
  double zero_sn_ratio = 0;
  PeakMatrixPtr raw_matrix_ptr = std::make_shared<PeakMatrix>(ms1_mzml_peaks, deconv_ms1_ptr_vec, 
                                                              score_para_ptr->bin_size_,
                                                              zero_sn_ratio); 
  SpecFeaturePtrVec ms2_features;

  Feature::assignFeatures(frac_features, env_coll_list, features, 
                          raw_matrix_ptr, deconv_ms1_ptr_vec, ms2_header_ptr_2d,
                          seed_ptr_2d, ms2_features, topfd_para_ptr, score_para_ptr); 

  std::cout << std::endl << "Number of proteoform features: " << features.size() << std::endl;
  /// output files
  std::string feat_file_name = output_base_name + "_ms1.csv";
  ecscore_write_feature::writeFeatures(feat_file_name, features);

  std::string output_file_name = output_base_name + "_" + "feature.xml";
  frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);
  std::string batmass_file_name = output_base_name + "_" + "frac.mzrt.csv";
  frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);

  SampleFeaturePtrVec sample_features;
  feature_util::getSampleFeatures(sample_features, frac_features, ms2_features);
  std::string sample_feature_file_name = output_base_name + "_"  + "ms1.feature";
  sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);
  output_file_name = output_base_name + "_"  + "ms2.feature";
  spec_feature_writer::writeFeatures(output_file_name, ms2_features);
}

}  // namespace

}  // namespace toppic
