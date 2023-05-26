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
#include "ms/spec/simple_msalign_reader.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
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
  SimpleMsAlignReader::readMsOneSpectra(ms1_file_name, ms1_ptr_vec);
  DeconvMsPtrVec ms2_ptr_vec;
  SimpleMsAlignReader::readMsTwoSpectra(ms2_file_name, ms2_ptr_vec);
  LOG_DEBUG("Processed msalign file."); 

  /// Read mzml data
  RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(mzml_file_name, 
                                                                topfd_para_ptr->getActivation(),
                                                                topfd_para_ptr->getPrecWindow());
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
      
      FracFeaturePtr feature_ptr = env_coll_util::getFracFeature(feat_id, ms1_ptr_vec, score_para_ptr->frac_id_,
                                                                 score_para_ptr->file_name_,
                                                                 env_coll_ptr, matrix_ptr,
                                                                 topfd_para_ptr->getMsOneSnRatio());
      /*
      feature_ptr->setPromexScore(feature.getScore());
      frac_features.push_back(feature_ptr);
      */
    }
  }
  /**
  // map MS2 features
  SpecFeaturePtrVec ms2_features;
  Feature::assign_features(ms1_ptr_vec, ms2_file_name, frac_features, env_coll_list, features, ms2_features,
                           prec_spectrum_Base_mono_mz, peak_matrix, model, model_escore, feature_para_ptr, para_ptr);
                           */
  std::cout << std::endl << "Number of proteoform features: " << features.size() << std::endl;

  /// output files
  std::string feat_file_name = base_file_name + "_ms1.csv";
  ecscore_write_feature::writeFeatures(feat_file_name, features);

  /*
  std::string output_file_name = base_name + "_" + file_num + "feature.xml";
  frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);

  std::string batmass_file_name = base_name + "_" + file_num + "frac.mzrt.csv";
  frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);

  SampleFeaturePtrVec sample_features;
  feature_detect_old::getSampleFeatures(sample_features, frac_features, ms2_features);
  std::string sample_feature_file_name = base_name + "_" + file_num + "ms1.feature";
  sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);
  output_file_name = base_name + "_" + file_num + "ms2.feature";
  spec_feature_writer::writeFeatures(output_file_name, ms2_features);
  */
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


  /**

  // read ms1 deconvoluted spectra
  if (!para_ptr->isMissingLevelOne()) {
    std::string file_num;
    for (size_t i = 0; i < voltage_vec.size(); i++) {
      if (isFaims) { file_num = str_util::toString(i) + "_"; } // if FAIME Data

      std::string ms1_file_name = base_name + "_" + file_num + "ms1.msalign";
      std::string ms2_file_name = base_name + "_" + file_num + "ms2.msalign";

      FracFeaturePtrVec frac_features;

      // EnvCNN
      std::string file_name = para_ptr->getResourceDir()
        + file_util::getFileSeparator() + "envcnn_models"
        + file_util::getFileSeparator() + "envcnn_2_block_model.json";
      fdeep::model model = fdeep::load_model(file_name, true, fdeep::dev_null_logger);

      // ECScore
      file_name = para_ptr->getResourceDir()
        + file_util::getFileSeparator() + "envcnn_models"
        + file_util::getFileSeparator() + "escore.json";
      fdeep::model model_escore = fdeep::load_model(file_name, true, fdeep::dev_null_logger);

      /// Read msalign file and get the seed envelopes.
      DeconvMsPtrVec ms1_ptr_vec;
      SimpleMsAlignReader::readMsOneSpectra(ms1_file_name, ms1_ptr_vec);
      std::cout << "Processed msalign file." << std::endl;

      /// Read mzml data
      double cur_voltage = voltage_vec[i].first;//if this is -1, it is non-FAIME data
      PeakPtrVec2D raw_peaks;
      std::vector<double> prec_spectrum_Base_mono_mz;
      RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(sp_file_name, para_ptr->getActivation(),
                                                                    para_ptr->getPrecWindow());
      raw_reader_ptr->getMsData(raw_peaks, prec_spectrum_Base_mono_mz, cur_voltage);
      raw_reader_ptr = nullptr;
      std::cout << "Processed mzML file." << std::endl;

      /// Prepare data -- seed envelopes
      std::vector<SeedEnvelope> seed_envs;
      for (auto &ms1_data: ms1_ptr_vec) {
        std::vector<DeconvPeakPtr> peaks = ms1_data->getPeakPtrVec();
        for (auto &peak: peaks)
          seed_envs.push_back(SeedEnvelope(peak));
      }
      std::sort(seed_envs.begin(), seed_envs.end(), SeedEnvelope::cmpInteDec);
      // write_out_files::write_seed_envelopes(seed_envs, "envs.csv");

      /// Prepare data -- Peak Matrix
      PeakMatrix peak_matrix = PeakMatrix(raw_peaks, ms1_ptr_vec, feature_para_ptr->bin_size_,
                                          para_ptr->getMsOneSnRatio());
      if (feature_para_ptr->filter_neighboring_peaks_)
        peak_matrix.find_remove_non_neighbors(feature_para_ptr->neighbor_mass_tole_);

      /// Extract Fetures
      std::cout << "Number of seed envelopes: " << seed_envs.size() << std::endl;
      int seed_num = seed_envs.size();
      int env_coll_num = 0;
      std::vector<EnvCollection> env_coll_list;
      std::vector<Feature> features;
      for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
        if (seed_env_idx % 10000 == 0)
          std::cout << "\r" << "Processing peak " << seed_env_idx << " and Features found " << env_coll_num << std::flush;
        SeedEnvelope env = seed_envs[seed_env_idx];
        bool valid = false;
        valid = evaluate_envelope::preprocess_env(peak_matrix, env, feature_para_ptr, para_ptr->getMsOneSnRatio());
        if (!valid) continue;
        EnvCollection env_coll = env_coll_util::find_env_collection(peak_matrix, env, feature_para_ptr, para_ptr->getMsOneSnRatio());
        if (!env_coll.isEmpty()) {
          if (env_coll_util::check_in_existing_features(peak_matrix, env_coll, env_coll_list, feature_para_ptr))
            continue;
          env_coll.refine_mono_mass();

          Feature feature = Feature(env_coll, peak_matrix, model, model_escore, env_coll_num,
                                    para_ptr->getMsOneSnRatio());
          if (feature.getScore() < para_ptr->getECScore()) continue;
          features.push_back(feature);
          env_coll.setEcscore(feature.getScore());
          env_coll.remove_peak_data(peak_matrix);
          env_coll_list.push_back(env_coll);
          FracFeaturePtr feature_ptr = Feature::getFeature(env_coll_num, ms1_ptr_vec, feature_para_ptr->frac_id_,
                                                           feature_para_ptr->file_name_, env_coll, peak_matrix,
                                                           para_ptr->getMsOneSnRatio());
          feature_ptr->setPromexScore(feature.getScore());
          frac_features.push_back(feature_ptr);
          env_coll_num = env_coll_num + 1;
        }
      }
      // map MS2 features
      SpecFeaturePtrVec ms2_features;
      Feature::assign_features(ms1_ptr_vec, ms2_file_name, frac_features, env_coll_list, features, ms2_features,
                               prec_spectrum_Base_mono_mz, peak_matrix, model, model_escore, feature_para_ptr, para_ptr);
      std::cout << std::endl << "Number of Envelope Collections: " << features.size() << std::endl;

      /// output files
      file_name = base_name + "_" + file_num + "ms1.csv";
      write_feature::writeFeatures(file_name, features);

      std::string output_file_name = base_name + "_" + file_num + "feature.xml";
      frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);

      std::string batmass_file_name = base_name + "_" + file_num + "frac.mzrt.csv";
      frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);

      SampleFeaturePtrVec sample_features;
      feature_detect_old::getSampleFeatures(sample_features, frac_features, ms2_features);
      std::string sample_feature_file_name = base_name + "_" + file_num + "ms1.feature";
      sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);
      output_file_name = base_name + "_" + file_num + "ms2.feature";
      spec_feature_writer::writeFeatures(output_file_name, ms2_features);

    }
  }
}
**/

