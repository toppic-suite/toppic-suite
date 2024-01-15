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

#include <ctime>
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_reader_util.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/sample_feature_writer.hpp"
#include "ms/feature/feature_util.hpp"
#include "ms/msmap/ms_map.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/score/ecscore.hpp"
#include "topfd/ecscore/score/ecscore_writer.hpp"
#include "topfd/ecscore/env_coll/env_coll_assign.hpp"
#include "topfd/ecscore/env_coll/env_coll_detect.hpp"

namespace toppic {

namespace env_coll_detect {

void process(TopfdParaPtr topfd_para_ptr) {
  if (topfd_para_ptr->isMissingLevelOne()) {
    return;
  }
  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(topfd_para_ptr->getFracId(), 
                                                                topfd_para_ptr->getMzmlFileName(),
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
  SeedEnvPtrVec seed_ptrs;
  SeedEnvPtr2D seed_ptr_2d;
  for (size_t i = 0; i < deconv_ms1_ptr_vec.size(); i++) {
    DeconvMsPtr ms1_data = deconv_ms1_ptr_vec[i]; 
    DeconvPeakPtrVec peaks = ms1_data->getPeakPtrVec();
    SeedEnvPtrVec one_spec_seed_ptrs;
    for (auto &peak: peaks) {
      SeedEnvPtr seed_ptr_1 = std::make_shared<SeedEnv>(i, peak);
      seed_ptrs.push_back(seed_ptr_1);
      SeedEnvPtr seed_ptr_2 = std::make_shared<SeedEnv>(i, peak);
      one_spec_seed_ptrs.push_back(seed_ptr_2);
    }
    seed_ptr_2d.push_back(one_spec_seed_ptrs);
  }

  std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnv::cmpSeedInteDec);
  // write_out_files::write_seed_envelopes(seed_envs, "envs.csv");

  double sn_ratio = topfd_para_ptr->getMsOneSnRatio();
  bool single_scan_noise = topfd_para_ptr->isUseSingleScanNoiseLevel();
  /// Prepare data -- Peak Matrix
  MsMapPtr matrix_ptr = std::make_shared<MsMap>(ms1_mzml_peaks, deconv_ms1_ptr_vec,
                                                score_para_ptr->bin_size_,
                                                sn_ratio, single_scan_noise);

  if (score_para_ptr->min_scan_num_ >= 2) {
    matrix_ptr->removeNonNeighbors(score_para_ptr->neighbor_mz_tole_);
  }

  /// Extract Fetures
  LOG_DEBUG("Number of seed envelopes: " << seed_ptrs.size());
  int seed_num = seed_ptrs.size();
  EnvCollPtrVec env_coll_list;
  int feat_id = 0;
  double perc = 0;
  ECScorePtrVec ecscore_list;
  FracFeaturePtrVec frac_features;
  for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
    int count = seed_env_idx + 1;
    if (count % 100 == 0 || count == seed_num) {
      perc = static_cast<int>(count * 100 / seed_num);
      std::cout << "\r" << "Processing feature " << count << " ...       " << perc << "\% finished." << std::flush;
    }
    SeedEnvPtr seed_ptr = seed_ptrs[seed_env_idx];
    seed_ptr = seed_env_util::preprocessSeedEnvPtr(seed_ptr, matrix_ptr,  
                                                   score_para_ptr, sn_ratio); 
    if (seed_ptr == nullptr) continue;
    EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(matrix_ptr, seed_ptr,
                                                         score_para_ptr, sn_ratio); 
    if (env_coll_ptr == nullptr) continue;
    if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                             env_coll_list, score_para_ptr)) {
      env_coll_ptr->removePeakData(matrix_ptr);
      continue;
    }
    env_coll_ptr->refineMonoMass();
    ECScorePtr ecscore_ptr = std::make_shared<ECScore>(env_coll_ptr, matrix_ptr,
                                                       feat_id, sn_ratio); 
    if (ecscore_ptr->getScore() < topfd_para_ptr->getEcscoreCutoff()) {
      continue;
    }
    ecscore_list.push_back(ecscore_ptr);
    env_coll_ptr->setEcscore(ecscore_ptr->getScore());
    env_coll_ptr->removePeakData(matrix_ptr);
    env_coll_list.push_back(env_coll_ptr);
    FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(feat_id, deconv_ms1_ptr_vec, 
                                                                 score_para_ptr->frac_id_,
                                                                 score_para_ptr->file_name_,
                                                                 env_coll_ptr, matrix_ptr, sn_ratio);
    frac_features.push_back(frac_feat_ptr);
    feat_id++;
  }
  std::cout << std::endl; 

  if (topfd_para_ptr->isSearchPrecWindow()) {
    // add ms1 feature based on precursor windows
    // set min match envelope to 1 to accept single scan features
    score_para_ptr->min_scan_num_ = 1;
    score_para_ptr->min_match_peak_ = 1;
    sn_ratio = 0;
    matrix_ptr->reconstruct(sn_ratio, single_scan_noise); 
    int ms1_spec_num = deconv_ms1_ptr_vec.size();
    for (size_t ms1_idx = 0; ms1_idx < deconv_ms1_ptr_vec.size(); ms1_idx++) {
      perc = static_cast<int>((ms1_idx + 1)* 100 / ms1_spec_num);
      int scan = deconv_ms1_ptr_vec[ms1_idx]->getMsHeaderPtr()->getFirstScanNum();
      std::cout << "\r" << "Additional feature search MS1 spectrum scan " 
        << scan << " ...       " << perc << "\% finished." << std::flush;
      for (size_t i = 0; i < ms2_header_ptr_2d[ms1_idx].size(); i++) {
        MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_idx][i];
        if (env_coll_assign::checkEnvColl(ms1_idx, header_ptr, env_coll_list)) {
          continue;
        }
        double prec_win_begin = header_ptr->getPrecWinBegin();
        double prec_win_end = header_ptr->getPrecWinEnd();
        SeedEnvPtrVec seed_ptr_list = seed_ptr_2d[ms1_idx];
        SeedEnvPtrVec selected_seed_list;
        LOG_DEBUG("ms1 id  " << ms1_idx << " seed number " << seed_ptr_list.size());
        for (size_t i = 0; i < seed_ptr_list.size(); i++) {
          double ref_mz = seed_ptr_list[i]->getReferMz();
          if (ref_mz > prec_win_begin && ref_mz < prec_win_end) {
            selected_seed_list.push_back(seed_ptr_list[i]);
          }
        }
        SeedEnvPtr seed_ptr;
        if (selected_seed_list.size() > 0) {
          // choose the highest intensity one
          std::sort(selected_seed_list.begin(), selected_seed_list.end(),
                    SeedEnv::cmpSeedInteDec);
          seed_ptr = selected_seed_list[0];  
          seed_ptr = seed_env_util::relaxProcessSeedEnvPtr(seed_ptr, matrix_ptr,  
                                                           score_para_ptr, sn_ratio); 
        }
        if (seed_ptr == nullptr) {
          continue;
        }

        EnvCollPtr env_coll_ptr = env_coll_util::findEnvCollWithSingleEnv(matrix_ptr, seed_ptr, score_para_ptr,
                                                                          sn_ratio); 
        if (env_coll_ptr == nullptr) continue;
        if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                                 env_coll_list, score_para_ptr)) {
          env_coll_ptr->removePeakData(matrix_ptr);
          continue;
        }
        env_coll_ptr->refineMonoMass();
        ECScorePtr ecscore_ptr = std::make_shared<ECScore>(env_coll_ptr, matrix_ptr,
                                                           feat_id, sn_ratio); 
        if (ecscore_ptr->getScore() < 0 || std::isnan(ecscore_ptr->getScore())) {
          continue;
        }
        ecscore_list.push_back(ecscore_ptr);
        env_coll_ptr->setEcscore(ecscore_ptr->getScore());
        env_coll_ptr->removePeakData(matrix_ptr);
        env_coll_list.push_back(env_coll_ptr);
        FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(feat_id, deconv_ms1_ptr_vec, 
                                                                     score_para_ptr->frac_id_,
                                                                     score_para_ptr->file_name_,
                                                                     env_coll_ptr, matrix_ptr, sn_ratio);
        frac_features.push_back(frac_feat_ptr);
        feat_id++;
      }
    }
    std::cout << std::endl; 
  }

  // map MS2 features
  SpecFeaturePtrVec ms2_features;
  env_coll_assign::assignEnvColls(frac_features, env_coll_list, ms2_header_ptr_2d,
                                  ms2_features, topfd_para_ptr->getEcscoreCutoff()); 
  std::cout << "Number of proteoform features: " << env_coll_list.size() << std::endl;
  /// output files
  if (topfd_para_ptr->isOutputCsvFeatureFile()) {
    std::string feat_file_name = output_base_name + "_ms1.csv";
    ecscore_writer::writeScores(feat_file_name, ecscore_list);
    std::string batmass_file_name = output_base_name + "_" + "frac.mzrt.csv";
    frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);
  }

  std::string output_file_name = output_base_name + "_" + "feature.xml";
  frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);
  output_file_name = output_base_name + "_" + "ms1.frac_feature";
  frac_feature_writer::writeFeatures(output_file_name, frac_features);

  SampleFeaturePtrVec sample_features;
  feature_util::getSampleFeatures(sample_features, frac_features, ms2_features);
  std::string sample_feature_file_name = output_base_name + "_"  + "ms1.feature";
  sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);
  output_file_name = output_base_name + "_"  + "ms2.feature";
  spec_feature_writer::writeFeatures(output_file_name, ms2_features);
}

void processMs2(TopfdParaPtr topfd_para_ptr, MsAlignWriterPtr writer_ptr, 
                FracFeaturePtrVec all_ms1_feat_list, 
                double mz_bgn, double mz_end, std::set<int> spec_id_set) {
  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(topfd_para_ptr->getFracId(), 
                                                                topfd_para_ptr->getMzmlFileName(),
                                                                topfd_para_ptr);
  score_para_ptr->min_match_peak_ = 1;
  score_para_ptr->min_scan_num_ = 1;
  // read deconvoluted MS2 peaks
  LOG_ERROR("Read ms2 file started");
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();
  std::string ms2_file_name = output_base_name + "_pre_ms2.msalign";
  DeconvMsPtrVec all_deconv_ms2_ptr_vec;
  msalign_reader_util::readAllSpectra(ms2_file_name, all_deconv_ms2_ptr_vec);
  DeconvMsPtrVec deconv_ms2_ptr_vec;
  LOG_ERROR("Total ms2 scan number " << all_deconv_ms2_ptr_vec.size());
  std::set<int> prec_feat_ids;
  for (size_t i = 0; i < all_deconv_ms2_ptr_vec.size(); i++) {
    int id = all_deconv_ms2_ptr_vec[i]->getMsHeaderPtr()->getSpecId();
    if (spec_id_set.find(id) != spec_id_set.end()) {
      deconv_ms2_ptr_vec.push_back(all_deconv_ms2_ptr_vec[i]);
      PrecursorPtrVec precursor_list = all_deconv_ms2_ptr_vec[i]->getMsHeaderPtr()->getPrecPtrVec();       
      for (size_t j = 0; j < precursor_list.size(); j++) {
        prec_feat_ids.insert(precursor_list[j]->getFeatureId());
      }
    }
  }
  FracFeaturePtrVec ms1_feat_list; 
  for (int feat_id : prec_feat_ids) {
    ms1_feat_list.push_back(all_ms1_feat_list[feat_id]);
  }
  LOG_ERROR("Selected ms2 scan number " << deconv_ms2_ptr_vec.size());
  LOG_ERROR("Selected ms1 feature number " << ms1_feat_list.size()); 

  // read ms2 raw peaks 
  PeakPtrVec2D ms2_mzml_peaks;
  MzmlMsGroupReaderPtr mzml_reader_ptr = 
    std::make_shared<MzmlMsGroupReader>(topfd_para_ptr->getMzmlFileName(), 
                                        topfd_para_ptr->getPrecWindowWidth(),
                                        topfd_para_ptr->getActivation(),
                                        topfd_para_ptr->getFracId(),
                                        topfd_para_ptr->isFaims(), 
                                        topfd_para_ptr->getFaimsVoltage(), 
                                        topfd_para_ptr->isMissingLevelOne());
  mzml_reader_ptr->getMs2Map(ms2_mzml_peaks, spec_id_set); 
  mzml_reader_ptr = nullptr;

  //Prepare seed envelopes
  SeedEnvPtrVec seed_ptrs;
  for (size_t i = 0; i < deconv_ms2_ptr_vec.size(); i++) {
    DeconvMsPtr ms2_data = deconv_ms2_ptr_vec[i];
    DeconvPeakPtrVec peaks = ms2_data->getPeakPtrVec();
    SeedEnvPtrVec one_spec_seed_ptrs;
    for (auto &peak: peaks) {
      SeedEnvPtr seed_ptr_1 = std::make_shared<SeedEnv>(i, peak);
      seed_ptrs.push_back(seed_ptr_1);
    }
  }

  std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnv::cmpSeedInteDec);
  // write_out_files::write_seed_envelopes(seed_envs, "envs.csv");

  double sn_ratio = topfd_para_ptr->getMsTwoSnRatio();
  bool single_scan_noise = true; 
  /// Prepare data -- Peak Matrix
  MsMapPtr matrix_ptr = std::make_shared<MsMap>(ms2_mzml_peaks, deconv_ms2_ptr_vec,
                                                score_para_ptr->bin_size_,
                                                sn_ratio, single_scan_noise);

  /// Extract Fetures
  LOG_ERROR("Number of seed envelopes: " << seed_ptrs.size());
  int seed_num = seed_ptrs.size();
  EnvCollPtrVec env_coll_list;
  int feat_id = 0;
  double perc = 0;
  ECScorePtrVec ecscore_list;
  FracFeaturePtrVec ms2_feat_list;
  for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
    int count = seed_env_idx + 1;
    if (count % 100 == 0 || count == seed_num) {
      perc = static_cast<int>(count * 100 / seed_num);
      std::cout << "\r" << "Processing feature " << count << " ...       " << perc << "\% finished." << std::flush;
    }
    SeedEnvPtr ori_seed_ptr = seed_ptrs[seed_env_idx];
    SeedEnvPtr seed_ptr = seed_env_util::preprocessSeedEnvPtr(ori_seed_ptr, matrix_ptr,  
                                                              score_para_ptr, sn_ratio); 
    
    EnvCollPtr env_coll_ptr;
    if (seed_ptr != nullptr) { 
      env_coll_ptr = env_coll_util::findEnvColl(matrix_ptr, seed_ptr,
                                                score_para_ptr, sn_ratio); 
    }
    if (seed_ptr == nullptr || env_coll_ptr == nullptr) {
      if (env_coll_util::checkSeedExistingFeatures(matrix_ptr, ori_seed_ptr,
                                                   env_coll_list, score_para_ptr)) {
        continue;
      }
      FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(feat_id,
                                                                   deconv_ms2_ptr_vec, 
                                                                   score_para_ptr->frac_id_,
                                                                   score_para_ptr->file_name_,
                                                                   ori_seed_ptr, matrix_ptr);
      ms2_feat_list.push_back(frac_feat_ptr);
      feat_id++;
    }
    else {
      if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                               env_coll_list, score_para_ptr)) {
        env_coll_ptr->removePeakData(matrix_ptr);
        continue;
      }
      env_coll_ptr->refineMonoMass();
      ECScorePtr ecscore_ptr = std::make_shared<ECScore>(env_coll_ptr, matrix_ptr,
                                                         feat_id, sn_ratio); 
      ecscore_list.push_back(ecscore_ptr);
      env_coll_ptr->setEcscore(ecscore_ptr->getScore());
      env_coll_ptr->removePeakData(matrix_ptr);
      env_coll_list.push_back(env_coll_ptr);
      FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(feat_id,
                                                                   deconv_ms2_ptr_vec, 
                                                                   score_para_ptr->frac_id_,
                                                                   score_para_ptr->file_name_,
                                                                   env_coll_ptr, matrix_ptr, sn_ratio);
      ms2_feat_list.push_back(frac_feat_ptr);
      feat_id++;
    }
  }
  std::cout << std::endl; 

  std::cout << "Number of fragment features: " << ms2_feat_list.size() << std::endl;
  std::string ms2_feature_file_name = output_base_name + "_" + str_util::toString(mz_bgn) + "_" + "ms2_feature.tsv";
  frac_feature_writer::writeFeatures(ms2_feature_file_name, ms2_feat_list);

  //get pseudo spectra
  int ms2_scan_num = 20;
  std::sort (ms1_feat_list.begin(), ms1_feat_list.end(), FracFeature::cmpInteDec);
  for (size_t i = 0; i < ms1_feat_list.size(); i++) {
    FracFeaturePtr ms1_feat = ms1_feat_list[i];
    int apex_ms1_scan = ms1_feat->getApexScan();
    int min_scan = apex_ms1_scan - ms2_scan_num;
    int max_scan = apex_ms1_scan + ms2_scan_num;
    FracFeaturePtrVec ms2_sele_list;
    FracFeaturePtrVec ms2_remain_feat_list;
    for (size_t j = 0; j < ms2_feat_list.size(); j++) {
      FracFeaturePtr ms2_feat = ms2_feat_list[j];
      if (ms2_feat->getMinMs1Id() == ms2_feat->getMaxMs1Id()) {
        if (ms2_feat->getApexScan() > apex_ms1_scan 
            && ms2_feat->getApexScan() <= max_scan 
            && ms2_feat->getMonoMass() < ms1_feat->getMonoMass()) {
          ms2_sele_list.push_back(ms2_feat);
        }
        else {
          ms2_remain_feat_list.push_back(ms2_feat);
        }
      }
      else {
        LOG_DEBUG("ms2 feat apex scan " << ms2_feat->getApexScan());
        if (ms2_feat->getApexScan() >= min_scan 
            && ms2_feat->getApexScan() <= max_scan 
            && ms2_feat->getMonoMass() < ms1_feat->getMonoMass()) {
          ms2_sele_list.push_back(ms2_feat);
        }
        else {
          ms2_remain_feat_list.push_back(ms2_feat);
        }
      }
    }
    if (ms2_sele_list.size() > 10) {
      ms2_feat_list = ms2_remain_feat_list;
      LOG_ERROR(i << " ms1 feature " << ms1_feat->getId() << " apex scan " << ms1_feat->getApexScan() <<  " frac number " << ms2_sele_list.size());
      DeconvMsPtr ms2_ptr;
      for (size_t j = 0; j < deconv_ms2_ptr_vec.size(); j++) {
        int ms2_scan = deconv_ms2_ptr_vec[j]->getMsHeaderPtr()->getFirstScanNum();
        if (ms2_scan > apex_ms1_scan && ms2_scan <= max_scan) {
          ms2_ptr = deconv_ms2_ptr_vec[j];
          break;
        }
      }
      if (ms2_ptr != nullptr) {
        LOG_ERROR("Ms 2 scan " << ms2_ptr->getMsHeaderPtr()->getFirstScanNum());
        MsHeaderPtr header_ptr = ms2_ptr->getMsHeaderPtr();
        header_ptr->setSpecId(topfd_para_ptr->getDiaSpecId());
        topfd_para_ptr->incDiaSpecId();
        double mono_mass = ms1_feat->getMonoMass();
        double prec_win_center = (header_ptr->getPrecWinBegin() + header_ptr->getPrecWinEnd())/2;
        int charge = std::round(mono_mass/prec_win_center);
        double mono_mz = peak_util::compMz(mono_mass, charge);
        PrecursorPtr prec_ptr = std::make_shared<Precursor>(0,
                                                            ms1_feat->getId(), 
                                                            mono_mz, 
                                                            charge, 
                                                            ms1_feat->getIntensity());
        PrecursorPtrVec prec_ptr_vec;
        prec_ptr_vec.push_back(prec_ptr);
        header_ptr->setPrecPtrVec(prec_ptr_vec);

        DeconvPeakPtrVec peak_list;
        int sp_id = header_ptr->getSpecId();
        for (size_t j = 0; j < ms2_sele_list.size(); j++) {
          FracFeaturePtr ms2_feat = ms2_sele_list[j];
          DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, j, 
                                                                ms2_feat->getMonoMass(),
                                                                ms2_feat->getIntensity(),
                                                                ms2_feat->getRepCharge(),
                                                                ms2_feat->getEcScore());
          peak_list.push_back(peak_ptr);
        }
        ms2_ptr->setPeakPtrVec(peak_list);
        writer_ptr->writeMs(ms2_ptr);
      }
    }
  }
}

}  // namespace

}  // namespace toppic
