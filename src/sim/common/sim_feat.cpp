//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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
#include <cstdlib> 
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
//#include "sim/png/pngwriter.hpp"

namespace toppic {

namespace sim_feat {

void writepng(MsMapPtr map_ptr, EnvSetPtr env_set_ptr, std::string png_file_name) {
  int scale =  100;
  int width = 512;
  int height = 512;
  int center = width/2;
  double mz_win_size = width/scale;

  // get the background row ranges
  int ref_spec_id = env_set_ptr->getSeedPtr()->getSpecId();
  int spec_range = env_set_ptr->getEndSpecId() - env_set_ptr->getStartSpecId() + 1;
  int direction = (rand() % 2 == 0) ? -1 : 1;
  int png_start_spec = ref_spec_id  + direction  * spec_range;
  int png_end_spec = png_start_spec + height - 1;
  if (direction == -1) {
    png_end_spec = ref_spec_id + direction * spec_range;
    png_start_spec = png_end_spec - height + 1;
  } 
  int png_ref_spec = png_start_spec + height/2;
  int spec_shift = png_ref_spec - ref_spec_id;

  if (png_start_spec < 0 || png_end_spec >= map_ptr->getRowNum()) {
    return;
  }

  // get the background mz range
  SeedEnvPtr seed_env = env_set_ptr->getSeedPtr(); 
  double ref_mz = seed_env->getReferMz(); 
  double left_mz = ref_mz - mz_win_size/2.0;
  double right_mz = ref_mz + mz_win_size/2.0;

  // get max intensity in m/z range
  double max_inte = 0;
  MsMapPeakPtr2D peak_2d_list = map_ptr->get2DPeaks();
  for (size_t i = 0; i < peak_2d_list.size(); i++) {
    for (size_t j = 0; j < peak_2d_list[i].size(); j++) {
      MsMapPeakPtr peak = peak_2d_list[i][j];
      double peak_pos = peak->getPosition();
      if (peak_pos >= left_mz && peak_pos < right_mz) {
        if (peak->getIntensity() > max_inte) {
          max_inte = peak->getIntensity();
        }
      }
    }
  }

  //pngwriter png(width, height, 1.0, (png_file_name +"_data.png").c_str());
  // add background
  for (size_t i = png_start_spec; i <= png_end_spec; i++) {
    for (size_t j = 0; j < peak_2d_list[i].size(); j++) {
      MsMapPeakPtr peak = peak_2d_list[i][j];
      double peak_pos = peak->getPosition();
      if (peak_pos >= left_mz && peak_pos < right_mz) {
        int x = (peak_pos - ref_mz)*scale + center;
        int y = i - png_start_spec;
        double inte = peak->getIntensity() / max_inte;
        //png.plot(x, y, inte, 0.0, 0.0);
      }
    }
  }

  MsMapEnvPtrVec ms_map_env_list = env_set_ptr->getMsMapEnvList();
  int start_spec_id = env_set_ptr->getStartSpecId();
  int end_spec_id = env_set_ptr->getEndSpecId();
  int peak_num = seed_env->getPeakPtrList().size();
  std::vector<int> y_start_list(peak_num, -1);
  std::vector<int> y_end_list (peak_num, -1); 
  for (int spec_id = start_spec_id; spec_id < end_spec_id; spec_id++) {
    MsMapEnvPtr env_ptr = ms_map_env_list[spec_id-start_spec_id];
    MsMapPeakPtrVec exp_peak_list = env_ptr->getMsMapPeakList();
    LOG_DEBUG("spec id " << spec_id << " peak num " << exp_peak_list.size());
    for (size_t i = 0; i < exp_peak_list.size(); i++) {
      MsMapPeakPtr peak = exp_peak_list[i];
      if (peak == nullptr) {
        continue;
      }
      int x = (peak->getPosition() - ref_mz)*scale + center;
      int y = spec_id + spec_shift - png_start_spec; 
      if (y_start_list[i] == -1) {
        y_start_list[i] = y;
        y_end_list[i] = y;
      }
      else {
        if (y < y_start_list[i]) {
          y_start_list[i] = y;
        }
        if (y > y_end_list[i]) {
          y_end_list[i] = y;
        }
      }
      double inte = peak->getIntensity() /max_inte;
      LOG_DEBUG("x: " << x << " pos: " << peak->getPosition() << " y: " << y << " peak intensity: " << peak->getIntensity() << " max inte " << max_inte);
      //png.plot(x, y, inte, 0.0, 0.0);
    }
  }
  //png.close();

  // update y_start_list and y_end_list
  int ref_idx = seed_env->getReferIdx();
  for (int i = 0; i < ref_idx; i++) {
    if (y_start_list[i] == -1) {
      continue;
    }
    if (y_start_list[i+1] == -1) {
      y_start_list[i+1] = y_start_list[i];
      y_end_list[i+1] = y_end_list[i];
    }
    else {
      if (y_start_list[i+1] > y_start_list[i]) {
        y_start_list[i+1] = y_start_list[i];
      }
      if (y_end_list[i+1] < y_end_list[i]) {
        y_end_list[i+1] = y_end_list[i];
      }
    }
  }
  for (int i = peak_num -1; i >ref_idx; i--) {
    if (y_start_list[i] == -1) {
      continue;
    }
    if (y_start_list[i-1] == -1) {
      y_start_list[i-1] = y_start_list[i];
      y_end_list[i-1] = y_end_list[i];
    }
    else {
      if (y_start_list[i-1] > y_start_list[i]) {
        y_start_list[i-1] = y_start_list[i];
      }
      if (y_end_list[i-1] < y_end_list[i]) {
        y_end_list[i-1] = y_end_list[i];
      }
    }
  }

  //pngwriter mask_png(width, height, 0.0, (png_file_name + "_mask.png").c_str());
  // plot mask 
  std::vector<int> pos_list;
  for (int i = 0; i < peak_num; i++) {
    pos_list.push_back((seed_env->getMz(i) - ref_mz)*scale + center);
    LOG_DEBUG(i << " start " << y_start_list[i] << " end " << y_end_list[i]); 
  }
  for (int i = 0; i < peak_num -1; i++) {
    if (y_start_list[i] == -1 or y_start_list[i+1] == -1) { 
      continue;
    }
    int left_pos = pos_list[i];
    int right_pos = pos_list[i+1];
    for (int j = left_pos; j <= right_pos; j++) {
      int start = y_start_list[i] + (y_start_list[i+1] - y_start_list[i]) * (j - left_pos) / (right_pos - left_pos);
      int end = y_end_list[i] + (y_end_list[i+1] - y_end_list[i]) * (j - left_pos) / (right_pos - left_pos);
      for (int k = start; k <= end; k++) {
        //mask_png.plot(j, k, 1.0, 1.0, 1.0);
      }
    }
  }
  /*
  for (int i = 0; i < peak_num; i++) {
    if (y_start_list[i] == -1 or y_start_list[i+1] == -1) { 
      continue;
    }
    int pos = pos_list[i];
    int start = y_start_list[i]; 
    int end = y_end_list[i];
    for (int k = start; k <= end; k++) {
      //mask_png.plot(pos, k, 0.0, 0.0, 1.0);
    }
  }
  */
  //mask_png.close();
}

void processMs1(TopfdParaPtr topfd_para_ptr) {
  if (topfd_para_ptr->isMissingLevelOne()) {
    return;
  }
  
  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(topfd_para_ptr->getFracId(), 
                                                                topfd_para_ptr->getMzmlFileName(),
                                                                topfd_para_ptr->getMaxCharge(),
                                                                topfd_para_ptr->getMs1MinScanNum());
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
  for (auto &ms1_data: deconv_ms1_ptr_vec) {
    DeconvPeakPtrVec peaks = ms1_data->getPeakPtrVec();
    SeedEnvPtrVec one_spec_seed_ptrs;
    for (auto &peak: peaks) {
      SeedEnvPtr seed_ptr_1 = std::make_shared<SeedEnv>(peak);
      seed_ptrs.push_back(seed_ptr_1);
      SeedEnvPtr seed_ptr_2 = std::make_shared<SeedEnv>(peak);
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
  seed_num = 1;
  EnvCollPtrVec env_coll_list;
  ECScorePtrVec ecscore_list;
  for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
    int count = seed_env_idx + 1;
    if (count % 100 == 0 || count == seed_num) {
      double perc = static_cast<int>(count * 100 / seed_num);
      std::cout << "\r" << "Processing feature " << count << " ...       " << perc << "\% finished." << std::flush;
    }
    SeedEnvPtr seed_ptr = seed_ptrs[seed_env_idx];
    seed_ptr = seed_env_util::preprocessSeedEnvPtr(seed_ptr, matrix_ptr,  
                                                   score_para_ptr, sn_ratio); 
    if (seed_ptr == nullptr) continue;
    EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(matrix_ptr, seed_ptr,
                                                         score_para_ptr, sn_ratio, 
                                                         topfd_para_ptr->getSplitIntensityRatio());
    if (env_coll_ptr == nullptr) continue;
    if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                             env_coll_list, score_para_ptr, sn_ratio)) {
      env_coll_ptr->removePeakData(matrix_ptr);
      continue;
    }
    env_coll_ptr->refineMonoMass();
    ECScorePtr ecscore_ptr = std::make_shared<ECScore>(env_coll_ptr, matrix_ptr,
                                                       sn_ratio); 
    if (ecscore_ptr->getScore() < topfd_para_ptr->getMs1EcscoreCutoff()) {
      continue;
    }
    env_coll_ptr->setEcscore(ecscore_ptr->getScore());
    env_coll_ptr->removePeakData(matrix_ptr);
    env_coll_list.push_back(env_coll_ptr);
    ecscore_list.push_back(ecscore_ptr);
  }
  std::cout << std::endl; 

  /*
  if (topfd_para_ptr->isSearchPrecWindow()) {
    // add ms1 feature based on precursor windows
    // set min match envelope to 1 to accept single scan features
    score_para_ptr->min_scan_num_ = 1;
    score_para_ptr->min_match_peak_ = 1;
    sn_ratio = 0;
    matrix_ptr->reconstruct(sn_ratio, single_scan_noise); 
    int ms1_spec_num = deconv_ms1_ptr_vec.size();
    for (std::size_t ms1_idx = 0; ms1_idx < deconv_ms1_ptr_vec.size(); ms1_idx++) {
      double perc = static_cast<int>((ms1_idx + 1)* 100 / ms1_spec_num);
      int scan = deconv_ms1_ptr_vec[ms1_idx]->getMsHeaderPtr()->getFirstScanNum();
      std::cout << "\r" << "Additional feature search MS1 spectrum scan " 
        << scan << " ...       " << perc << "\% finished." << std::flush;
      for (std::size_t i = 0; i < ms2_header_ptr_2d[ms1_idx].size(); i++) {
        MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_idx][i];
        if (env_coll_assign::checkEnvColl(header_ptr, env_coll_list)) {
          continue;
        }
        double prec_win_begin = header_ptr->getPrecWinBegin();
        double prec_win_end = header_ptr->getPrecWinEnd();
        SeedEnvPtrVec seed_ptr_list = seed_ptr_2d[ms1_idx];
        SeedEnvPtrVec selected_seed_list;
        LOG_DEBUG("ms1 id  " << ms1_idx << " seed number " << seed_ptr_list.size());
        for (std::size_t i = 0; i < seed_ptr_list.size(); i++) {
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
                                                                          sn_ratio, topfd_para_ptr->getSplitIntensityRatio());
        if (env_coll_ptr == nullptr) continue;
        if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                                 env_coll_list, score_para_ptr, sn_ratio)) {
          env_coll_ptr->removePeakData(matrix_ptr);
          continue;
        }
        env_coll_ptr->refineMonoMass();
        ECScorePtr ecscore_ptr = std::make_shared<ECScore>(env_coll_ptr, matrix_ptr,
                                                           sn_ratio); 
        if (ecscore_ptr->getScore() < 0 || std::isnan(ecscore_ptr->getScore())) {
          continue;
        }
        env_coll_ptr->setEcscore(ecscore_ptr->getScore());
        env_coll_ptr->removePeakData(matrix_ptr);
        env_coll_list.push_back(env_coll_ptr);
        ecscore_list.push_back(ecscore_ptr);
      }
    }
    std::cout << std::endl; 
  }
  */

  /*
  FracFeaturePtrVec frac_features;
  for (std::size_t i = 0; i < env_coll_list.size(); i++) {
    EnvCollPtr env_coll_ptr = env_coll_list[i];
    FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(i, deconv_ms1_ptr_vec, 
                                                                 score_para_ptr->frac_id_,
                                                                 score_para_ptr->file_name_,
                                                                 env_coll_ptr, matrix_ptr, sn_ratio);
    frac_features.push_back(frac_feat_ptr);
  }

  /// output files
  if (topfd_para_ptr->isOutputCsvFeatureFile()) {
    std::string feat_file_name = output_base_name + "_ms1.csv";
    ecscore_writer::writeScores(feat_file_name, ecscore_list);
    std::string batmass_file_name = output_base_name + "_" + "frac_ms1.mzrt.csv";
    int num_spec = ms1_mzml_peaks.size();
    frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features, num_spec);
  }

  std::string frac_feat_xml_file_name = output_base_name + "_feature.xml";
  frac_feature_writer::writeXmlFeatures(frac_feat_xml_file_name, frac_features);
  std::string frac_feat_file_name = output_base_name + "_" + "ms1.feature";
  frac_feature_writer::writeFeatures(frac_feat_file_name, frac_features);
  */

  PeakPtrVec2D ms1_mzml_peaks_new;
  MsHeaderPtr2D ms2_header_ptr_2d_new;
  MzmlMsGroupReaderPtr mzml_reader_ptr_new = 
    std::make_shared<MzmlMsGroupReader>(topfd_para_ptr->getMzmlFileName(), 
                                        topfd_para_ptr->getPrecWindowWidth(),
                                        topfd_para_ptr->getActivation(),
                                        topfd_para_ptr->getFracId(),
                                        topfd_para_ptr->isFaims(), 
                                        topfd_para_ptr->getFaimsVoltage(), 
                                        topfd_para_ptr->isMissingLevelOne());
  mzml_reader_ptr_new->getMs1Map(ms1_mzml_peaks_new, ms2_header_ptr_2d_new); 
  mzml_reader_ptr_new = nullptr;

  sn_ratio = 0;
  MsMapPtr sim_matrix_ptr = std::make_shared<MsMap>(ms1_mzml_peaks_new, deconv_ms1_ptr_vec,
                                                    score_para_ptr->bin_size_,
                                                    sn_ratio, single_scan_noise);

  //for (std::size_t i = 0; i < env_coll_list.size(); i++) {
  for (std::size_t i = 0; i < 1; i++) {
    EnvSetPtrVec env_set_list = env_coll_list[i]->getEnvSetList();
    for (std::size_t j = 0; j < env_set_list.size(); j++) {
    //for (std::size_t j = 0; j < 1; j++) {
      EnvSetPtr env_set_ptr = env_set_list[j];
      std::cout << "feature " << i << " charge " << env_set_ptr->getCharge() << " mono mass " << env_set_ptr->getMonoMass() << std::endl;
      if (env_set_ptr != nullptr) {
        std::string png_file_name = topfd_para_ptr->getMzmlFileName()+ "_feature_" + std::to_string(i) + "_charge_" + std::to_string(env_set_ptr->getCharge());
        writepng(sim_matrix_ptr, env_set_ptr, png_file_name); 
      }
    }
  }

}

}  // namespace

}  // namespace toppic
