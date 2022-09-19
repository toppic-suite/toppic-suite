//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <ctime>

#include "common/util/file_util.hpp"
#include "seq/fasta_util.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/env/env_para.hpp"
#include "ms/feature/frac_feature.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "topfd/feature_detect/feature_para.hpp"
#include "topfd/feature_detect/feature_detect.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/env_collection/env_coll_util.hpp"
#include "topfd/feature_detect/feature/feature.hpp"
#include "topfd/feature_detect/util/write_feature.hpp"
#include "topfd/feature_detect/envelope/evaluate_envelope.hpp"
#include "topfd/feature_detect/test_output_functions/write_out_files.hpp"

namespace toppic {

namespace feature_detect {

void process(int frac_id, const std::string &sp_file_name,
             bool missing_level_one, const std::string &resource_dir, const std::string &activation,
             bool isFaims, const std::vector<std::pair<double, int>> voltage_vec) {
  std::clock_t start;
  double duration;

  start = std::clock();

  //logger::setLogLevel(2);
  FeatureParaPtr para_ptr
      = std::make_shared<FeaturePara>(frac_id, sp_file_name, resource_dir);
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
  std::string base_name = file_util::basename(sp_file_name);

  ///// params
  double bin_size = 0.1;
  double neighbor_mass_tole = 0.01;
  double mass_tole = 0.008;
  double corr_tole = 0.05;
  int max_miss_env = 2;
  int max_miss_charge = 2;
  int max_miss_peak = 2;
  int para_max_charge = 30;
  double ratio_multi = 2.0;
  int match_peak_tole = 3;

  // read ms1 deconvoluted spectra
  if (!missing_level_one) {
    std::string file_num;
    //for (size_t i = 0; i < voltage_vec.size(); i++) {
    for (size_t i = 0; i <1; i++) {
      if (isFaims) {file_num = str_util::toString(i) + "_";} // if FAIME Data

      FracFeaturePtrVec frac_features;
      std::string file_name = resource_dir
                              + file_util::getFileSeparator() + "envcnn_models"
                              + file_util::getFileSeparator() + "envcnn_2_block_model.json";
      fdeep::model model = fdeep::load_model(file_name, true, fdeep::dev_null_logger);

      /// Read msalign file and get the seed envelopes.
      DeconvMsPtrVec ms1_ptr_vec;
      std::string ms1_file_name = base_name + "_" + file_num + "ms1.msalign";
      SimpleMsAlignReader::readMsOneSpectra(ms1_file_name, ms1_ptr_vec);
      std::cout << "Processed msalign file." << std::endl;

      /// Read mzml data
      //double cur_voltage = voltage_vec[i].first;//if this is -1, it is non-FAIME data
      double cur_voltage = -1;
      PeakPtrVec2D raw_peaks;
      RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(sp_file_name, activation, env_para_ptr->prec_deconv_interval_);
      raw_reader_ptr->getMs1Peaks(raw_peaks, cur_voltage);
      raw_reader_ptr = nullptr;
      std::cout << "Processed mzML file." << std::endl;

      /// Prepare data -- seed envelopes
      std::vector<SeedEnvelope> seed_envs;
      for (auto &ms1_data : ms1_ptr_vec) {
        std::vector<DeconvPeakPtr> peaks =  ms1_data->getPeakPtrVec();
        for (auto &peak : peaks)
          seed_envs.push_back(SeedEnvelope(peak));
      }
      std::sort(seed_envs.begin(), seed_envs.end(), SeedEnvelope::cmpInteDec);
//      write_out_files::write_seed_envelopes(seed_envs,  base_name + "_" + "envs.txt");

      /// Prepare data -- Peak Matrix
      PeakMatrix peak_matrix = PeakMatrix(raw_peaks, ms1_ptr_vec, bin_size, env_para_ptr->ms_one_sn_ratio_);
//      write_out_files::write_peak_matrix(peak_matrix, base_name + "_" + "matrix.txt");
      peak_matrix.find_remove_non_neighbors(neighbor_mass_tole);
//      write_out_files::write_peak_matrix(peak_matrix, base_name + "_" + "matrix_labeled.txt");

      std::cout << "Number of seed envelopes: " << seed_envs.size() << std::endl;

      /// Extract Fetures
      int seed_num = seed_envs.size();
      int env_coll_num = 0;
      std::vector<EnvCollection> env_coll_list;
      std::vector<Feature> features;
      for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++){
        if (seed_env_idx % 500 == 0)
          std::cout << "Processing peak " << seed_env_idx << " and Features found " << env_coll_num << std::endl;
        SeedEnvelope env = seed_envs[seed_env_idx];
        bool valid = false;
        valid = evaluate_envelope::preprocess_env(peak_matrix, env, mass_tole, corr_tole, valid);
        if (!valid) continue;
        EnvCollection env_coll = env_coll_util::find_env_collection(peak_matrix, env, mass_tole, max_miss_env, max_miss_charge, max_miss_peak,
                                                                    para_max_charge, ratio_multi, match_peak_tole, env_para_ptr->ms_one_sn_ratio_,
                                                                    base_name + "_" + str_util::toString(env_coll_num) +"_details.txt");
        if (!env_coll.isEmpty()) {
          if (env_coll_util::check_in_existing_features(peak_matrix, env_coll, env_coll_list, 10E-6, 0.8))
            continue;
          env_coll.refine_mono_mass();
          env_coll_list.push_back(env_coll);
          features.push_back(Feature(env_coll, peak_matrix, model, env_coll_num, env_para_ptr->ms_one_sn_ratio_));
          env_coll.remove_peak_data(peak_matrix);
          env_coll_num = env_coll_num + 1;
        }
      }
      std::cout << "Number of Envelope Collections: " << features.size() << std::endl;
      file_name = base_name + "_" + file_num + "ms1.feature";
      write_feature::writeFeatures(file_name, features);
      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      std::cout<<"Duration (seconds): "<< duration <<'\n';
    }
  }
}
}  // namespace
}  // namespace toppic
