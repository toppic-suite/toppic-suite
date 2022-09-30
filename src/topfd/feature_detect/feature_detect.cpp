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

#include "topfd/feature_detect/feature_detect.hpp"
#include "topfd/feature_detect/feature_detect_old.hpp"
#include "common/util/file_util.hpp"
#include "seq/fasta_util.hpp"
#include "ms/spec/peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/env/env_para.hpp"
#include "ms/feature/frac_feature.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "topfd/feature_detect/feature_para.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/env_collection/env_coll_util.hpp"
#include "topfd/feature_detect/feature/feature.hpp"
#include "topfd/feature_detect/util/write_feature.hpp"
#include "topfd/feature_detect/envelope/evaluate_envelope.hpp"
#include "topfd/feature_detect/test_output_functions/write_out_files.hpp"
#include "ms/feature/spec_feature.hpp"

namespace toppic {

namespace feature_detect {

  FracFeaturePtr getFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, FeatureParaPtr para_ptr,
                            EnvCollection& env_coll, PeakMatrix &peak_matrix) {
    double snr = 3.0;
    double noise_inte = peak_matrix.get_min_inte();
    spec_list spectra_list = peak_matrix.get_spectra_list();
    int ms1_id_begin = env_coll.getStartSpecId();
    int ms1_id_end = env_coll.getEndSpecId();
    double feat_inte = env_coll.get_intensity(snr, peak_matrix.get_min_inte());
    double feat_mass = env_coll.getMass();
    int min_charge = env_coll.getMinCharge();
    int max_charge = env_coll.getMaxCharge();
    double time_apex = env_coll.get_apex_elution_time(spectra_list);
    double ms1_time_begin = env_coll.get_min_elution_time(spectra_list);
    double ms1_time_end = env_coll.get_max_elution_time(spectra_list);
    int ms1_scan_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getFirstScanNum();
    int ms1_scan_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getFirstScanNum();
    FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id,
                                                               para_ptr->frac_id_,
                                                               para_ptr->file_name_,
                                                               feat_mass,
                                                               feat_inte,
                                                               ms1_id_begin,
                                                               ms1_id_end,
                                                               ms1_time_begin,
                                                               ms1_time_end,
                                                               ms1_scan_begin,
                                                               ms1_scan_end,
                                                               min_charge,
                                                               max_charge,
                                                               0,
                                                               time_apex);
    SingleChargeFeaturePtrVec single_features;
    for (EnvSet &es: env_coll.getEnvSetList()) {
        int id_begin = es.getStartSpecId();
        int id_end = es.getEndSpecId();
        double time_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getRetentionTime();
        double time_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getRetentionTime();
        int scan_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getFirstScanNum();
        int scan_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getFirstScanNum();
        double inte = es.comp_intensity(snr, noise_inte);
        int charge = es.getCharge();
//        std::vector<std::vector<double>> map = es.get_map(snr, noise_inte);
//        int peak_num = component_score::get_num_theo_peaks(map);
        SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(charge,
                                                                                      time_begin, time_end,
                                                                                      scan_begin, scan_end,
                                                                                      inte, 0);
        single_features.push_back(single_feature);
    }
    feature_ptr->setSingleFeatures(single_features);
    return feature_ptr;
  }

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

      file_name = resource_dir
                              + file_util::getFileSeparator() + "envcnn_models"
                              + file_util::getFileSeparator() + "escore.json";
      fdeep::model model_escore = fdeep::load_model(file_name, true, fdeep::dev_null_logger);

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
          Feature feature = Feature(env_coll, peak_matrix, model, model_escore, env_coll_num, env_para_ptr->ms_one_sn_ratio_);
		  env_coll.remove_peak_data(peak_matrix);
          if (feature.getScore() < 0.5)
            continue;
          features.push_back(feature);
		  env_coll_list.push_back(env_coll);
          FracFeaturePtr feature_ptr = getFeature(env_coll_num, ms1_ptr_vec, para_ptr, env_coll, peak_matrix);
          feature_ptr->setPromexScore(feature.getScore());
          frac_features.push_back(feature_ptr);
          env_coll_num = env_coll_num + 1;
        }
      }

      std::string ms2_file_name = base_name + "_" + file_num + "ms2.msalign";
      MsHeaderPtrVec header_ptr_vec;
      SpecFeaturePtrVec ms2_features;
      feature_detect_old::readHeaders(ms2_file_name, header_ptr_vec);
      feature_detect_old::getMs2Features(ms1_ptr_vec, header_ptr_vec, frac_features, para_ptr, ms2_features);
      SampleFeaturePtrVec sample_features;
      feature_detect_old::getSampleFeatures(sample_features, frac_features, ms2_features);

      std::cout << "Number of Envelope Collections: " << features.size() << std::endl;
      file_name = base_name + "_" + file_num + "ms1.csv";
      write_feature::writeFeatures(file_name, features);


      std::string output_file_name = base_name + "_" + file_num + "feature.xml";
      frac_feature_writer::writeXmlFeatures(output_file_name, frac_features);
      std::string batmass_file_name = base_name + "_" + file_num + "frac.mzrt.csv";
      frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features);
      std::string sample_feature_file_name = base_name + "_" + file_num + "ms1.feature";
      sample_feature_writer::writeFeatures(sample_feature_file_name, sample_features);

      output_file_name = base_name + "_" + file_num + "ms2.feature";
      spec_feature_writer::writeFeatures(output_file_name, ms2_features);

      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      std::cout<<"Duration (seconds): "<< duration <<'\n';
    }
  }
}
}  // namespace
}  // namespace toppic
