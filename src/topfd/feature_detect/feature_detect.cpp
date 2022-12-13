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

#include "topfd/feature_detect/feature_detect.hpp"
#include "ms/feature/sample_feature.hpp"
#include "feature_detect_old.hpp"

namespace toppic {
  namespace feature_detect {
    void process(int frac_id, const std::string &sp_file_name,
                 bool missing_level_one, const std::string &resource_dir, const std::string &activation,
                 bool isFaims, const std::vector<std::pair<double, int>> voltage_vec, double score_cutoff) {
      std::clock_t start;
      start = std::clock();
      //logger::setLogLevel(2);
      FeatureParaPtr para_ptr = std::make_shared<FeaturePara>(frac_id, sp_file_name, resource_dir);
      EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
      std::string base_name = file_util::basename(sp_file_name);

      // read ms1 deconvoluted spectra
      if (!missing_level_one) {
        std::string file_num;
        for (size_t i = 0; i < 1; i++) { //////////////////////////////////////////////////////
//        for (size_t i = 0; i < voltage_vec.size(); i++) {
          if (isFaims) { file_num = str_util::toString(i) + "_"; } // if FAIME Data
          FracFeaturePtrVec frac_features;

          std::string ms1_file_name = base_name + "_" + file_num + "ms1.msalign";
          std::string ms2_file_name = base_name + "_" + file_num + "ms2.msalign";

          // EnvCNN
          std::string file_name = resource_dir
                                  + file_util::getFileSeparator() + "envcnn_models"
                                  + file_util::getFileSeparator() + "envcnn_2_block_model.json";
          fdeep::model model = fdeep::load_model(file_name, true, fdeep::dev_null_logger);

          // ECScore
          file_name = resource_dir
                      + file_util::getFileSeparator() + "envcnn_models"
                      + file_util::getFileSeparator() + "escore.json";
          fdeep::model model_escore = fdeep::load_model(file_name, true, fdeep::dev_null_logger);

          /// Read msalign file and get the seed envelopes.
          DeconvMsPtrVec ms1_ptr_vec;
          SimpleMsAlignReader::readMsOneSpectra(ms1_file_name, ms1_ptr_vec);
          std::cout << "Processed msalign file." << std::endl;

          /// Read mzml data
//          double cur_voltage = voltage_vec[i].first;//if this is -1, it is non-FAIME data
          double cur_voltage = 0; //////////////////////////////////////////////////////
          PeakPtrVec2D raw_peaks;
          std::vector<double> prec_spectrum_Base_mono_mz;
          RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(sp_file_name, activation,
                                                                        env_para_ptr->prec_deconv_interval_);
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

          /// Prepare data -- Peak Matrix
          PeakMatrix peak_matrix = PeakMatrix(raw_peaks, ms1_ptr_vec, para_ptr->bin_size_,
                                              env_para_ptr->ms_one_sn_ratio_);
          peak_matrix.find_remove_non_neighbors(para_ptr->neighbor_mass_tole_);

          /// Extract Fetures
          std::cout << "Number of seed envelopes: " << seed_envs.size() << std::endl;
          int seed_num = seed_envs.size();
          int env_coll_num = 0;
          std::vector<EnvCollection> env_coll_list;
          std::vector<Feature> features;
          for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
            if (seed_env_idx % 10000 == 0)
              std::cout << "\r" << "Processing peak " << seed_env_idx << " and Features found " << env_coll_num
                        << std::flush;
            SeedEnvelope env = seed_envs[seed_env_idx];
            bool valid = false;
            valid = evaluate_envelope::preprocess_env(peak_matrix, env, para_ptr, env_para_ptr->ms_one_sn_ratio_);
            if (!valid) continue;
            EnvCollection env_coll = env_coll_util::find_env_collection(peak_matrix, env, para_ptr,
                                                                        env_para_ptr->ms_one_sn_ratio_);
            if (!env_coll.isEmpty()) {
              if (env_coll_util::check_in_existing_features(peak_matrix, env_coll, env_coll_list, para_ptr))
                continue;
              env_coll.refine_mono_mass();

              Feature feature = Feature(env_coll, peak_matrix, model, model_escore, env_coll_num,
                                        env_para_ptr->ms_one_sn_ratio_);
              if (feature.getScore() < score_cutoff) continue;
              features.push_back(feature);
              env_coll.setEcscore(feature.getScore());
              env_coll.remove_peak_data(peak_matrix);
              env_coll_list.push_back(env_coll);
              FracFeaturePtr feature_ptr = Feature::getFeature(env_coll_num, ms1_ptr_vec, para_ptr->frac_id_,
                                                               para_ptr->file_name_, env_coll, peak_matrix,
                                                               env_para_ptr->ms_one_sn_ratio_);
              feature_ptr->setPromexScore(feature.getScore());
              frac_features.push_back(feature_ptr);
              env_coll_num = env_coll_num + 1;
            }
          }
          // map MS2 features
          SpecFeaturePtrVec ms2_features;
          Feature::assign_features(ms1_ptr_vec, ms2_file_name, frac_features, env_coll_list, features, ms2_features,
                                   prec_spectrum_Base_mono_mz, peak_matrix, model, model_escore, para_ptr,
                                   env_para_ptr, score_cutoff);
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

          double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
          std::cout << "Feature Extraction duration (minutes): " << duration / 60.0 << '\n';
        }
      }
    }
  }  // namespace
}  // namespace toppic
