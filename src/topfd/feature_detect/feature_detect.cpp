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

namespace toppic {
  namespace feature_detect {
    void process_ms2(int frac_id, const std::string &sp_file_name, const std::string &resource_dir, const std::string &activation, double score_cutoff) {
      std::clock_t start;
      start = std::clock();
      //logger::setLogLevel(2);
      FeatureParaPtr para_ptr = std::make_shared<FeaturePara>(frac_id, sp_file_name, resource_dir);
      EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
      std::string base_name = file_util::basename(sp_file_name);

      FracFeaturePtrVec frac_features;

      std::string ms2_file_name = base_name + "_ms2.msalign";

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
      DeconvMsPtrVec ms2_ptr_vec;
      SimpleMsAlignReader::readMsTwoSpectra(ms2_file_name, ms2_ptr_vec);
      std::cout << "Processed msalign file." << std::endl;

      /// Prepare data -- seed envelopes
      std::vector<SeedEnvelope> seed_envs;
      for (auto &ms2_data: ms2_ptr_vec) {
        std::vector<DeconvPeakPtr> peaks = ms2_data->getPeakPtrVec();
        for (auto &peak: peaks)
          seed_envs.push_back(SeedEnvelope(peak));
      }
      std::sort(seed_envs.begin(), seed_envs.end(), SeedEnvelope::cmpInteDec);
      write_out_files::write_seed_envelopes(seed_envs, "envs_828_832.txt");

      /// Read mzml data
      PeakPtrVec2D raw_peaks;
      RawMsReaderPtr raw_reader_ptr = std::make_shared<RawMsReader>(sp_file_name, activation,
                                                                    env_para_ptr->prec_deconv_interval_);
      raw_reader_ptr->getMs2Peaks(raw_peaks);
      raw_reader_ptr = nullptr;
      std::cout << "Processed mzML file." << std::endl;

      /// Prepare data -- Peak Matrix
      PeakMatrix peak_matrix = PeakMatrix(raw_peaks, ms2_ptr_vec, para_ptr->bin_size_,
                                          env_para_ptr->ms_two_sn_ratio_);
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
        bool valid = true;
        valid = evaluate_envelope::preprocess_env(peak_matrix, env, para_ptr, env_para_ptr->ms_two_sn_ratio_);
        if (!valid) continue;
        EnvCollection env_coll = env_coll_util::find_env_collection(peak_matrix, env, para_ptr,
                                                                    env_para_ptr->ms_two_sn_ratio_);
        if (!env_coll.isEmpty()) {
          if (env_coll_util::check_in_existing_features(peak_matrix, env_coll, env_coll_list, para_ptr))
            continue;
          env_coll.refine_mono_mass();

          Feature feature = Feature(env_coll, peak_matrix, model, model_escore, env_coll_num,
                                    env_para_ptr->ms_two_sn_ratio_);
          if (feature.getScore() < score_cutoff) continue;
          features.push_back(feature);
          env_coll.setEcscore(feature.getScore());
          env_coll.remove_peak_data(peak_matrix);
          env_coll_list.push_back(env_coll);
          FracFeaturePtr feature_ptr = Feature::getFeature(env_coll_num, ms2_ptr_vec, para_ptr->frac_id_,
                                                           para_ptr->file_name_, env_coll, peak_matrix,
                                                           env_para_ptr->ms_two_sn_ratio_);
          feature_ptr->setPromexScore(feature.getScore());
          frac_features.push_back(feature_ptr);
          env_coll_num = env_coll_num + 1;
        }
      }
      std::cout << std::endl << "Number of Envelope Collections: " << features.size() << std::endl;

      /// output temporary files
      file_name = base_name + "_ms2.csv";
      write_feature::writeFeatures(file_name, features);

      std::string batmass_file_name = base_name + "_frac_ms2.mzrt.csv";
      frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features, ms2_ptr_vec.size());

      /// Write pseudo spectra
      file_name = base_name + "_pseudo_ms2.msalign";
      MsAlignWriterPtr ms2_ptr = std::make_shared<MsAlignWriter>(file_name);
      MsHeaderPtr header_ptr = ms2_ptr_vec[0]->getMsHeaderPtr();
      DeconvPeakPtrVec peak_list;
      int sp_id = header_ptr->getId();
      for (size_t i = 0; i < features.size(); i++) {
        DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, i, features[i].getMonoMass(), features[i].getAbundance(),
                                                              features[i].getRepCharge(), features[i].getScore());
        peak_list.push_back(peak_ptr);
      }
      DeconvMsPtr ms_ptr = std::make_shared<DeconvMs>(header_ptr, peak_list);
      ms2_ptr->write(ms_ptr);

      double duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
      std::cout << "Feature Extraction duration (minutes): " << duration / 60.0 << '\n';
    }

  }  // namespace
}  // namespace toppic
