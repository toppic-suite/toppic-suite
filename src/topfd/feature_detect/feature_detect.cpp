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

namespace toppic {

namespace feature_detect {
std::vector<double> getSpectrumNoiseIntensities(const PeakPtrVec2D& raw_peaks){
  std::vector<double> spectrum_noise_levels;
  for (const auto& peaks : raw_peaks) {
    std::vector<double> intes;
    for (auto & peak : peaks) {
      intes.push_back(peak->getIntensity());
    }
    double noise = baseline_util::getBaseLine(intes);
    spectrum_noise_levels.push_back(noise);
  }
  return spectrum_noise_levels;
}

void process(int frac_id, const std::string &sp_file_name,
             bool missing_level_one, const std::string &resource_dir, const std::string &activation,
             bool isFaims, const std::vector<std::pair<double, int>> voltage_vec) {
  //logger::setLogLevel(2);
  FeatureParaPtr para_ptr
      = std::make_shared<FeaturePara>(frac_id, sp_file_name, resource_dir);
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
  std::string base_name = file_util::basename(sp_file_name);

  ///// params
  double bin_size = 0.1;
  double neighbor_mass_tole = 0.01;
  double mass_tole = 0.01;
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

      /// Prepare data
      std::vector<SeedEnvelope> seed_envs;
      for (auto & ms1_data : ms1_ptr_vec) {
        std::vector<DeconvPeakPtr> peaks =  ms1_data->getPeakPtrVec();
        for (const auto& peak : peaks)
          seed_envs.push_back(SeedEnvelope(peak));
      }
      std::sort(seed_envs.begin(), seed_envs.end(), SeedEnvelope::cmpInteDec);

//      PeakMatrix peak_matrix = PeakMatrix();
      PeakMatrix peak_matrix = PeakMatrix(raw_peaks, ms1_ptr_vec, bin_size, env_para_ptr->ms_one_sn_ratio_);
      peak_matrix.find_remove_non_neighbors(neighbor_mass_tole);

      /// Noise Intensity level
//      std::vector<double> spectrum_noise_levels;
      std::vector<double> spectrum_noise_levels = getSpectrumNoiseIntensities(raw_peaks);
      std::cout << "Number of spectra and seed envelopes: " << spectrum_noise_levels.size() << " , " << seed_envs.size() << std::endl;


      ///
      int seed_num = seed_envs.size();
      int env_coll_num = 0;
//      std::vector<EnvCollection> env_coll_list;
      std::vector<Feature> features;
      std::vector<bool> used;
      for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++){
        SeedEnvelope env = seed_envs[i];
        EnvCollection env_coll = env_coll_util::find_env_collection(peak_matrix, env, mass_tole, max_miss_env,
                                                                            max_miss_charge, max_miss_peak, para_max_charge,
                                                                            ratio_multi, match_peak_tole, env_para_ptr->ms_one_sn_ratio_);
        if (!env_coll.isEmpty()) {
          features.push_back(Feature(env_coll, peak_matrix, env_coll_num));
          env_coll.remove_matrix_peaks(peak_matrix);
//          env_coll_list.push_back(env_coll);
          env_coll_num = env_coll_num + 1;
//          break;
        }
      }
      std::cout << "Number of Envelope Collections: " << features.size() << std::endl;
      std::string file_name = base_name + "_" + file_num + "ms1.feature";
      write_feature::writeFeatures(file_name, features);
    }
  }
}

}  // namespace

}  // namespace toppic
