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

#ifndef TOPPIC_TOPDIA_COMMON_TOPDIA_PARA_HPP_
#define TOPPIC_TOPDIA_COMMON_TOPDIA_PARA_HPP_

#include <memory>
#include <string>

#include "topfd/common/topfd_para.hpp"

namespace toppic {

class TopdiaPara {
 public:
  TopdiaPara() {};

  std::string getParaStr(const std::string &prefix,
                         const std::string &sep,
                         TopfdParaPtr topfd_para);

  double getPseudoScoreCutoff() const { return pseudo_score_cutoff_; }
  int getPseudoMinPeaks() const { return pseudo_min_peaks_; }
  double getMs1SeedEnvInteCorrToleCutoff() const { return ms1_seed_env_inte_corr_tole_cutoff_; }
  double getMs2SeedEnvInteCorrToleCutoff() const { return ms2_seed_env_inte_corr_tole_cutoff_; }

  void setPseudoScoreCutoff(double pseudoScoreCutoff) { pseudo_score_cutoff_ = pseudoScoreCutoff; }
  void setPseudoMinPeaks(int pseudoMinPeaks) { pseudo_min_peaks_ = pseudoMinPeaks; }
  void setMs1SeedEnvInteCorrToleCutoff(double ms1SeedEnvInteCorrToleCutoff) { ms1_seed_env_inte_corr_tole_cutoff_ = ms1SeedEnvInteCorrToleCutoff; }
  void setMs2SeedEnvInteCorrToleCutoff(double ms2SeedEnvInteCorrToleCutoff) { ms2_seed_env_inte_corr_tole_cutoff_ = ms2SeedEnvInteCorrToleCutoff; }

private:
  double pseudo_score_cutoff_ = 0.55;
  int pseudo_min_peaks_ = 25;
  double ms1_seed_env_inte_corr_tole_cutoff_ = 0.5;
  double ms2_seed_env_inte_corr_tole_cutoff_ = 0;
};

typedef std::shared_ptr<TopdiaPara> TopdiaParaPtr;

}  // namespace toppic

#endif
