// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "topfd/deconv/deconv_single_sp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <mutex>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

#include "common/util/logger.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_filter.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/match_env_filter.hpp"
#include "ms/env/match_env_refine.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/spec/peak_list_util.hpp"
#include "topfd/deconv/deconv_data_util.hpp"
#include "topfd/dp/co_table.hpp"
#include "topfd/dp/dp_a.hpp"
#include "topfd/dp/dp_assign.hpp"
#include "topfd/dp/dp_para.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

// File-private helpers: the --output-dp-envs debugging dumps (win_envs.txt /
// dp_envs.txt) and the MS/MS harmonic check. They live in an unnamed
// namespace so they have internal
// linkage: they are private to this file, stay out of the class interface in
// the header, and their generic names (writeEnv, g_*_count, ...) cannot clash
// with same-named symbols elsewhere in the single toppic_common library.
namespace {

// Serializes the debug dumps and numbers the spectra across worker threads.
std::mutex g_dp_env_mutex;
int g_dp_env_spec_count = 0;
int g_cand_env_spec_count = 0;
int g_filtered_env_spec_count = 0;

// Write one envelope's information: a summary line followed by its peaks.
void writeEnv(std::ostream& os, const MatchEnvPtr& match_env) {
  EnvPtr theo_env = match_env->getTheoEnvPtr();
  int peak_num = theo_env->getPeakNum();
  os << "  env mono_mass " << theo_env->getMonoNeutralMass() << " charge "
     << theo_env->getCharge() << " intensity " << theo_env->compInteSum()
     << " peak_num " << peak_num << "\n";
  for (int k = 0; k < peak_num; k++) {
    os << "    " << theo_env->getMz(k) << " " << theo_env->getInte(k) << "\n";
  }
}

// Dump the non-empty windowed candidate envelopes (win_envs) and the
// DP-selected envelopes (dp_envs) of one spectrum to win_envs.txt /
// dp_envs.txt. The first spectrum of the run truncates the files; later spectra
// append.
void outputDpEnvs(const MatchEnvPtr2D& win_envs, const MatchEnvPtrVec& dp_envs,
                  int ms_level) {
  std::lock_guard<std::mutex> lock(g_dp_env_mutex);
  std::ios::openmode mode =
      (g_dp_env_spec_count == 0) ? std::ios::out : std::ios::app;
  int spec_index = g_dp_env_spec_count++;

  std::ofstream win_out("win_envs.txt", mode);
  win_out << "# spectrum " << spec_index << " ms_level " << ms_level << "\n";
  for (size_t w = 0; w < win_envs.size(); w++) {
    if (win_envs[w].empty()) {
      continue;
    }
    win_out << "window " << w << " " << win_envs[w].size() << " envelopes\n";
    for (size_t i = 0; i < win_envs[w].size(); i++) {
      writeEnv(win_out, win_envs[w][i]);
    }
  }

  std::ofstream dp_out("dp_envs.txt", mode);
  dp_out << "# spectrum " << spec_index << " ms_level " << ms_level << " "
         << dp_envs.size() << " envelopes\n";
  for (size_t i = 0; i < dp_envs.size(); i++) {
    writeEnv(dp_out, dp_envs[i]);
  }
}

// Dump the non-null candidate envelopes (per peak and charge) of one spectrum
// to cand_envs.txt, before any filtering. The first spectrum truncates the
// file; later spectra append.
void outputCandEnvs(const MatchEnvPtr2D& cand_envs, int ms_level) {
  std::lock_guard<std::mutex> lock(g_dp_env_mutex);
  std::ios::openmode mode =
      (g_cand_env_spec_count == 0) ? std::ios::out : std::ios::app;
  int spec_index = g_cand_env_spec_count++;

  std::ofstream out("cand_envs.txt", mode);
  out << "# spectrum " << spec_index << " ms_level " << ms_level << "\n";
  for (size_t i = 0; i < cand_envs.size(); i++) {
    for (size_t j = 0; j < cand_envs[i].size(); j++) {
      if (cand_envs[i][j] != nullptr) {
        out << "peak " << i << " charge " << (j + 1) << "\n";
        writeEnv(out, cand_envs[i][j]);
      }
    }
  }
}

// Dump the non-null candidate envelopes (per peak and charge) of one spectrum
// to filtered_envs.txt, after env_filter::filter has run. The first spectrum
// truncates the file; later spectra append.
void outputFilteredEnvs(const MatchEnvPtr2D& cand_envs, int ms_level) {
  std::lock_guard<std::mutex> lock(g_dp_env_mutex);
  std::ios::openmode mode =
      (g_filtered_env_spec_count == 0) ? std::ios::out : std::ios::app;
  int spec_index = g_filtered_env_spec_count++;

  std::ofstream out("filtered_envs.txt", mode);
  out << "# spectrum " << spec_index << " ms_level " << ms_level << "\n";
  for (size_t i = 0; i < cand_envs.size(); i++) {
    for (size_t j = 0; j < cand_envs[i].size(); j++) {
      if (cand_envs[i][j] != nullptr) {
        out << "peak " << i << " charge " << (j + 1) << "\n";
        writeEnv(out, cand_envs[i][j]);
      }
    }
  }
}

// The isotope cosine of a matched envelope: the cosine between the
// theoretical peak intensities and the observed ones over the envelope's
// positions, a missing peak contributing zero to the observed vector. This is
// what FLASHDeconv tests a harmonic parent with (its isotope cosine against
// the averagine distribution, kMinIsotopeCosine = 0.85); TopFD's theoretical
// envelope is the averagine distribution scaled to the observed intensities,
// so the same quantity is the cosine between the two vectors here.
double compIsotopeCosine(const MatchEnvPtr& env) {
  EnvPtr theo_env = env->getTheoEnvPtr();
  ExpEnvPtr exp_env = env->getExpEnvPtr();
  int peak_num = theo_env->getPeakNum();
  if (peak_num != exp_env->getPeakNum() || peak_num == 0) {
    return 0;
  }
  double dot = 0;
  double theo_norm = 0;
  double exp_norm = 0;
  for (int i = 0; i < peak_num; i++) {
    double theo_inte = theo_env->getInte(i);
    double exp_inte = exp_env->isExist(i) ? exp_env->getInte(i) : 0;
    dot += theo_inte * exp_inte;
    theo_norm += theo_inte * theo_inte;
    exp_norm += exp_inte * exp_inte;
  }
  if (theo_norm <= 0 || exp_norm <= 0) {
    return 0;
  }
  return dot / std::sqrt(theo_norm * exp_norm);
}

// The signal-to-noise ratio of a matched envelope, in the sense
// PeakGroup::updateSnr uses: of the observed power, the part the theoretical
// shape explains (cos^2) against the part it does not ((1 - cos^2)) plus the
// noise floor of the spectrum over the envelope's positions.
double compEnvSnr(const MatchEnvPtr& env, double noise_inte) {
  double cosine = compIsotopeCosine(env);
  ExpEnvPtr exp_env = env->getExpEnvPtr();
  int peak_num = exp_env->getPeakNum();
  double signal_power = 0;
  for (int i = 0; i < peak_num; i++) {
    if (exp_env->isExist(i)) {
      double inte = exp_env->getInte(i);
      signal_power += inte * inte;
    }
  }
  double cos_squared = cosine * cosine;
  double noise_power = peak_num * noise_inte * noise_inte;
  double denominator = (1 - cos_squared) * signal_power + noise_power;
  if (denominator <= 0) {
    return 0;
  }
  return cos_squared * signal_power / denominator;
}

// The harmonic check, ported from the FLASHDeconv implementation in this
// repository (FlashDeconvSingleSp::removeHarmonicPeakGroups and the harmonic
// vetoes of updateCandidateMassBins), applied to the envelopes the dynamic
// programming selected.
//
// An envelope of mass M at charge z whose peaks are every k-th isotopic
// position of a real envelope of mass k*M at charge k*z is a harmonic
// artefact: the peaks in between belong to the parent, and TopFD's own filters
// cannot see them, because they compare envelopes only at a shared reference
// peak. For each selected envelope the parents of k = 2 and 3 are built from
// the same peak list and scored; a parent that matches more peaks with a
// better MS-Deconv score replaces the harmonic. The parent is not required to
// pass envelope validation - that validation is what hides these masses - but
// it must match more peaks than the validation of its own mass group requires.
//
// The parent is built by detectEnvByRefPeak from the peak nearest its
// theoretical reference position, the way every other selected envelope is
// built: match_env_refine::mzRefine, which runs next, re-derives the reference
// distribution from the experimental reference peak and reads as many peaks
// from it as the experimental envelope has on each side, which an envelope
// built around its monoisotopic peak does not guarantee.
void replaceHarmonicEnvs(MatchEnvPtrVec& envs, const PeakPtrVec& peak_list,
                         const DeconvDataPtr& data_ptr,
                         const EnvParaPtr& env_para_ptr) {
  const std::vector<int> harmonic_factors = {2, 3};
  // FLASHDeconv's kMinIsotopeCosine, the floor its own harmonic removal uses
  const double min_isotope_cosine = 0.85;
  double min_inte = data_ptr->getMinInte();
  double min_ref_inte = data_ptr->getMinRefInte();
  double bound_inte = data_ptr->getNoiseInte();
  double max_mass = data_ptr->getMaxMass();
  int max_charge = data_ptr->getMaxCharge();
  bool is_full = true;
  bool check_real_env = false;
  // the masses already reported, so that a parent is not added twice
  std::set<std::pair<long, int> > reported;
  auto mass_key = [](double mass, int charge) {
    return std::make_pair(std::lround(mass * 10), charge);
  };
  for (const MatchEnvPtr& env : envs) {
    reported.insert(mass_key(env->getTheoEnvPtr()->getMonoNeutralMass(),
                             env->getTheoEnvPtr()->getCharge()));
  }
  for (std::size_t i = 0; i < envs.size(); i++) {
    double mass = envs[i]->getTheoEnvPtr()->getMonoNeutralMass();
    int charge = envs[i]->getTheoEnvPtr()->getCharge();
    int match_num = envs[i]->getExpEnvPtr()->getMatchPeakNum();
    double score = envs[i]->getMsdeconvScore();
    MatchEnvPtr best_parent = nullptr;
    for (int factor : harmonic_factors) {
      double parent_mass = mass * factor;
      int parent_charge = charge * factor;
      if (parent_charge > max_charge || parent_mass >= max_mass) {
        continue;
      }
      // a cheap pre-check on the nominal mass; the mass actually built is
      // checked against the same set below, because the two can fall in
      // different buckets
      if (reported.count(mass_key(parent_mass, parent_charge)) > 0) {
        continue;
      }
      // where the parent's reference (most intense) peak would be
      MatchEnvPtr probe = env_detect::detectEnvByMonoMass(
          peak_list, parent_mass, parent_charge, min_inte, bound_inte,
          env_para_ptr, is_full);
      if (probe == nullptr || probe->getTheoEnvPtr()->getPeakNum() == 0) {
        continue;
      }
      double tolerance = env_para_ptr->getMzTolerance(parent_charge);
      int ref_peak = peak_list_util::getNearPeakIdx(
          peak_list, probe->getTheoEnvPtr()->getReferMz(), tolerance);
      if (ref_peak < 0) {
        continue;
      }
      MatchEnvPtr parent = env_detect::detectEnvByRefPeak(
          peak_list, ref_peak, parent_charge, max_mass, min_inte, min_ref_inte,
          bound_inte, env_para_ptr, is_full, check_real_env);
      if (parent == nullptr || parent->getExpEnvPtr()->getPeakNum() == 0) {
        continue;
      }
      // the envelope built from the reference peak must be the parent mass
      double built_mass = parent->getTheoEnvPtr()->getMonoNeutralMass();
      if (std::abs(built_mass - parent_mass) >
          std::max(parent_mass * 2e-5, 0.05)) {
        continue;
      }
      if (reported.count(mass_key(built_mass, parent_charge)) > 0) {
        continue;
      }
      int mass_group = env_para_ptr->getMassGroup(built_mass);
      if (mass_group < 0) {
        continue;
      }
      if (parent->getExpEnvPtr()->getMatchPeakNum() <=
          std::max(match_num, env_para_ptr->min_match_peak_num_[mass_group])) {
        continue;
      }
      parent->compMsdeconvScr(env_para_ptr);
      if (parent->getMsdeconvScore() <= score) {
        continue;
      }
      // The shape guard of the reference implementation
      // (FlashDeconvSingleSp::removeHarmonicPeakGroups requires the parent to
      // score under kMinIsotopeCosine and to have the higher SNR). Without it
      // the two tests above are nearly one test: both the matched peak count
      // and the MS-Deconv score grow with the number of matched peaks, and a
      // high charge envelope in a dense region matches many peaks by chance.
      if (compIsotopeCosine(parent) < min_isotope_cosine) {
        continue;
      }
      if (compEnvSnr(parent, bound_inte) <= compEnvSnr(envs[i], bound_inte)) {
        continue;
      }
      if (best_parent == nullptr ||
          parent->getMsdeconvScore() > best_parent->getMsdeconvScore()) {
        best_parent = parent;
      }
    }
    if (best_parent != nullptr) {
      reported.erase(mass_key(mass, charge));
      reported.insert(
          mass_key(best_parent->getTheoEnvPtr()->getMonoNeutralMass(),
                   best_parent->getTheoEnvPtr()->getCharge()));
      envs[i] = best_parent;
    }
  }
}

}  // namespace

DeconvSingleSp::DeconvSingleSp(const TopfdParaPtr& topfd_para_ptr,
                               PeakPtrVec& peak_list, int ms_level,
                               double max_mass, int max_charge) {
  topfd_para_ptr_ = topfd_para_ptr;
  env_para_ptr_ = std::make_shared<EnvPara>(
      topfd_para_ptr->getMzError(), topfd_para_ptr->getMaxMissPeakNum(),
      topfd_para_ptr->isRunFilterByMz(),
      topfd_para_ptr->getMassGroupBoundary());
  // EnvPara is shared by the MS1 and the MS/MS deconvolution, so the window of
  // the filtering-by-mz step is set for this level here rather than through the
  // EnvPara constructor
  env_para_ptr_->rank_peak_distance_ =
      (ms_level == 2) ? topfd_para_ptr->getMsTwoRankPeakDistance()
                      : topfd_para_ptr->getMsOneRankPeakDistance();
  dp_para_ptr_ = std::make_shared<DpPara>(topfd_para_ptr->getMzError());
  dp_para_ptr_->max_peak_sharing_ratio_ =
      topfd_para_ptr->getMaxPeakSharingRatio();
  ms_level_ = ms_level;
  double sn_ratio;
  if (ms_level == 1) {
    sn_ratio = topfd_para_ptr_->getMsOneSnRatio();
  } else {
    sn_ratio = topfd_para_ptr_->getMsTwoSnRatio();
  }
  data_ptr_ = deconv_data_util::getDataPtr(
      peak_list, max_mass, max_charge, dp_para_ptr_->dp_window_size_,
      topfd_para_ptr_->isEstimateMinInte(), sn_ratio);
}

void DeconvSingleSp::postprocess(MatchEnvPtrVec& dp_envs, int ms_level) {
  // assign intensity
  PeakPtrVec peak_list = data_ptr_->getPeakList();

  // replace harmonic envelopes by their parent masses
  if (ms_level == 2 && topfd_para_ptr_->isMsTwoHarmonicCheck()) {
    replaceHarmonicEnvs(dp_envs, peak_list, data_ptr_, env_para_ptr_);
  }
  match_env_util::assignIntensity(peak_list, dp_envs);

  // refinement
  match_env_refine::mzRefine(dp_envs, env_para_ptr_->refine_core_ratio_);
  result_envs_ = dp_envs;

  if (topfd_para_ptr_->isKeepUnusedPeaks()) {
    // all added envelopes have charge 1, so we use charge 1 to get mz_tolerance
    double mz_tole = env_para_ptr_->getMzTolerance(1);
    result_envs_ =
        match_env_util::addUnusedMasses(result_envs_, peak_list, mz_tole);
  }

  // Obtain EnvCNN Scores for envelopes
  onnx_env_cnn::computeEnvScores(peak_list, result_envs_);
  // filtering by EnvCNN score
  if (ms_level == 2) {
    result_envs_ = match_env_filter::filterByEnvCnnScore(
        result_envs_, topfd_para_ptr_->getMs2EnvCnnScoreCutoff());
  }
  // filtering by estimated number of amino acids
  if (topfd_para_ptr_->isAANumBasedFilter()) {
    result_envs_ = match_env_filter::filterByAANum(
        result_envs_, data_ptr_->getMaxMass(),
        topfd_para_ptr_->isSortUseMsDeconv(), env_para_ptr_);
  }
  // sorting
  if (topfd_para_ptr_->isSortUseMsDeconv()) {
    std::sort(result_envs_.begin(), result_envs_.end(),
              MatchEnv::cmpMsdeconvScoreDec);
  } else {
    std::sort(result_envs_.begin(), result_envs_.end(),
              MatchEnv::cmpEnvcnnScoreDec);
  }
  // reassign intensity
  match_env_util::assignIntensity(peak_list, result_envs_);
}

MatchEnvPtrVec DeconvSingleSp::deconv() {
  if (data_ptr_ == nullptr) {
    return result_envs_;
  }
  PeakPtrVec peak_list = data_ptr_->getPeakList();
  // envelope detection; the candidates that env_filter::filter's first
  // test would remove are not built (check_real_env)
  bool check_real_env = true;
  MatchEnvPtr2D cand_envs = env_detect::getCandidateEnv(
      peak_list, data_ptr_->getMaxCharge(), data_ptr_->getMaxMass(),
      data_ptr_->getMinInte(), data_ptr_->getMinRefInte(),
      data_ptr_->getNoiseInte(), env_para_ptr_, check_real_env);
  LOG_DEBUG("candidate complete");

  if (topfd_para_ptr_->isOutputDpEnvs()) {
    outputCandEnvs(cand_envs, ms_level_);
  }

  // envelope filter
  env_filter::filter(cand_envs, peak_list, env_para_ptr_);

  if (topfd_para_ptr_->isOutputDpEnvs()) {
    outputFilteredEnvs(cand_envs, ms_level_);
  }

  // prepare for dp
  MatchEnvPtr2D win_envs = dp_assign::assignWinEnv(
      cand_envs, data_ptr_->getWinNum(), data_ptr_->getWinIdVec(),
      dp_para_ptr_->env_num_per_win_);

  // prepare table for dp
  if (dp_para_ptr_->check_double_increase_) {
    LOG_DEBUG("Generating coexisting table...");
    dp_para_ptr_->coexist_table_ = CoTable::initCoexistTable(
        win_envs, dp_para_ptr_->mz_tolerance_,
        dp_para_ptr_->max_peak_sharing_ratio_);
  }
  // dp
  LOG_DEBUG("Generating Graph and DP...");
  DpA dp(data_ptr_, win_envs, dp_para_ptr_, dp_para_ptr_->mz_tolerance_);
  MatchEnvPtrVec dp_envs = dp.getResult();

  if (topfd_para_ptr_->isOutputDpEnvs()) {
    outputDpEnvs(win_envs, dp_envs, ms_level_);
  }

  postprocess(dp_envs, ms_level_);

  return result_envs_;
}

}  // namespace toppic
