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

#include "topfd/deconv/greedy_ms1_deconv.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/exp_env.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

namespace greedy_ms1_deconv {

namespace {

// The parameters of the search. The defaults are the setting that gave the
// most identifications on 480_F3_02 in toppic_test experiment 104 (runs g3
// and g6, with the filters of g5 and g6 that report fewer envelopes at the
// same identifications). They can be changed for research with the
// environment variable TOPFD_GREEDY_PARA, a comma-separated list of
// name=value pairs (e.g. "min_cos=0.7,keep_ratio=2").
struct GreedyPara {
  // isotope cosine an envelope of the seed charge needs
  double min_cos = 0.6;
  // isotope cosine an envelope of a neighbouring charge needs
  double nb_min_cos = 0.6;
  // charges in a row that may fail inside a charge series
  int nb_max_gap = 1;
  // isotopes with theoretical intensity >= core_ratio x the apex are used to
  // fit the scale of an envelope
  double core_ratio = 0.5;
  // isotopes below this fraction of the apex are never used
  double min_rel_inte = 0.01;
  // a matched peak with residual intensity > keep_ratio x its theoretical
  // intensity keeps the difference after the subtraction
  double keep_ratio = 1.5;
  // weight of the theoretical intensity predicted where there is no peak
  double miss_weight = 1.0;
  // the score of an envelope: explained intensity minus the weighted
  // intensity predicted where there is no peak (0), or that times the
  // squared isotope cosine (1; the cos^2 weighting of FLASHDeconv's SNR), so
  // that a charge 2z envelope that matches the peaks of a charge z envelope
  // and unrelated peaks between them in a dense region loses to charge z
  int score_mode = 1;
  // monoisotopic offsets (isotopes) searched around the apex alignment
  int max_offset = 2;
  // times a peak may be the seed
  int max_seed_num = 2;
  // whether a peak matched in an envelope found before (and kept, with the
  // theoretical intensity deducted) may still be a seed (1) or not (0)
  int seed_matched = 1;
  // the seed envelope must explain at least this fraction of the residual
  // intensity in its m/z window (0 = no test)
  double min_explained_frac = 0;
  // the monoisotopic mass: the offset whose series scores best (0), or the
  // offset whose averagine distribution has the best cosine with the
  // isotope intensities summed over the series' charges, with the observed
  // intensity beyond the distribution counted (1; FLASHDeconv's
  // getIsotopeCosineAndIsoOffset)
  int mono_by = 1;
  // the envelope boundary of the subtraction found from the observed peaks
  // (1) or the isotopes fitted (0): from the apex outwards, the envelope goes
  // on while the next isotope has a peak whose residual intensity is at most
  // dyn_rise x that of the isotope before it; inside the boundary the tail
  // isotopes (below core_ratio) and the isotopes beyond the fitted ones are
  // removed whole, so that a tail wider than the averagine distribution
  // leaves no residual to seed an isotope-shifted copy of the mass
  int dyn_bound = 1;
  double dyn_rise = 1.2;
  // isotopes beyond the fitted ones the boundary may reach
  int dyn_max_ext = 6;
  // a charge z is tried for a seed only when a peak lies one isotope
  // spacing (getIsotopeMass() / z) from it on either side (1; DeconTools'
  // charge check) or always (0)
  int precheck = 1;
  // harmonic test of the seed charge z (0 = none): for f = 2, 3 dividing z,
  // the observed / theoretical intensity ratio of the isotopes off the grid
  // of every f-th isotope through the apex, over that of the isotopes on it,
  // must be at least harm_ratio; an envelope of charge z that is really
  // one of charge z / f matches the peaks on the grid and only unrelated
  // peaks between them (FLASHDeconv's and kDecon's harmonic checks)
  double harm_ratio = 0.4;
  // a mass found at one charge only is reported (and subtracted) only when
  // its envelope's isotope cosine reaches single_min_cos (0 = no test)
  double single_min_cos = 0.8;
  // a seed whose mass lies within merge_window Da of a mass found before
  // gives that mass's envelopes at the seed series' charges the earlier
  // mass lacks, when they fit, and nothing otherwise (0 = no merge): the
  // peaks of a charge state missing from a mass's series must not become a
  // near-duplicate mass
  double merge_window = 0;
  // envelopes with an EnvCNN score below min_envcnn are not reported (0 =
  // all reported): they are a third of the envelopes, seed the features
  // that the feature detection would otherwise spend 50% more time on, and
  // gave no identification (experiment 104, run g6)
  double min_envcnn = 0.2;
  // research aid for timing: 0 skips the EnvCNN scores
  int envcnn = 1;
  // main's precursor-window step first (1) or not (0)
  int prec_win = 0;
  // an envelope within dup_window Da of a mass already found at the same
  // charge is not used (0 = no test): the residual of a subtracted envelope
  // must not come back as an isotope-shifted copy of its mass
  double dup_window = 3.5;
  // the postprocessing step (1) or not (0): every precursor window of the
  // following MS/MS spectra in which no reported envelope has its apex gets
  // the envelope of its most intense peak (addPrecWinEnvs; toppic_test
  // experiment 113)
  int prec_win_post = 1;
};

GreedyPara readPara() {
  GreedyPara para;
  const char* str = std::getenv("TOPFD_GREEDY_PARA");
  if (str == nullptr) {
    return para;
  }
  std::stringstream ss(str);
  std::string item;
  while (std::getline(ss, item, ',')) {
    std::size_t pos = item.find('=');
    if (pos == std::string::npos) {
      continue;
    }
    std::string name = item.substr(0, pos);
    double value = std::atof(item.substr(pos + 1).c_str());
    if (name == "min_cos") {
      para.min_cos = value;
    } else if (name == "nb_min_cos") {
      para.nb_min_cos = value;
    } else if (name == "nb_max_gap") {
      para.nb_max_gap = static_cast<int>(value);
    } else if (name == "core_ratio") {
      para.core_ratio = value;
    } else if (name == "min_rel_inte") {
      para.min_rel_inte = value;
    } else if (name == "keep_ratio") {
      para.keep_ratio = value;
    } else if (name == "miss_weight") {
      para.miss_weight = value;
    } else if (name == "score_mode") {
      para.score_mode = static_cast<int>(value);
    } else if (name == "max_offset") {
      para.max_offset = static_cast<int>(value);
    } else if (name == "max_seed_num") {
      para.max_seed_num = static_cast<int>(value);
    } else if (name == "prec_win") {
      para.prec_win = static_cast<int>(value);
    } else if (name == "dup_window") {
      para.dup_window = value;
    } else if (name == "seed_matched") {
      para.seed_matched = static_cast<int>(value);
    } else if (name == "min_explained_frac") {
      para.min_explained_frac = value;
    } else if (name == "mono_by") {
      para.mono_by = static_cast<int>(value);
    } else if (name == "single_min_cos") {
      para.single_min_cos = value;
    } else if (name == "merge_window") {
      para.merge_window = value;
    } else if (name == "min_envcnn") {
      para.min_envcnn = value;
    } else if (name == "envcnn") {
      para.envcnn = static_cast<int>(value);
    } else if (name == "harm_ratio") {
      para.harm_ratio = value;
    } else if (name == "precheck") {
      para.precheck = static_cast<int>(value);
    } else if (name == "dyn_bound") {
      para.dyn_bound = static_cast<int>(value);
    } else if (name == "dyn_rise") {
      para.dyn_rise = value;
    } else if (name == "dyn_max_ext") {
      para.dyn_max_ext = static_cast<int>(value);
    } else if (name == "prec_win_post") {
      para.prec_win_post = static_cast<int>(value);
    } else {
      LOG_WARN("Unknown TOPFD_GREEDY_PARA entry " << name);
    }
  }
  return para;
}

const GreedyPara& getPara() {
  static const GreedyPara para = readPara();
  return para;
}

// The averagine distribution of a mass bin of EnvBase: per isotope, the mass
// offset from the monoisotopic peak and the intensity relative to the apex.
struct Distr {
  std::vector<double> offsets;
  std::vector<double> intes;
  int apex = 0;
  // the isotopes above GreedyPara::min_rel_inte
  int first = 0;
  int last = -1;
};

// The theoretical envelope of one mass at one charge, fitted to the residual
// spectrum.
struct Fit {
  bool valid = false;
  double mono_mass = 0;
  int charge = 0;
  // theoretical intensity of the apex isotope (scale x relative intensity)
  double scale = 0;
  double cos = 0;
  int match_num = 0;
  int miss_num = 0;
  double explained = 0;
  double penalty = 0;
  bool apex_match = false;
  // the residual intensity of all peaks in the m/z window of the isotopes
  // used (half an isotope spacing beyond both ends)
  double window_inte = 0;
  // the isotopes used: [first, last] of the distribution
  int first = 0;
  int last = -1;
  // per isotope of the distribution: the matched peak or -1
  std::vector<int> peak_idxes;
  double score() const {
    double score = explained - getPara().miss_weight * penalty;
    return getPara().score_mode == 1 ? score * cos * cos : score;
  }
};

using FitVec = std::vector<Fit>;

// The envelope detection redone in a precursor window (fitWindow).
struct WinFit {
  // the envelope of the window's charge; not valid when there is none
  Fit fit;
  // the peaks in the window and the most intense one (-1: none)
  int peak_num = 0;
  int seed = -1;
  // the score of the mass's series and its number of charges
  double score = 0;
  int charge_num = 0;
  // the monoisotopic offset applied to the apex alignment (isotopes)
  int offset = 0;
};

// the envelope of a fit as it is reported (defined below)
MatchEnvPtr getMatchEnv(const PeakPtrVec& peak_list, const Fit& fit,
                        const EnvParaPtr& env_para_ptr, double min_inte);

// One spectrum's search.
class GreedySearch {
 public:
  // sn_ratio: the signal-to-noise ratio of the seed threshold; max_mass,
  // max_charge: the bounds of the masses and charges; bound (may be empty):
  // the charges a mass may have (deconvSeries); max_shared_peak_num: the
  // limit of the peaks an envelope may share with each envelope reported
  // before (-1: no limit; see removeSharedEnvs)
  GreedySearch(const PeakPtrVec& peak_list, const TopfdParaPtr& topfd_para_ptr,
               const EnvParaPtr& env_para_ptr, double sn_ratio,
               double max_mass, int max_charge, ChargeBound bound,
               int max_shared_peak_num, std::ostream* debug_os);

  // the charge series of the masses found, in the order they were found
  std::vector<FitVec> run();

  // the envelope of the most intense peak in an m/z window, on a search
  // that has not subtracted anything
  WinFit fitWindow(double bgn, double end);

 private:
  const GreedyPara& para_;
  EnvParaPtr env_para_ptr_;
  ChargeBound bound_;
  std::ostream* debug_os_;
  int max_charge_;
  double max_mass_;
  double min_inte_ = 0;
  double min_ref_inte_ = 0;
  std::vector<double> mzs_;
  std::vector<double> intes_;
  std::vector<double> residuals_;
  std::unordered_map<int, Distr> distr_cache_;
  // the masses found so far, per charge
  std::vector<std::vector<double>> found_masses_;
  PeakPtrVec peak_list_;
  int max_shared_peak_num_;
  // per peak of the spectrum, the reported envelopes that match it (indexes
  // into reported_envs_); the mass and charge of each reported envelope
  std::vector<std::vector<int>> peak_envs_;
  std::vector<std::pair<double, int>> reported_envs_;

  const Distr& getDistr(double mono_mass);

  int findPeak(double mz, double tolerance) const;

  Fit fitEnv(double mono_mass, int charge);

  bool passes(const Fit& fit, double min_cos) const;

  FitVec getSeries(double mono_mass, int charge, const Fit& seed_fit);

  static double seriesScore(const FitVec& series);

  int getAggregateOffset(const FitVec& series);

  std::vector<int> getTailPeaks(const Fit& fit);

  bool isHarmonic(const Fit& fit);

  void subtract(const FitVec& series);

  FitVec removeSharedEnvs(FitVec series);
};

GreedySearch::GreedySearch(const PeakPtrVec& peak_list,
                           const TopfdParaPtr& topfd_para_ptr,
                           const EnvParaPtr& env_para_ptr, double sn_ratio,
                           double max_mass, int max_charge, ChargeBound bound,
                           int max_shared_peak_num, std::ostream* debug_os)
    : para_(getPara()),
      env_para_ptr_(env_para_ptr),
      bound_(std::move(bound)),
      debug_os_(debug_os),
      max_charge_(max_charge),
      max_mass_(max_mass),
      peak_list_(peak_list),
      max_shared_peak_num_(max_shared_peak_num) {
  mzs_.reserve(peak_list.size());
  intes_.reserve(peak_list.size());
  for (const PeakPtr& peak : peak_list) {
    mzs_.push_back(peak->getPosition());
    intes_.push_back(peak->getIntensity());
  }
  residuals_ = intes_;
  found_masses_.resize(max_charge_ + 1);
  if (max_shared_peak_num_ >= 0) {
    peak_envs_.resize(peak_list.size());
  }
  // the noise level and the thresholds of DeconvData::initMinInte
  double noise_inte = baseline_util::getBaseLine(intes_);
  min_ref_inte_ = noise_inte * sn_ratio;
  min_inte_ = sn_ratio < 1.0 ? noise_inte * sn_ratio : noise_inte;
}

const Distr& GreedySearch::getDistr(double mono_mass) {
  int key = static_cast<int>(mono_mass / 10);
  auto it = distr_cache_.find(key);
  if (it != distr_cache_.end()) {
    return it->second;
  }
  Distr distr;
  EnvPtr env = EnvBase::getEnvByMonoMass(mono_mass);
  if (env != nullptr && env->getPeakNum() > 0) {
    double mono = env->getMonoMz();
    double max_inte = env->getReferInte();
    distr.apex = env->getReferIdx();
    for (int i = 0; i < env->getPeakNum(); i++) {
      distr.offsets.push_back(env->getMz(i) - mono);
      distr.intes.push_back(env->getInte(i) / max_inte);
    }
    distr.first = distr.apex;
    distr.last = distr.apex;
    while (distr.first > 0 &&
           distr.intes[distr.first - 1] >= para_.min_rel_inte) {
      distr.first--;
    }
    while (distr.last + 1 < static_cast<int>(distr.intes.size()) &&
           distr.intes[distr.last + 1] >= para_.min_rel_inte) {
      distr.last++;
    }
  }
  return distr_cache_.emplace(key, std::move(distr)).first->second;
}

// the peak nearest to mz within the tolerance, or -1
int GreedySearch::findPeak(double mz, double tolerance) const {
  auto it = std::lower_bound(mzs_.begin(), mzs_.end(), mz);
  int idx = -1;
  double best = tolerance;
  if (it != mzs_.end()) {
    int i = static_cast<int>(it - mzs_.begin());
    if (mzs_[i] - mz <= best) {
      best = mzs_[i] - mz;
      idx = i;
    }
  }
  if (it != mzs_.begin()) {
    int i = static_cast<int>(it - mzs_.begin()) - 1;
    if (mz - mzs_[i] <= best) {
      idx = i;
    }
  }
  return idx;
}

Fit GreedySearch::fitEnv(double mono_mass, int charge) {
  Fit fit;
  fit.mono_mass = mono_mass;
  fit.charge = charge;
  if (mono_mass < env_para_ptr_->min_mass_ || mono_mass >= max_mass_) {
    return fit;
  }
  // the charges the bound excludes (seeds, series and offsets alike)
  if (bound_ && !bound_(mono_mass, charge)) {
    return fit;
  }
  const Distr& distr = getDistr(mono_mass);
  if (distr.last < distr.first) {
    return fit;
  }
  double mono_mz = mono_mass / charge + mass_constant::getProtonMass();
  double tolerance = env_para_ptr_->getMzTolerance(charge);
  int n = static_cast<int>(distr.intes.size());
  fit.peak_idxes.assign(n, -1);
  std::vector<double> obs(n, 0.0);
  for (int i = distr.first; i <= distr.last; i++) {
    int idx = findPeak(mono_mz + distr.offsets[i] / charge, tolerance);
    if (idx >= 0 && intes_[idx] >= min_inte_) {
      // one peak serves one isotope: keep the nearer one
      if (i > distr.first && fit.peak_idxes[i - 1] == idx) {
        double prev_err = std::abs(mzs_[idx] - mono_mz -
                                   distr.offsets[i - 1] / charge);
        double cur_err =
            std::abs(mzs_[idx] - mono_mz - distr.offsets[i] / charge);
        if (cur_err >= prev_err) {
          continue;
        }
        fit.peak_idxes[i - 1] = -1;
        obs[i - 1] = 0;
      }
      fit.peak_idxes[i] = idx;
      obs[i] = residuals_[idx];
    }
  }
  // the scale from the core isotopes (least squares)
  double tt = 0;
  double to = 0;
  for (int i = distr.first; i <= distr.last; i++) {
    if (distr.intes[i] >= para_.core_ratio && obs[i] > 0) {
      tt += distr.intes[i] * distr.intes[i];
      to += distr.intes[i] * obs[i];
    }
  }
  if (to <= 0) {
    return fit;
  }
  fit.scale = to / tt;
  // the isotopes whose theoretical intensity is above the noise level
  fit.first = distr.apex;
  fit.last = distr.apex;
  while (fit.first > distr.first &&
         fit.scale * distr.intes[fit.first - 1] >= min_inte_) {
    fit.first--;
  }
  while (fit.last < distr.last &&
         fit.scale * distr.intes[fit.last + 1] >= min_inte_) {
    fit.last++;
  }
  double oo = 0;
  tt = 0;
  to = 0;
  for (int i = fit.first; i <= fit.last; i++) {
    double theo = fit.scale * distr.intes[i];
    tt += distr.intes[i] * distr.intes[i];
    to += distr.intes[i] * obs[i];
    oo += obs[i] * obs[i];
    if (fit.peak_idxes[i] < 0) {
      fit.miss_num++;
      fit.penalty += theo;
    } else if (obs[i] > 0) {
      fit.match_num++;
      fit.explained += std::min(obs[i], theo);
    }
  }
  fit.apex_match = obs[distr.apex] > 0;
  if (para_.min_explained_frac > 0) {
    double half = 0.5 * mass_constant::getIsotopeMass() / charge;
    double bgn = mono_mz + distr.offsets[fit.first] / charge - half;
    double end = mono_mz + distr.offsets[fit.last] / charge + half;
    for (auto it = std::lower_bound(mzs_.begin(), mzs_.end(), bgn);
         it != mzs_.end() && *it <= end; ++it) {
      fit.window_inte += residuals_[it - mzs_.begin()];
    }
  }
  fit.cos = (oo > 0 && tt > 0) ? to / std::sqrt(tt * oo) : 0;
  fit.valid = true;
  return fit;
}

bool GreedySearch::passes(const Fit& fit, double min_cos) const {
  if (!fit.valid || !fit.apex_match || fit.cos < min_cos) {
    return false;
  }
  int mass_group = env_para_ptr_->getMassGroup(fit.mono_mass);
  if (mass_group < 0 ||
      fit.match_num < env_para_ptr_->min_match_peak_num_[mass_group]) {
    return false;
  }
  if (para_.dup_window > 0) {
    for (double mass : found_masses_[fit.charge]) {
      if (std::abs(mass - fit.mono_mass) < para_.dup_window) {
        return false;
      }
    }
  }
  return true;
}

// The envelopes of one mass: the seed charge and the neighbouring charges
// that fit, going outwards until nb_max_gap + 1 charges in a row fail.
FitVec GreedySearch::getSeries(double mono_mass, int charge,
                               const Fit& seed_fit) {
  FitVec series = {seed_fit};
  for (int dir = -1; dir <= 1; dir += 2) {
    int gap = 0;
    for (int z = charge + dir; z >= 1 && z <= max_charge_; z += dir) {
      Fit fit = fitEnv(mono_mass, z);
      if (passes(fit, para_.nb_min_cos)) {
        series.push_back(std::move(fit));
        gap = 0;
      } else {
        gap++;
        if (gap > para_.nb_max_gap) {
          break;
        }
      }
    }
  }
  return series;
}

double GreedySearch::seriesScore(const FitVec& series) {
  double score = 0;
  for (const Fit& fit : series) {
    score += fit.score();
  }
  return score;
}

// The isotope offset k in [-max_offset, max_offset] of the series' mass
// whose averagine distribution, shifted by k isotopes, has the best cosine
// with the residual intensities at the isotope positions summed over the
// series' charges. The positions run max_offset isotopes beyond both ends of
// the distribution, so intensity the shifted distribution does not explain
// lowers its cosine.
int GreedySearch::getAggregateOffset(const FitVec& series) {
  double base_mass = series[0].mono_mass;
  const Distr& distr = getDistr(base_mass);
  int k_max = para_.max_offset;
  int bgn = series[0].first - k_max;
  int end = series[0].last + k_max;
  for (const Fit& fit : series) {
    bgn = std::min(bgn, fit.first - k_max);
    end = std::max(end, fit.last + k_max);
  }
  double iso_mass = mass_constant::getIsotopeMass();
  int n = static_cast<int>(distr.intes.size());
  std::vector<double> aggr(end - bgn + 1, 0.0);
  for (const Fit& fit : series) {
    int z = fit.charge;
    double mono_mz = base_mass / z + mass_constant::getProtonMass();
    double tolerance = env_para_ptr_->getMzTolerance(z);
    for (int j = bgn; j <= end; j++) {
      double offset = (j >= 0 && j < n) ? distr.offsets[j] : j * iso_mass;
      int idx = findPeak(mono_mz + offset / z, tolerance);
      if (idx >= 0 && intes_[idx] >= min_inte_) {
        aggr[j - bgn] += residuals_[idx];
      }
    }
  }
  double aa = 0;
  for (double v : aggr) {
    aa += v * v;
  }
  int best_k = 0;
  double best_cos = -1;
  for (int k = -k_max; k <= k_max; k++) {
    double ab = 0;
    double bb = 0;
    for (int j = bgn; j <= end; j++) {
      int i = j - k;
      if (i < distr.first || i > distr.last) {
        continue;
      }
      ab += aggr[j - bgn] * distr.intes[i];
      bb += distr.intes[i] * distr.intes[i];
    }
    double cos = (aa > 0 && bb > 0) ? ab / std::sqrt(aa * bb) : 0;
    // ties keep the apex alignment
    if (cos > best_cos + 1e-9 || (k == 0 && cos >= best_cos - 1e-9)) {
      best_cos = cos;
      best_k = k;
    }
  }
  return best_k;
}

// Whether the envelope of a fit looks like an envelope of charge z / f,
// f = 2 or 3 (see GreedyPara::harm_ratio).
bool GreedySearch::isHarmonic(const Fit& fit) {
  if (para_.harm_ratio <= 0) {
    return false;
  }
  const Distr& distr = getDistr(fit.mono_mass);
  for (int f = 2; f <= 3; f++) {
    if (fit.charge % f != 0) {
      continue;
    }
    double on_obs = 0;
    double on_theo = 0;
    double off_obs = 0;
    double off_theo = 0;
    int off_num = 0;
    for (int i = fit.first; i <= fit.last; i++) {
      int idx = fit.peak_idxes[i];
      double obs = idx >= 0 ? residuals_[idx] : 0;
      if (((i - distr.apex) % f + f) % f == 0) {
        on_obs += obs;
        on_theo += distr.intes[i];
      } else {
        off_obs += obs;
        off_theo += distr.intes[i];
        off_num++;
      }
    }
    if (off_num < 2 || on_obs <= 0 || off_theo <= 0) {
      continue;
    }
    double ratio = (off_obs / off_theo) / (on_obs / on_theo);
    if (ratio < para_.harm_ratio) {
      return true;
    }
  }
  return false;
}

void GreedySearch::subtract(const FitVec& series) {
  for (const Fit& fit : series) {
    const Distr& distr = getDistr(fit.mono_mass);
    std::vector<int> whole;
    if (para_.dyn_bound != 0) {
      whole = getTailPeaks(fit);
    }
    for (int i = fit.first; i <= fit.last; i++) {
      int idx = fit.peak_idxes[i];
      if (idx < 0) {
        continue;
      }
      double theo = fit.scale * distr.intes[i];
      if (residuals_[idx] > para_.keep_ratio * theo) {
        residuals_[idx] -= theo;
      } else {
        residuals_[idx] = 0;
      }
    }
    for (int idx : whole) {
      residuals_[idx] = 0;
    }
  }
}

// The peaks the subtraction of a fit removes whole (dyn_bound): from the
// apex outwards on each side, the isotopes up to the boundary where the
// observed series stops (no peak) or rises again (residual > dyn_rise x the
// residual of the isotope before it); of those, the tail isotopes (relative
// intensity below core_ratio) and the isotopes beyond the fitted ones.
std::vector<int> GreedySearch::getTailPeaks(const Fit& fit) {
  std::vector<int> result;
  const Distr& distr = getDistr(fit.mono_mass);
  int n = static_cast<int>(distr.intes.size());
  int z = fit.charge;
  double mono_mz = fit.mono_mass / z + mass_constant::getProtonMass();
  double tolerance = env_para_ptr_->getMzTolerance(z);
  double iso_mass = mass_constant::getIsotopeMass();
  auto peakAt = [&](int i) {
    if (i >= fit.first && i <= fit.last) {
      return fit.peak_idxes[i];
    }
    double offset = (i >= 0 && i < n) ? distr.offsets[i]
                                      : distr.offsets[0] + i * iso_mass;
    int idx = findPeak(mono_mz + offset / z, tolerance);
    return (idx >= 0 && intes_[idx] >= min_inte_) ? idx : -1;
  };
  int apex_idx = fit.peak_idxes[distr.apex];
  if (apex_idx < 0) {
    return result;
  }
  for (int dir = -1; dir <= 1; dir += 2) {
    double prev = residuals_[apex_idx];
    int stop = dir > 0 ? fit.last + para_.dyn_max_ext
                       : fit.first - para_.dyn_max_ext;
    for (int i = distr.apex + dir; dir > 0 ? i <= stop : i >= stop; i += dir) {
      if (i < 0) {
        break;
      }
      int idx = peakAt(i);
      if (idx < 0 || residuals_[idx] <= 0 ||
          residuals_[idx] > para_.dyn_rise * prev) {
        break;
      }
      prev = residuals_[idx];
      bool tail = i < fit.first || i > fit.last ||
                  distr.intes[i] < para_.core_ratio;
      if (tail) {
        result.push_back(idx);
      }
    }
  }
  return result;
}

std::vector<FitVec> GreedySearch::run() {
  std::vector<FitVec> result;
  std::priority_queue<std::pair<double, int>> queue;
  for (std::size_t i = 0; i < residuals_.size(); i++) {
    if (residuals_[i] >= min_ref_inte_) {
      queue.emplace(residuals_[i], static_cast<int>(i));
    }
  }
  std::vector<int> seed_nums(residuals_.size(), 0);
  std::vector<bool> matched(residuals_.size(), false);
  double proton = mass_constant::getProtonMass();
  double iso_mass = mass_constant::getIsotopeMass();
  while (!queue.empty()) {
    auto [inte, seed] = queue.top();
    queue.pop();
    if (residuals_[seed] < min_ref_inte_ ||
        seed_nums[seed] >= para_.max_seed_num ||
        (para_.seed_matched == 0 && matched[seed])) {
      continue;
    }
    if (inte != residuals_[seed]) {
      // the intensity changed since the peak was queued
      queue.emplace(residuals_[seed], seed);
      continue;
    }
    seed_nums[seed]++;
    // 1. the charge: the series of every charge with the seed at its apex
    FitVec best_series;
    double best_score = 0;
    int best_charge = 0;
    for (int z = 1; z <= max_charge_; z++) {
      double ref_mass = (mzs_[seed] - proton) * z;
      if (ref_mass >= max_mass_ || ref_mass < env_para_ptr_->min_mass_) {
        continue;
      }
      if (para_.precheck != 0) {
        double tolerance = env_para_ptr_->getMzTolerance(z);
        int left = findPeak(mzs_[seed] - iso_mass / z, tolerance);
        int right = findPeak(mzs_[seed] + iso_mass / z, tolerance);
        if ((left < 0 || intes_[left] < min_inte_) &&
            (right < 0 || intes_[right] < min_inte_)) {
          continue;
        }
      }
      EnvPtr ref_env = EnvBase::getEnvByRefMass(ref_mass);
      if (ref_env == nullptr) {
        continue;
      }
      double mono_mass =
          ref_mass - (ref_env->getReferMz() - ref_env->getMonoMz());
      Fit fit = fitEnv(mono_mass, z);
      if (!passes(fit, para_.min_cos) || isHarmonic(fit)) {
        continue;
      }
      if (para_.min_explained_frac > 0 &&
          fit.explained < para_.min_explained_frac * fit.window_inte) {
        continue;
      }
      FitVec series = getSeries(mono_mass, z, fit);
      double score = seriesScore(series);
      if (debug_os_ != nullptr) {
        *debug_os_ << "  seed " << mzs_[seed] << " z " << z << " mass "
                   << mono_mass << " cos " << fit.cos << " charges "
                   << series.size() << " score " << score << "\n";
      }
      if (score > best_score) {
        best_score = score;
        best_series = std::move(series);
        best_charge = z;
      }
    }
    if (best_charge == 0) {
      continue;
    }
    // 2. the monoisotopic mass: the best series over the isotope offsets
    double base_mass = best_series[0].mono_mass;
    int aggr_k = para_.mono_by == 1 ? getAggregateOffset(best_series) : 0;
    for (int k = -para_.max_offset; k <= para_.max_offset; k++) {
      if (k == 0 || (para_.mono_by == 1 && k != aggr_k)) {
        continue;
      }
      double mono_mass = base_mass + k * iso_mass;
      Fit fit = fitEnv(mono_mass, best_charge);
      // the seed must stay matched in the seed charge's envelope
      bool has_seed = false;
      for (int i = fit.first; fit.valid && i <= fit.last; i++) {
        if (fit.peak_idxes[i] == seed) {
          has_seed = true;
        }
      }
      if (!has_seed || !passes(fit, para_.min_cos)) {
        continue;
      }
      FitVec series = getSeries(mono_mass, best_charge, fit);
      double score = seriesScore(series);
      if (score > best_score || para_.mono_by == 1) {
        best_score = score;
        best_series = std::move(series);
      }
    }
    if (debug_os_ != nullptr) {
      *debug_os_ << "mass " << best_series[0].mono_mass << " seed "
                 << mzs_[seed] << " inte " << inte << " charge "
                 << best_charge << " charges";
      for (const Fit& fit : best_series) {
        *debug_os_ << " " << fit.charge;
      }
      *debug_os_ << " score " << best_score << "\n";
    }
    int merge_into = -1;
    if (para_.merge_window > 0) {
      double best_diff = para_.merge_window;
      for (std::size_t k = 0; k < result.size(); k++) {
        double diff =
            std::abs(result[k][0].mono_mass - best_series[0].mono_mass);
        if (diff < best_diff) {
          best_diff = diff;
          merge_into = static_cast<int>(k);
        }
      }
    }
    if (merge_into >= 0) {
      double mass = result[merge_into][0].mono_mass;
      FitVec merged;
      for (const Fit& fit : best_series) {
        bool has_charge = false;
        for (const Fit& f : result[merge_into]) {
          has_charge = has_charge || f.charge == fit.charge;
        }
        if (has_charge) {
          continue;
        }
        Fit refit = fitEnv(mass, fit.charge);
        if (passes(refit, para_.nb_min_cos)) {
          merged.push_back(std::move(refit));
        }
      }
      if (merged.empty()) {
        continue;
      }
      best_series = std::move(merged);
    } else if (best_series.size() == 1 &&
               best_series[0].cos < para_.single_min_cos) {
      continue;
    }
    best_series = removeSharedEnvs(std::move(best_series));
    if (best_series.empty()) {
      continue;
    }
    // 3. subtract the envelopes of the mass at all its charges
    subtract(best_series);
    for (const Fit& fit : best_series) {
      found_masses_[fit.charge].push_back(fit.mono_mass);
      for (int i = fit.first; i <= fit.last; i++) {
        if (fit.peak_idxes[i] >= 0) {
          matched[fit.peak_idxes[i]] = true;
        }
      }
    }
    for (const Fit& fit : best_series) {
      for (int i = fit.first; i <= fit.last; i++) {
        int idx = fit.peak_idxes[i];
        if (idx >= 0 && residuals_[idx] >= min_ref_inte_) {
          queue.emplace(residuals_[idx], idx);
        }
      }
    }
    if (merge_into >= 0) {
      for (Fit& fit : best_series) {
        result[merge_into].push_back(std::move(fit));
      }
    } else {
      result.push_back(std::move(best_series));
    }
  }
  return result;
}

// The test of the shared peaks (TopfdPara::getMs1MaxSharedPeakNum, hidden
// option --ms1-max-shared-peak-num; toppic_test experiment 157): the
// envelopes of a mass, in the order of its charge series (the seed charge
// first; with merge_window the envelopes merged into an earlier mass), are
// built as they would be reported (getMatchEnv) and compared with every
// envelope reported before in the spectrum, those of the earlier masses and
// those of the mass kept before. The peaks two envelopes share are the
// spectrum peaks both match, as the DP counts them
// (CoTable::cntSharedPeakNum). An envelope that shares more than
// max_shared_peak_num_ peaks with one of them is not reported (nor
// subtracted); the others are kept and their peaks recorded. Not tested:
// the envelopes of the precursor windows (addPrecWinEnvs) and the PTCR
// search (deconvSeries, no limit).
FitVec GreedySearch::removeSharedEnvs(FitVec series) {
  if (max_shared_peak_num_ < 0) {
    return series;
  }
  FitVec kept;
  for (Fit& fit : series) {
    std::vector<int> peak_idxes;
    MatchEnvPtr env = getMatchEnv(peak_list_, fit, env_para_ptr_, min_inte_);
    if (env != nullptr) {
      ExpEnvPtr exp_env = env->getExpEnvPtr();
      for (int i = 0; i < exp_env->getPeakNum(); i++) {
        if (exp_env->isExist(i)) {
          peak_idxes.push_back(exp_env->getPeakIdx(i));
        }
      }
    }
    // the number of peaks shared with each envelope reported before
    std::unordered_map<int, int> shared_nums;
    int max_num = 0;
    int max_env = -1;
    for (int idx : peak_idxes) {
      for (int e : peak_envs_[idx]) {
        int num = ++shared_nums[e];
        if (num > max_num) {
          max_num = num;
          max_env = e;
        }
      }
    }
    if (max_num > max_shared_peak_num_) {
      if (debug_os_ != nullptr) {
        *debug_os_ << "shared mass " << fit.mono_mass << " charge "
                   << fit.charge << " shares " << max_num
                   << " peaks with mass " << reported_envs_[max_env].first
                   << " charge " << reported_envs_[max_env].second
                   << ": not reported\n";
      }
      continue;
    }
    int env_id = static_cast<int>(reported_envs_.size());
    for (int idx : peak_idxes) {
      if (peak_envs_[idx].empty() || peak_envs_[idx].back() != env_id) {
        peak_envs_[idx].push_back(env_id);
      }
    }
    reported_envs_.emplace_back(fit.mono_mass, fit.charge);
    kept.push_back(std::move(fit));
  }
  return kept;
}

// The envelope of the most intense peak in the m/z window [bgn, end], for a
// precursor window without an envelope (GreedyPara::prec_win_post). It is
// run on a search that has not subtracted anything, so on the original
// intensities. With one peak in the window, the peak is of charge 1. With
// more, every charge is tried as run() tries a seed (the averagine apex
// aligned with the peak, the series of the mass at the neighbouring
// charges), without the acceptance thresholds (isotope cosine, matched
// peaks, pre-check, single_min_cos) but with the harmonic test, and the
// best series score wins even when it is negative. The monoisotopic offset
// is chosen as in run() and kept only when the apex stays in the window.
// A peak below the noise level lowers the level to its intensity for this
// envelope, so that it is matched.
WinFit GreedySearch::fitWindow(double bgn, double end) {
  WinFit result;
  for (auto it = std::lower_bound(mzs_.begin(), mzs_.end(), bgn);
       it != mzs_.end() && *it <= end; ++it) {
    int i = static_cast<int>(it - mzs_.begin());
    if (intes_[i] <= 0) {
      continue;
    }
    result.peak_num++;
    if (result.seed < 0 || intes_[i] > intes_[result.seed]) {
      result.seed = i;
    }
  }
  if (result.seed < 0) {
    return result;
  }
  int seed = result.seed;
  double noise_inte = min_inte_;
  min_inte_ = std::min(min_inte_, intes_[seed]);
  double proton = mass_constant::getProtonMass();
  int max_charge = result.peak_num == 1 ? 1 : max_charge_;
  FitVec best_series;
  for (int z = 1; z <= max_charge; z++) {
    double ref_mass = (mzs_[seed] - proton) * z;
    if (ref_mass >= max_mass_ || ref_mass < env_para_ptr_->min_mass_) {
      continue;
    }
    EnvPtr ref_env = EnvBase::getEnvByRefMass(ref_mass);
    if (ref_env == nullptr) {
      continue;
    }
    double mono_mass =
        ref_mass - (ref_env->getReferMz() - ref_env->getMonoMz());
    Fit fit = fitEnv(mono_mass, z);
    if (!fit.valid || !fit.apex_match || isHarmonic(fit)) {
      continue;
    }
    FitVec series = getSeries(mono_mass, z, fit);
    double score = seriesScore(series);
    if (best_series.empty() || score > result.score) {
      result.score = score;
      best_series = std::move(series);
    }
  }
  if (!best_series.empty()) {
    result.fit = best_series[0];
    result.charge_num = static_cast<int>(best_series.size());
    int k = getAggregateOffset(best_series);
    int z = result.fit.charge;
    double mono_mass =
        result.fit.mono_mass + k * mass_constant::getIsotopeMass();
    Fit fit = k == 0 ? Fit() : fitEnv(mono_mass, z);
    bool has_seed = false;
    for (int i = fit.first; fit.valid && i <= fit.last; i++) {
      if (fit.peak_idxes[i] == seed) {
        has_seed = true;
      }
    }
    if (has_seed && fit.apex_match) {
      const Distr& distr = getDistr(mono_mass);
      double apex_mz = mono_mass / z + proton + distr.offsets[distr.apex] / z;
      if (apex_mz > bgn && apex_mz < end) {
        result.fit = std::move(fit);
        result.offset = k;
      }
    }
  }
  min_inte_ = noise_inte;
  return result;
}

// The envelope of a fit as main's MatchEnv: the EnvBase distribution of the
// mass at the charge, scaled to the fit and cut to the isotopes used,
// matched to the peaks of the spectrum.
MatchEnvPtr getMatchEnv(const PeakPtrVec& peak_list, const Fit& fit,
                        const EnvParaPtr& env_para_ptr, double min_inte) {
  EnvPtr theo_env = EnvBase::getEnvByMonoMass(fit.mono_mass, fit.charge);
  if (theo_env == nullptr) {
    return nullptr;
  }
  theo_env->changeIntensity(fit.scale / theo_env->getReferInte());
  int apex = theo_env->getReferIdx();
  theo_env = theo_env->getSubEnv(apex - fit.first, fit.last - apex);
  ExpEnvPtr exp_env = std::make_shared<ExpEnv>(
      peak_list, theo_env, env_para_ptr->getMzTolerance(fit.charge), min_inte);
  int mass_group = env_para_ptr->getMassGroup(fit.mono_mass);
  MatchEnvPtr env = std::make_shared<MatchEnv>(mass_group, theo_env, exp_env);
  env->compMsdeconvScr(env_para_ptr);
  return env;
}

// The window log (TOPFD_GREEDY_PREC_WIN_LOG=<file>): one line per precursor
// window, appended by the worker threads under a lock.
std::mutex prec_win_log_mutex;
bool prec_win_log_started = false;

void writePrecWinLog(const std::string& lines) {
  const char* name = std::getenv("TOPFD_GREEDY_PREC_WIN_LOG");
  if (name == nullptr || lines.empty()) {
    return;
  }
  std::lock_guard<std::mutex> lock(prec_win_log_mutex);
  std::ofstream os(name, std::ios::app);
  if (!prec_win_log_started) {
    os << "ms2_scan\twin_begin\twin_end\tstatus\tpeak_num\tseed_mz\t"
       << "seed_inte\tnoise_inte\tref_inte\tseed_in_env\tcharge\t"
       << "mono_mass\toffset\tcos\tmatch_num\tcharge_num\tscore\t"
       << "envcnn\tapex_mz\tapex_in_win\n";
    prec_win_log_started = true;
  }
  os << lines;
}

// The postprocessing step (GreedyPara::prec_win_post): every precursor
// window in which no envelope of envs has its reference (apex) m/z, the
// test of the MS1 feature detection's search for precursor windows without
// a feature (env_coll_detect, SeedEnv::getReferMz), gets the envelope of its
// most intense peak (GreedySearch::fitWindow), with its EnvCNN score and
// whatever the score; a window without a peak gets nothing. The envelopes
// are added to envs, so an envelope added for one window counts for the
// next windows.
void addPrecWinEnvs(const PeakPtrVec& peak_list,
                    const TopfdParaPtr& topfd_para_ptr,
                    const EnvParaPtr& env_para_ptr, double min_inte,
                    const std::vector<PrecWin>& prec_wins,
                    MatchEnvPtrVec& envs) {
  bool write_log = std::getenv("TOPFD_GREEDY_PREC_WIN_LOG") != nullptr;
  std::ostringstream log;
  log.precision(10);
  std::unique_ptr<GreedySearch> search;
  for (const PrecWin& win : prec_wins) {
    log << win.scan << "\t" << win.begin << "\t" << win.end << "\t";
    MatchEnvPtr win_env = nullptr;
    for (const MatchEnvPtr& env : envs) {
      double mz = env->getTheoEnvPtr()->getReferMz();
      if (mz > win.begin && mz < win.end) {
        win_env = env;
        break;
      }
    }
    if (win_env != nullptr) {
      EnvPtr theo = win_env->getTheoEnvPtr();
      log << "env\t\t\t\t\t\t\t" << theo->getCharge() << "\t"
          << theo->getMonoNeutralMass() << "\t\t\t\t\t\t"
          << win_env->getEnvcnnScore() << "\t" << theo->getReferMz()
          << "\t1\n";
      continue;
    }
    if (search == nullptr) {
      search = std::make_unique<GreedySearch>(
          peak_list, topfd_para_ptr, env_para_ptr,
          topfd_para_ptr->getMsOneSnRatio(), topfd_para_ptr->getMaxMass(),
          topfd_para_ptr->getMaxCharge(), ChargeBound(), -1, nullptr);
    }
    WinFit win_fit = search->fitWindow(win.begin, win.end);
    if (win_fit.seed < 0) {
      log << "no_peak\t0\n";
      continue;
    }
    MatchEnvPtr env = nullptr;
    if (win_fit.fit.valid) {
      double seed_inte = peak_list[win_fit.seed]->getIntensity();
      env = getMatchEnv(peak_list, win_fit.fit, env_para_ptr,
                        std::min(min_inte, seed_inte));
    }
    if (env != nullptr) {
      MatchEnvPtrVec win_envs = {env};
      if (getPara().envcnn != 0) {
        onnx_env_cnn::computeEnvScores(peak_list, win_envs);
      }
    }
    if (write_log) {
      // whether the seed is matched in an envelope reported before
      int seed_in_env = 0;
      for (const MatchEnvPtr& e : envs) {
        ExpEnvPtr exp_env = e->getExpEnvPtr();
        for (int p = 0; p < exp_env->getPeakNum(); p++) {
          if (exp_env->isExist(p) && exp_env->getPeakIdx(p) == win_fit.seed) {
            seed_in_env = 1;
          }
        }
      }
      std::vector<double> intes;
      for (const PeakPtr& peak : peak_list) {
        intes.push_back(peak->getIntensity());
      }
      double noise_inte = baseline_util::getBaseLine(intes);
      log << (env != nullptr ? "new" : "none") << "\t" << win_fit.peak_num
          << "\t" << peak_list[win_fit.seed]->getPosition() << "\t"
          << peak_list[win_fit.seed]->getIntensity() << "\t" << noise_inte
          << "\t" << noise_inte * topfd_para_ptr->getMsOneSnRatio() << "\t"
          << seed_in_env << "\t";
      if (env != nullptr) {
        const Fit& fit = win_fit.fit;
        double apex_mz = env->getTheoEnvPtr()->getReferMz();
        log << fit.charge << "\t" << fit.mono_mass << "\t" << win_fit.offset
            << "\t" << fit.cos << "\t" << fit.match_num << "\t"
            << win_fit.charge_num << "\t" << win_fit.score << "\t"
            << env->getEnvcnnScore() << "\t" << apex_mz << "\t"
            << (apex_mz > win.begin && apex_mz < win.end ? 1 : 0);
      }
      log << "\n";
    }
    if (env != nullptr) {
      envs.push_back(env);
    }
  }
  if (write_log) {
    writePrecWinLog(log.str());
  }
}

}  // namespace

bool isPrecWinFirst() { return getPara().prec_win != 0; }

MatchEnvPtrVec deconv(const PeakPtrVec& peak_list,
                      const TopfdParaPtr& topfd_para_ptr,
                      const std::vector<PrecWin>& prec_wins,
                      const std::string& debug_name) {
  MatchEnvPtrVec envs;
  if (peak_list.empty()) {
    return envs;
  }
  // the envelope parameters of DeconvSingleSp for MS1
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>(
      topfd_para_ptr->getMzError(), topfd_para_ptr->getMaxMissPeakNum(),
      topfd_para_ptr->isRunFilterByMz(),
      topfd_para_ptr->getMassGroupBoundary());
  std::unique_ptr<std::ofstream> debug_os;
  if (!debug_name.empty()) {
    debug_os = std::make_unique<std::ofstream>(debug_name);
  }
  GreedySearch search(peak_list, topfd_para_ptr, env_para_ptr,
                      topfd_para_ptr->getMsOneSnRatio(),
                      topfd_para_ptr->getMaxMass(),
                      topfd_para_ptr->getMaxCharge(), ChargeBound(),
                      topfd_para_ptr->getMs1MaxSharedPeakNum(),
                      debug_os.get());
  std::vector<FitVec> masses = search.run();

  std::vector<double> intes;
  intes.reserve(peak_list.size());
  for (const PeakPtr& peak : peak_list) {
    intes.push_back(peak->getIntensity());
  }
  double noise_inte = baseline_util::getBaseLine(intes);
  double sn_ratio = topfd_para_ptr->getMsOneSnRatio();
  double min_inte = sn_ratio < 1.0 ? noise_inte * sn_ratio : noise_inte;
  for (const FitVec& series : masses) {
    for (const Fit& fit : series) {
      MatchEnvPtr env = getMatchEnv(peak_list, fit, env_para_ptr, min_inte);
      if (env != nullptr) {
        envs.push_back(env);
      }
    }
  }
  // as DeconvSingleSp::postprocess
  PeakPtrVec peaks = peak_list;
  match_env_util::assignIntensity(peaks, envs);
  if (getPara().envcnn != 0) {
    onnx_env_cnn::computeEnvScores(peak_list, envs);
  }
  if (getPara().min_envcnn > 0) {
    MatchEnvPtrVec kept;
    for (const MatchEnvPtr& env : envs) {
      if (env->getEnvcnnScore() >= getPara().min_envcnn) {
        kept.push_back(env);
      }
    }
    envs = std::move(kept);
  }
  if (getPara().prec_win_post != 0 && !prec_wins.empty()) {
    addPrecWinEnvs(peak_list, topfd_para_ptr, env_para_ptr, min_inte,
                   prec_wins, envs);
  }
  if (topfd_para_ptr->isSortUseMsDeconv()) {
    std::sort(envs.begin(), envs.end(), MatchEnv::cmpMsdeconvScoreDec);
  } else {
    std::sort(envs.begin(), envs.end(), MatchEnv::cmpEnvcnnScoreDec);
  }
  match_env_util::assignIntensity(peaks, envs);
  return envs;
}

std::vector<MatchEnvPtrVec> deconvSeries(const PeakPtrVec& peak_list,
                                         const TopfdParaPtr& topfd_para_ptr,
                                         double sn_ratio, double max_mass,
                                         int max_charge,
                                         const ChargeBound& bound,
                                         double min_envcnn,
                                         const std::string& debug_name) {
  std::vector<MatchEnvPtrVec> result;
  if (peak_list.empty()) {
    return result;
  }
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>(
      topfd_para_ptr->getMzError(), topfd_para_ptr->getMaxMissPeakNum(),
      topfd_para_ptr->isRunFilterByMz(),
      topfd_para_ptr->getMassGroupBoundary());
  std::unique_ptr<std::ofstream> debug_os;
  if (!debug_name.empty()) {
    debug_os = std::make_unique<std::ofstream>(debug_name);
  }
  GreedySearch search(peak_list, topfd_para_ptr, env_para_ptr, sn_ratio,
                      max_mass, max_charge, bound, -1, debug_os.get());
  std::vector<FitVec> masses = search.run();

  std::vector<double> intes;
  intes.reserve(peak_list.size());
  for (const PeakPtr& peak : peak_list) {
    intes.push_back(peak->getIntensity());
  }
  double noise_inte = baseline_util::getBaseLine(intes);
  double min_inte = sn_ratio < 1.0 ? noise_inte * sn_ratio : noise_inte;
  // the envelopes of all masses, with the mass each belongs to
  MatchEnvPtrVec envs;
  std::vector<int> owners;
  for (std::size_t k = 0; k < masses.size(); k++) {
    for (const Fit& fit : masses[k]) {
      MatchEnvPtr env = getMatchEnv(peak_list, fit, env_para_ptr, min_inte);
      if (env != nullptr) {
        envs.push_back(env);
        owners.push_back(static_cast<int>(k));
      }
    }
  }
  // the intensities shared between envelopes and the EnvCNN scores, as in
  // deconv; the envelopes below min_envcnn are dropped and the intensities
  // assigned again among the kept ones
  PeakPtrVec peaks = peak_list;
  match_env_util::assignIntensity(peaks, envs);
  if (getPara().envcnn != 0) {
    onnx_env_cnn::computeEnvScores(peak_list, envs);
  }
  MatchEnvPtrVec kept;
  result.resize(masses.size());
  for (std::size_t i = 0; i < envs.size(); i++) {
    if (envs[i]->getEnvcnnScore() >= min_envcnn) {
      kept.push_back(envs[i]);
      result[owners[i]].push_back(envs[i]);
    }
  }
  peaks = peak_list;
  match_env_util::assignIntensity(peaks, kept);
  result.erase(std::remove_if(result.begin(), result.end(),
                              [](const MatchEnvPtrVec& v) { return v.empty(); }),
               result.end());
  return result;
}

}  // namespace greedy_ms1_deconv

}  // namespace toppic
