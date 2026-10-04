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

#include "topfd/deconv/greedy_ms2_deconv.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/exp_env.hpp"
#include "ms/env/match_env_filter.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

namespace greedy_ms2_deconv {

namespace {

// The parameters of the search. The search started as a copy of
// greedy_ms1_deconv's (toppic_test experiment 104) with that experiment's
// defaults; the MS/MS values are from toppic_test experiment 106. They can be
// changed for research with the environment variable TOPFD_GREEDY_MS2_PARA, a
// comma-separated list of name=value pairs (e.g. "min_cos=0.7,keep_ratio=2").
struct GreedyMs2Para {
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
  // squared isotope cosine (1), or (2) as 1 with the theoretical intensity
  // a matched peak falls short of also subtracted (2 found fewer fragments
  // in experiment 106)
  int score_mode = 1;
  // monoisotopic offsets (isotopes) searched around the apex alignment
  int max_offset = 2;
  // times a peak may be the seed
  int max_seed_num = 2;
  // the monoisotopic mass: the offset whose series scores best (0), or the
  // offset whose averagine distribution has the best cosine with the
  // isotope intensities summed over the series' charges (1; FLASHDeconv's
  // getIsotopeCosineAndIsoOffset; better for MS1, worse for fragments in
  // toppic_test experiment 106)
  int mono_by = 0;
  // the envelope boundary of the subtraction found from the observed peaks
  // (1) or the isotopes fitted (0), as in greedy_ms1_deconv
  int dyn_bound = 1;
  double dyn_rise = 1.2;
  int dyn_max_ext = 6;
  // a charge z is tried for a seed only when a peak lies one isotope
  // spacing (getIsotopeMass() / z) from it on either side (1), also at the
  // averagine distribution's own spacing around its apex (2), or always (0)
  int precheck = 2;
  // harmonic test of the seed charge z (0 = none), as in greedy_ms1_deconv;
  // off for MS/MS: with 0.4 it cost fragments in toppic_test experiment 106
  double harm_ratio = 0;
  // a mass found at one charge only is reported (and subtracted) only when
  // its envelope's isotope cosine reaches single_min_cos (0 = no test); most
  // fragments are found at one charge only, and MS1's 0.8 lost some of them
  double single_min_cos = 0.7;
  // envelopes with an EnvCNN score below min_envcnn are not reported (0 =
  // all reported, as main's MS/MS cutoff; TopPIC filters with -F)
  double min_envcnn = 0;
  // an envelope within dup_window Da of a mass already found at the same
  // charge is not used (0 = no test); main's MS/MS filtering-by-mz window
  // is 2 Da
  double dup_window = 2;
  // with dup_window: only a mass an integer number of isotope spacings
  // (within dup_ppm ppm) from the earlier one is a duplicate (1), or any
  // mass inside the window (0). Distinct masses of one charge can lie
  // closer than 2 Da (in one 480_F3_02 spectrum the b/y fragments 1,951.94
  // Da at charge 2 and 627.33 Da at charge 1 lie 1.93 and 0.95 Da from more
  // intense masses of the same charge); an isotope-shifted copy of a mass
  // is k x getIsotopeMass() away.
  int dup_iso = 1;
  double dup_ppm = 20;
  // every charge competes with the best of its isotope offsets (1), or with
  // the apex alignment only, the offset being chosen after the charge (0,
  // greedy_ms1_deconv). Below about 5,000 Da the monoisotopic peak and its
  // neighbour are of similar height, so the apex alignment of the right
  // charge is often one isotope off, fits badly and loses to a wrong charge
  // (a 1,823.84 Da charge 2 fragment read as 911.92 Da at charge 1).
  int charge_offsets = 1;
  // the m/z tolerance of the search in ppm (0 = TopFD's m/z error, -e,
  // 0.02 by default). On 480_F3_02 the peaks of envelopes with an EnvCNN
  // score >= 0.9 lie within 0.6 ppm of their isotopes (median; 4.1 ppm at
  // the 90th percentile), so 0.02 m/z (25-40 ppm) lets a high-charge
  // envelope collect the peaks of other species. The envelopes reported are
  // still matched to the peaks with TopFD's tolerance (for EnvCNN).
  double ppm_tol = 5;
  // with ppm_tol: an isotope without a peak within ppm_tol but with one
  // within TopFD's m/z error is not counted as missing (1; no penalty, no
  // explained intensity), or it is (0). Overlapping envelopes' peaks merge
  // into one centroid between their two positions, several ppm from both.
  int near_no_pen = 1;
  // how the isotope offset of a charge is chosen with charge_offsets: by
  // the score of the whole charge series at each offset (0), or by the
  // score of the seed charge's envelope alone (1), the series being built
  // at that offset. A neighbouring charge can be the envelope of another
  // mass one isotope away, which then decides the offset of the seed.
  int offset_by_seed = 1;
  // a neighbouring charge joins a series only when its envelope fits the
  // series' mass better than the masses one isotope lower and higher (1),
  // or whenever it passes (0)
  int nb_consistent = 1;
  // rounds of the refit pass after the greedy search (0 = none): every mass
  // is fitted again at +-1 isotope on the spectrum without the envelopes of
  // all the other masses found, and moved when its seed charge's envelope
  // scores better there. The greedy order decides a mass before the
  // weaker species whose peaks lie in its envelope, which then lend
  // intensity to a wrong isotope offset (a 2,538.39 Da charge 4 fragment
  // read as 2,537.39 Da through a charge 5 species' peak at its M-1). Off:
  // it moved few masses, and with alt_ratio 0.9 and miss_weight 0.5 it gave
  // no more identifications at 1% FDR (toppic_test experiment 106, run g4).
  int refit_rounds = 0;
  // a mass below alt_max_mass whose seed charge's envelope one isotope
  // lower or higher scores at least alt_ratio x the chosen one is also
  // reported at that mass (0 = never). TopPIC tolerates a +-1 isotope error
  // only above 5,000 Da (its fragment mass extension); below, about 11% of
  // the fragments found are one isotope off (toppic_test experiment 106).
  // TopFD's former addMultiMass reported such alternatives too. Off: 0.9
  // recovered fragments but added masses, and gave no more identifications
  // at 1% FDR (experiment 106, run g4).
  double alt_ratio = 0;
  double alt_max_mass = 5000;
};

GreedyMs2Para readPara() {
  GreedyMs2Para para;
  const char* str = std::getenv("TOPFD_GREEDY_MS2_PARA");
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
    } else if (name == "mono_by") {
      para.mono_by = static_cast<int>(value);
    } else if (name == "dyn_bound") {
      para.dyn_bound = static_cast<int>(value);
    } else if (name == "dyn_rise") {
      para.dyn_rise = value;
    } else if (name == "dyn_max_ext") {
      para.dyn_max_ext = static_cast<int>(value);
    } else if (name == "precheck") {
      para.precheck = static_cast<int>(value);
    } else if (name == "harm_ratio") {
      para.harm_ratio = value;
    } else if (name == "single_min_cos") {
      para.single_min_cos = value;
    } else if (name == "min_envcnn") {
      para.min_envcnn = value;
    } else if (name == "dup_window") {
      para.dup_window = value;
    } else if (name == "dup_iso") {
      para.dup_iso = static_cast<int>(value);
    } else if (name == "dup_ppm") {
      para.dup_ppm = value;
    } else if (name == "offset_by_seed") {
      para.offset_by_seed = static_cast<int>(value);
    } else if (name == "alt_ratio") {
      para.alt_ratio = value;
    } else if (name == "alt_max_mass") {
      para.alt_max_mass = value;
    } else if (name == "refit_rounds") {
      para.refit_rounds = static_cast<int>(value);
    } else if (name == "nb_consistent") {
      para.nb_consistent = static_cast<int>(value);
    } else if (name == "near_no_pen") {
      para.near_no_pen = static_cast<int>(value);
    } else if (name == "ppm_tol") {
      para.ppm_tol = value;
    } else if (name == "charge_offsets") {
      para.charge_offsets = static_cast<int>(value);
    } else {
      LOG_WARN("Unknown TOPFD_GREEDY_MS2_PARA entry " << name);
    }
  }
  return para;
}

const GreedyMs2Para& getPara() {
  static const GreedyMs2Para para = readPara();
  return para;
}

// The averagine distribution of a mass bin of EnvBase: per isotope, the mass
// offset from the monoisotopic peak and the intensity relative to the apex.
struct Distr {
  std::vector<double> offsets;
  std::vector<double> intes;
  int apex = 0;
  // the isotopes above GreedyMs2Para::min_rel_inte
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
  // the theoretical intensity matched peaks fall short of
  double shortfall = 0;
  bool apex_match = false;
  // the isotopes used: [first, last] of the distribution
  int first = 0;
  int last = -1;
  // per isotope of the distribution: the matched peak or -1
  std::vector<int> peak_idxes;
  double score() const {
    double score = explained - getPara().miss_weight * penalty;
    if (getPara().score_mode == 2) {
      score -= getPara().miss_weight * shortfall;
    }
    return getPara().score_mode >= 1 ? score * cos * cos : score;
  }
};

using FitVec = std::vector<Fit>;

// One spectrum's search.
class GreedySearch {
 public:
  GreedySearch(const PeakPtrVec& peak_list, const TopfdParaPtr& topfd_para_ptr,
               const EnvParaPtr& env_para_ptr, double max_mass,
               int max_charge, std::ostream* debug_os);

  // the charge series of the masses found, in the order they were found
  std::vector<FitVec> run();

  double getMinInte() const { return min_inte_; }

  // the alternative masses (GreedyMs2Para::alt_ratio)
  const FitVec& getAlts() const { return alts_; }

 private:
  const GreedyMs2Para& para_;
  EnvParaPtr env_para_ptr_;
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
  // the alternative masses one isotope away (alt_ratio)
  FitVec alts_;

  const Distr& getDistr(double mono_mass);

  int findPeak(double mz, double tolerance) const;

  // the matching tolerance of the search at an m/z and charge
  double getTolerance(double mz, int charge) const;

  Fit fitEnv(double mono_mass, int charge);

  bool passes(const Fit& fit, double min_cos) const;

  static bool hasPeak(const Fit& fit, int peak_idx);

  bool hasNeighbourPeak(int seed, int charge, double ref_mass);

  bool isBestOffset(const Fit& fit);

  FitVec getSeries(double mono_mass, int charge, const Fit& seed_fit);

  static double seriesScore(const FitVec& series);

  int getAggregateOffset(const FitVec& series);

  std::vector<int> getTailPeaks(const Fit& fit);

  bool isHarmonic(const Fit& fit);

  void subtract(const FitVec& series);

  void addContrib(const Fit& fit, std::vector<double>& contrib,
                  double sign);

  void refit(std::vector<FitVec>& result);
};

GreedySearch::GreedySearch(const PeakPtrVec& peak_list,
                           const TopfdParaPtr& topfd_para_ptr,
                           const EnvParaPtr& env_para_ptr, double max_mass,
                           int max_charge, std::ostream* debug_os)
    : para_(getPara()),
      env_para_ptr_(env_para_ptr),
      debug_os_(debug_os),
      max_charge_(max_charge),
      max_mass_(max_mass) {
  mzs_.reserve(peak_list.size());
  intes_.reserve(peak_list.size());
  for (const PeakPtr& peak : peak_list) {
    mzs_.push_back(peak->getPosition());
    intes_.push_back(peak->getIntensity());
  }
  residuals_ = intes_;
  found_masses_.resize(max_charge_ + 1);
  // the noise level and the thresholds of DeconvData::initMinInte with the
  // MS/MS signal-to-noise ratio
  double noise_inte = baseline_util::getBaseLine(intes_);
  double sn_ratio = topfd_para_ptr->getMsTwoSnRatio();
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

double GreedySearch::getTolerance(double mz, int charge) const {
  if (para_.ppm_tol > 0) {
    return mz * para_.ppm_tol * 1e-6;
  }
  return env_para_ptr_->getMzTolerance(charge);
}

Fit GreedySearch::fitEnv(double mono_mass, int charge) {
  Fit fit;
  fit.mono_mass = mono_mass;
  fit.charge = charge;
  if (mono_mass < env_para_ptr_->min_mass_ || mono_mass >= max_mass_) {
    return fit;
  }
  const Distr& distr = getDistr(mono_mass);
  if (distr.last < distr.first) {
    return fit;
  }
  double mono_mz = mono_mass / charge + mass_constant::getProtonMass();
  double tolerance = getTolerance(mono_mz, charge);
  int n = static_cast<int>(distr.intes.size());
  fit.peak_idxes.assign(n, -1);
  std::vector<double> obs(n, 0.0);
  // the isotopes with a peak only within TopFD's m/z error (near_no_pen)
  std::vector<bool> near(n, false);
  for (int i = distr.first; i <= distr.last; i++) {
    int idx = findPeak(mono_mz + distr.offsets[i] / charge, tolerance);
    if (idx < 0 && para_.near_no_pen != 0 && para_.ppm_tol > 0) {
      int near_idx = findPeak(mono_mz + distr.offsets[i] / charge,
                              env_para_ptr_->getMzTolerance(charge));
      near[i] = near_idx >= 0 && intes_[near_idx] >= min_inte_;
    }
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
    if (fit.peak_idxes[i] < 0 && !near[i]) {
      fit.miss_num++;
      fit.penalty += theo;
    } else if (obs[i] > 0) {
      fit.match_num++;
      fit.explained += std::min(obs[i], theo);
      fit.shortfall += std::max(theo - obs[i], 0.0);
    }
  }
  fit.apex_match = obs[distr.apex] > 0;
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
    double iso_mass = mass_constant::getIsotopeMass();
    double tole = fit.mono_mass * para_.dup_ppm * 1e-6;
    for (double mass : found_masses_[fit.charge]) {
      double diff = std::abs(mass - fit.mono_mass);
      if (diff >= para_.dup_window) {
        continue;
      }
      if (para_.dup_iso == 0 ||
          std::abs(diff - std::round(diff / iso_mass) * iso_mass) <= tole) {
        return false;
      }
    }
  }
  return true;
}

// The charge pre-check (GreedyMs2Para::precheck): whether a peak lies one
// isotope spacing from the seed at the charge, on either side. The spacing
// is getIsotopeMass() or, with precheck 2, also the spacing of the
// averagine distribution around its apex at the seed's mass: the M+1 peak
// of averagine is 1.00286 Da from its monoisotopic peak at every mass
// (toppic_test experiment 107), 3 ppm from the constant at m/z 203
// and charge 1, so that a few ppm of measurement error put it outside a
// 5 ppm window (6 ppm in a 480_F3_02 spectrum).
bool GreedySearch::hasNeighbourPeak(int seed, int charge, double ref_mass) {
  double mz = mzs_[seed];
  double tolerance = getTolerance(mz, charge);
  std::vector<double> spacings = {mass_constant::getIsotopeMass()};
  if (para_.precheck == 2) {
    const Distr& distr = getDistr(ref_mass);
    int n = static_cast<int>(distr.offsets.size());
    if (distr.apex + 1 < n) {
      spacings.push_back(distr.offsets[distr.apex + 1] -
                         distr.offsets[distr.apex]);
    }
    if (distr.apex >= 1) {
      spacings.push_back(distr.offsets[distr.apex] -
                         distr.offsets[distr.apex - 1]);
    }
  }
  for (double spacing : spacings) {
    for (int dir = -1; dir <= 1; dir += 2) {
      int idx = findPeak(mz + dir * spacing / charge, tolerance);
      if (idx >= 0 && intes_[idx] >= min_inte_) {
        return true;
      }
    }
  }
  return false;
}

// whether the peak is matched in the isotopes the fit uses
bool GreedySearch::hasPeak(const Fit& fit, int peak_idx) {
  for (int i = fit.first; fit.valid && i <= fit.last; i++) {
    if (fit.peak_idxes[i] == peak_idx) {
      return true;
    }
  }
  return false;
}

// whether the fit's mass fits its charge's peaks better than the masses one
// isotope lower and higher (nb_consistent)
bool GreedySearch::isBestOffset(const Fit& fit) {
  double iso_mass = mass_constant::getIsotopeMass();
  double score = fit.score();
  for (int k = -1; k <= 1; k += 2) {
    Fit alt = fitEnv(fit.mono_mass + k * iso_mass, fit.charge);
    if (alt.valid && alt.apex_match && alt.score() > score) {
      return false;
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
      if (passes(fit, para_.nb_min_cos) &&
          (para_.nb_consistent == 0 || isBestOffset(fit))) {
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
// series' charges (greedy_ms1_deconv's getAggregateOffset).
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
    double tolerance = getTolerance(mono_mz, z);
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
// f = 2 or 3 (greedy_ms1_deconv's isHarmonic).
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

// The peaks the subtraction of a fit removes whole (dyn_bound), as
// greedy_ms1_deconv's getTailPeaks: from the apex outwards, the isotopes up
// to where the observed series stops or rises again; of those, the tail
// isotopes and the isotopes beyond the fitted ones.
std::vector<int> GreedySearch::getTailPeaks(const Fit& fit) {
  std::vector<int> result;
  const Distr& distr = getDistr(fit.mono_mass);
  int n = static_cast<int>(distr.intes.size());
  int z = fit.charge;
  double mono_mz = fit.mono_mass / z + mass_constant::getProtonMass();
  double tolerance = getTolerance(mono_mz, z);
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
    if (residuals_[i] >= min_ref_inte_ && residuals_[i] > 0) {
      queue.emplace(residuals_[i], static_cast<int>(i));
    }
  }
  std::vector<int> seed_nums(residuals_.size(), 0);
  double proton = mass_constant::getProtonMass();
  double iso_mass = mass_constant::getIsotopeMass();
  while (!queue.empty()) {
    auto [inte, seed] = queue.top();
    queue.pop();
    if (residuals_[seed] < min_ref_inte_ || residuals_[seed] <= 0 ||
        seed_nums[seed] >= para_.max_seed_num) {
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
      if (para_.precheck != 0 && !hasNeighbourPeak(seed, z, ref_mass)) {
        continue;
      }
      EnvPtr ref_env = EnvBase::getEnvByRefMass(ref_mass);
      if (ref_env == nullptr) {
        continue;
      }
      double apex_mass =
          ref_mass - (ref_env->getReferMz() - ref_env->getMonoMz());
      int k_max = para_.charge_offsets != 0 ? para_.max_offset : 0;
      // offset_by_seed: only the offset whose seed envelope scores best
      double best_seed_score = 0;
      int best_seed_k = 0;
      if (para_.offset_by_seed != 0) {
        for (int k = -k_max; k <= k_max; k++) {
          Fit fit = fitEnv(apex_mass + k * iso_mass, z);
          if (passes(fit, para_.min_cos) && hasPeak(fit, seed) &&
              fit.score() > best_seed_score) {
            best_seed_score = fit.score();
            best_seed_k = k;
          }
        }
      }
      for (int k = -k_max; k <= k_max; k++) {
        if (para_.offset_by_seed != 0 && k != best_seed_k) {
          continue;
        }
        double mono_mass = apex_mass + k * iso_mass;
        Fit fit = fitEnv(mono_mass, z);
        if (!passes(fit, para_.min_cos) || !hasPeak(fit, seed) ||
            isHarmonic(fit)) {
          continue;
        }
        FitVec series = getSeries(mono_mass, z, fit);
        double score = seriesScore(series);
        if (debug_os_ != nullptr) {
          *debug_os_ << "  seed " << mzs_[seed] << " z " << z << " mass "
                     << mono_mass << " cos " << fit.cos << " charges "
                     << series.size() << " score " << score << " |";
          for (const Fit& f : series) {
            *debug_os_ << " z" << f.charge << ":cos=" << f.cos
                       << ",expl=" << f.explained << ",pen=" << f.penalty
                       << ",m=" << f.match_num;
          }
          *debug_os_ << "\n";
        }
        if (score > best_score) {
          best_score = score;
          best_series = std::move(series);
          best_charge = z;
        }
      }
    }
    if (best_charge == 0) {
      continue;
    }
    // 2. the monoisotopic mass: the best series over the isotope offsets
    // (already compared with charge_offsets unless mono_by is 1)
    double base_mass = best_series[0].mono_mass;
    int aggr_k = para_.mono_by == 1 ? getAggregateOffset(best_series) : 0;
    for (int k = -para_.max_offset; k <= para_.max_offset; k++) {
      if (k == 0 || (para_.mono_by == 1 && k != aggr_k) ||
          (para_.mono_by == 0 && para_.charge_offsets != 0)) {
        continue;
      }
      double mono_mass = base_mass + k * iso_mass;
      Fit fit = fitEnv(mono_mass, best_charge);
      // the seed must stay matched in the seed charge's envelope
      if (!hasPeak(fit, seed) || !passes(fit, para_.min_cos)) {
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
    if (best_series.size() == 1 &&
        best_series[0].cos < para_.single_min_cos) {
      continue;
    }
    // the alternative mass one isotope away (alt_ratio)
    if (para_.alt_ratio > 0 && best_series[0].mono_mass < para_.alt_max_mass) {
      const Fit& chosen = best_series[0];
      double min_score = para_.alt_ratio * chosen.score();
      Fit best_alt;
      for (int k = -1; k <= 1; k += 2) {
        Fit alt = fitEnv(chosen.mono_mass + k * iso_mass, chosen.charge);
        if (passes(alt, para_.min_cos) && alt.score() >= min_score &&
            (!best_alt.valid || alt.score() > best_alt.score())) {
          best_alt = std::move(alt);
        }
      }
      if (best_alt.valid) {
        alts_.push_back(std::move(best_alt));
      }
    }
    // 3. subtract the envelopes of the mass at all its charges
    subtract(best_series);
    for (const Fit& fit : best_series) {
      found_masses_[fit.charge].push_back(fit.mono_mass);
    }
    for (const Fit& fit : best_series) {
      for (int i = fit.first; i <= fit.last; i++) {
        int idx = fit.peak_idxes[i];
        if (idx >= 0 && residuals_[idx] >= min_ref_inte_ &&
            residuals_[idx] > 0) {
          queue.emplace(residuals_[idx], idx);
        }
      }
    }
    result.push_back(std::move(best_series));
  }
  refit(result);
  return result;
}

// adds (sign 1) or removes (sign -1) the theoretical intensities of a fit's
// matched isotopes to or from their peaks
void GreedySearch::addContrib(const Fit& fit, std::vector<double>& contrib,
                              double sign) {
  const Distr& distr = getDistr(fit.mono_mass);
  for (int i = fit.first; i <= fit.last; i++) {
    int idx = fit.peak_idxes[i];
    if (idx >= 0) {
      contrib[idx] += sign * fit.scale * distr.intes[i];
    }
  }
}

// The refit pass (GreedyMs2Para::refit_rounds).
void GreedySearch::refit(std::vector<FitVec>& result) {
  double iso_mass = mass_constant::getIsotopeMass();
  auto eraseMass = [this](int charge, double mass) {
    std::vector<double>& masses = found_masses_[charge];
    auto it = std::find(masses.begin(), masses.end(), mass);
    if (it != masses.end()) {
      masses.erase(it);
    }
  };
  for (int round = 0; round < para_.refit_rounds; round++) {
    std::vector<double> contrib(intes_.size(), 0.0);
    for (const FitVec& series : result) {
      for (const Fit& fit : series) {
        addContrib(fit, contrib, 1);
      }
    }
    int changed = 0;
    for (FitVec& series : result) {
      for (const Fit& fit : series) {
        addContrib(fit, contrib, -1);
        eraseMass(fit.charge, fit.mono_mass);
      }
      // the spectrum without the other masses' envelopes
      for (std::size_t p = 0; p < intes_.size(); p++) {
        residuals_[p] = std::max(intes_[p] - contrib[p], 0.0);
      }
      double mass = series[0].mono_mass;
      int charge = series[0].charge;
      Fit cur = fitEnv(mass, charge);
      double best_score = cur.valid ? cur.score() : 0;
      int best_k = 0;
      for (int k = -1; k <= 1; k += 2) {
        Fit alt = fitEnv(mass + k * iso_mass, charge);
        if (passes(alt, para_.min_cos) && alt.score() > best_score) {
          best_score = alt.score();
          best_k = k;
        }
      }
      if (best_k != 0) {
        double new_mass = mass + best_k * iso_mass;
        FitVec moved;
        for (const Fit& fit : series) {
          Fit alt = fitEnv(new_mass, fit.charge);
          double min_cos =
              fit.charge == charge ? para_.min_cos : para_.nb_min_cos;
          if (passes(alt, min_cos)) {
            moved.push_back(std::move(alt));
          }
        }
        if (!moved.empty() && moved[0].charge == charge) {
          if (debug_os_ != nullptr) {
            *debug_os_ << "refit mass " << mass << " charge " << charge
                       << " -> " << new_mass << "\n";
          }
          series = std::move(moved);
          changed++;
        }
      }
      for (const Fit& fit : series) {
        addContrib(fit, contrib, 1);
        found_masses_[fit.charge].push_back(fit.mono_mass);
      }
    }
    if (changed == 0) {
      break;
    }
  }
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

}  // namespace

MatchEnvPtrVec deconv(const PeakPtrVec& peak_list,
                      const TopfdParaPtr& topfd_para_ptr, double max_mass,
                      int max_charge, const std::string& debug_name) {
  MatchEnvPtrVec envs;
  if (peak_list.empty()) {
    return envs;
  }
  // the envelope parameters of DeconvSingleSp
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>(
      topfd_para_ptr->getMzError(), topfd_para_ptr->getMaxMissPeakNum(),
      topfd_para_ptr->isRunFilterByMz(),
      topfd_para_ptr->getMassGroupBoundary());
  std::unique_ptr<std::ofstream> debug_os;
  if (!debug_name.empty()) {
    debug_os = std::make_unique<std::ofstream>(debug_name);
  }
  GreedySearch search(peak_list, topfd_para_ptr, env_para_ptr, max_mass,
                      max_charge, debug_os.get());
  std::vector<FitVec> masses = search.run();

  double min_inte = search.getMinInte();
  for (const FitVec& series : masses) {
    for (const Fit& fit : series) {
      MatchEnvPtr env = getMatchEnv(peak_list, fit, env_para_ptr, min_inte);
      if (env != nullptr) {
        envs.push_back(env);
      }
    }
  }
  for (const Fit& fit : search.getAlts()) {
    MatchEnvPtr env = getMatchEnv(peak_list, fit, env_para_ptr, min_inte);
    if (env != nullptr) {
      envs.push_back(env);
    }
  }
  // as DeconvSingleSp::postprocess
  PeakPtrVec peaks = peak_list;
  match_env_util::assignIntensity(peaks, envs);
  onnx_env_cnn::computeEnvScores(peak_list, envs);
  envs = match_env_filter::filterByEnvCnnScore(
      envs, std::max(getPara().min_envcnn,
                     topfd_para_ptr->getMs2EnvCnnScoreCutoff()));
  if (topfd_para_ptr->isAANumBasedFilter()) {
    envs = match_env_filter::filterByAANum(
        envs, max_mass, topfd_para_ptr->isSortUseMsDeconv(), env_para_ptr);
  }
  if (topfd_para_ptr->isSortUseMsDeconv()) {
    std::sort(envs.begin(), envs.end(), MatchEnv::cmpMsdeconvScoreDec);
  } else {
    std::sort(envs.begin(), envs.end(), MatchEnv::cmpEnvcnnScoreDec);
  }
  match_env_util::assignIntensity(peaks, envs);
  return envs;
}

}  // namespace greedy_ms2_deconv

}  // namespace toppic
