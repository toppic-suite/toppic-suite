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

#include "topfd/deconv/prec_pair_refine.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/env/env_base.hpp"
#include "ms/spec/peak_util.hpp"

namespace toppic {

namespace prec_pair_refine {

namespace {

struct RefinePara {
  // pair sums within +- window Da of a precursor mass are used
  double window = 3.0;
  // the width of a cluster of pair sums, in ppm of the precursor mass
  double tol_ppm = 10.0;
  // pairs a cluster needs (with 200 fragment masses about 0.1 random pair
  // sums fall into a 10 ppm window at 20 kDa)
  int min_pairs = 6;
  // 1: a later precursor with a cluster replaces a first one without
  int promote = 1;
  // high-resolution MS1 (correctIsotopes)
  // isotope offsets tested on each side of the precursor mass
  int iso_range = 1;
  // the window of a pair sum around the precursor mass + j isotopes, ppm
  double iso_tol_ppm = 10.0;
  // fragments with a lower EnvCNN score are not paired (TopPIC's -F default)
  double iso_min_envcnn = 0.4;
  // pairs the best offset needs, and its count against the count at P
  int iso_min_pairs = 3;
  double iso_ratio = 2.0;
  // 1: the MS1 averagine fit must not prefer P by iso_ms1_margin or more
  int iso_ms1_check = 1;
  double iso_ms1_margin = 0.05;
  // MS1 fit: m/z tolerance (ppm) and the isotopes of the averagine used
  // (relative intensity to the apex)
  double iso_ms1_tol_ppm = 10.0;
  double iso_ms1_min_rel = 0.05;
};

RefinePara readPara() {
  RefinePara para;
  const char* str = std::getenv("TOPFD_PREC_REFINE_PARA");
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
    if (name == "window") {
      para.window = value;
    } else if (name == "tol_ppm") {
      para.tol_ppm = value;
    } else if (name == "min_pairs") {
      para.min_pairs = static_cast<int>(value);
    } else if (name == "promote") {
      para.promote = static_cast<int>(value);
    } else if (name == "iso_range") {
      para.iso_range = static_cast<int>(value);
    } else if (name == "iso_tol_ppm") {
      para.iso_tol_ppm = value;
    } else if (name == "iso_min_envcnn") {
      para.iso_min_envcnn = value;
    } else if (name == "iso_min_pairs") {
      para.iso_min_pairs = static_cast<int>(value);
    } else if (name == "iso_ratio") {
      para.iso_ratio = value;
    } else if (name == "iso_ms1_check") {
      para.iso_ms1_check = static_cast<int>(value);
    } else if (name == "iso_ms1_margin") {
      para.iso_ms1_margin = value;
    } else if (name == "iso_ms1_tol_ppm") {
      para.iso_ms1_tol_ppm = value;
    } else if (name == "iso_ms1_min_rel") {
      para.iso_ms1_min_rel = value;
    } else {
      LOG_WARN("Unknown TOPFD_PREC_REFINE_PARA entry " << name);
    }
  }
  return para;
}

const RefinePara& getPara() {
  static const RefinePara para = readPara();
  return para;
}

struct Cluster {
  int count = 0;
  double mass = 0;
};

// The densest cluster of pair sums of the sorted masses within +- window
// of center.
Cluster findCluster(const std::vector<double>& masses, double center,
                    const RefinePara& para) {
  std::vector<double> sums;
  for (size_t i = 0; i < masses.size(); i++) {
    double lo = center - para.window - masses[i];
    double hi = center + para.window - masses[i];
    if (hi < masses[i]) {
      break;
    }
    auto first = std::lower_bound(masses.begin() + i + 1, masses.end(), lo);
    for (auto it = first; it != masses.end() && *it <= hi; ++it) {
      sums.push_back(masses[i] + *it);
    }
  }
  Cluster best;
  if (sums.empty()) {
    return best;
  }
  std::sort(sums.begin(), sums.end());
  double tol = center * para.tol_ppm * 1e-6;
  size_t begin = 0;
  size_t best_begin = 0;
  size_t best_end = 0;
  for (size_t end = 0; end < sums.size(); end++) {
    while (sums[end] - sums[begin] > tol) {
      begin++;
    }
    if (static_cast<int>(end - begin + 1) > best.count) {
      best.count = end - begin + 1;
      best_begin = begin;
      best_end = end;
    }
  }
  size_t mid = best_begin + (best_end - best_begin) / 2;
  best.mass = (best_end - best_begin) % 2 == 0
                  ? sums[mid]
                  : (sums[mid] + sums[mid + 1]) / 2;
  return best;
}

// The number of pairs of the sorted masses whose sum is within tol of
// target.
int countPairs(const std::vector<double>& masses, double target, double tol) {
  int count = 0;
  for (size_t i = 0; i < masses.size(); i++) {
    if (2 * masses[i] > target + tol) {
      break;
    }
    auto lo = std::lower_bound(masses.begin() + i + 1, masses.end(),
                               target - tol - masses[i]);
    auto hi = std::upper_bound(masses.begin() + i + 1, masses.end(),
                               target + tol - masses[i]);
    count += hi - lo;
  }
  return count;
}

// The averagine fit of a precursor (mono mass, charge) in the MS1 spectrum
// (m/z and intensities sorted by m/z) at its mass + j isotopes, j =
// -range..range (index j + range): the cosine between the theoretical
// intensities and the most intense peak within tol_ppm of each theoretical
// m/z (0 without one), summed over the charges z-1..z+1. The observed and
// the theoretical intensities are compared on one grid of isotope
// positions around the reported envelope (range isotopes beyond it on each
// side, so the scores depend on range), and a placement is penalized both
// for theoretical peaks without an observed one and for observed peaks it
// leaves out.
std::vector<double> ms1Fit(const std::vector<double>& mzs,
                           const std::vector<double>& intes, double mono_mass,
                           int charge, int range, const RefinePara& para) {
  std::vector<double> scores(2 * range + 1, 0);
  EnvPtr env = EnvBase::getEnvByMonoMass(mono_mass);
  if (env == nullptr || env->getPeakNum() == 0) {
    return scores;
  }
  int peak_num = env->getPeakNum();
  std::vector<double> offsets;
  std::vector<double> rel;
  double max_inte = 0;
  for (int i = 0; i < peak_num; i++) {
    offsets.push_back(env->getMz(i) - env->getMonoMz());
    rel.push_back(env->getInte(i));
    max_inte = std::max(max_inte, env->getInte(i));
  }
  int first = -1;
  int last = -1;
  for (int i = 0; i < peak_num; i++) {
    rel[i] /= max_inte;
    if (rel[i] >= para.iso_ms1_min_rel) {
      if (first < 0) {
        first = i;
      }
      last = i;
    }
  }
  double iso_mass = mass_constant::getIsotopeMass();
  // grid positions: isotope index first - range .. last + range
  std::vector<int> idxes;
  std::vector<double> grid;
  for (int i = first - range; i <= last + range; i++) {
    idxes.push_back(i);
    if (i < 0) {
      grid.push_back(i * iso_mass);
    } else if (i < peak_num) {
      grid.push_back(offsets[i]);
    } else {
      grid.push_back(offsets[peak_num - 1] + (i - peak_num + 1) * iso_mass);
    }
  }
  for (int z = charge - 1; z <= charge + 1; z++) {
    if (z < 1) {
      continue;
    }
    std::vector<double> obs(grid.size(), 0);
    double obs_norm = 0;
    for (size_t t = 0; t < grid.size(); t++) {
      double mz = peak_util::compMz(mono_mass + grid[t], z);
      double tol = mz * para.iso_ms1_tol_ppm * 1e-6;
      auto lo = std::lower_bound(mzs.begin(), mzs.end(), mz - tol);
      auto hi = std::upper_bound(mzs.begin(), mzs.end(), mz + tol);
      for (auto it = lo; it != hi; ++it) {
        obs[t] = std::max(obs[t], intes[it - mzs.begin()]);
      }
      obs_norm += obs[t] * obs[t];
    }
    if (obs_norm == 0) {
      continue;
    }
    obs_norm = std::sqrt(obs_norm);
    for (int j = -range; j <= range; j++) {
      double dot = 0;
      double theo_norm = 0;
      for (size_t t = 0; t < grid.size(); t++) {
        int k = idxes[t] - j;
        if (k >= first && k <= last) {
          dot += obs[t] * rel[k];
          theo_norm += rel[k] * rel[k];
        }
      }
      scores[j + range] += dot / (obs_norm * std::sqrt(theo_norm));
    }
  }
  return scores;
}

}  // namespace

void refine(const MatchEnvPtrVec& frag_envs, PrecursorPtrVec& precs) {
  if (precs.empty() || frag_envs.size() < 2) {
    return;
  }
  const RefinePara& para = getPara();
  std::vector<double> masses;
  for (const MatchEnvPtr& env : frag_envs) {
    masses.push_back(env->getExpEnvPtr()->getMonoNeutralMass());
  }
  std::sort(masses.begin(), masses.end());
  std::vector<Cluster> clusters;
  for (const PrecursorPtr& prec : precs) {
    clusters.push_back(findCluster(masses, prec->getMonoMass(), para));
  }
  int chosen = -1;
  if (clusters[0].count >= para.min_pairs) {
    chosen = 0;
  } else if (para.promote) {
    for (size_t i = 1; i < precs.size(); i++) {
      if (clusters[i].count >= para.min_pairs &&
          (chosen < 0 || clusters[i].count > clusters[chosen].count)) {
        chosen = i;
      }
    }
  }
  if (chosen < 0) {
    return;
  }
  const PrecursorPtr& old = precs[chosen];
  double mono_mz =
      clusters[chosen].mass / old->getCharge() + mass_constant::getProtonMass();
  PrecursorPtr refined = std::make_shared<Precursor>(
      old->getPrecId(), old->getFeatureId(), mono_mz, old->getCharge(),
      old->getInte());
  precs.erase(precs.begin() + chosen);
  precs.insert(precs.begin(), refined);
  for (size_t i = 0; i < precs.size(); i++) {
    precs[i]->setPrecId(i);
  }
}

void correctIsotopes(const MatchEnvPtrVec& frag_envs,
                     const PeakPtrVec& ms1_peaks, PrecursorPtrVec& precs) {
  if (precs.empty() || frag_envs.size() < 2) {
    return;
  }
  const RefinePara& para = getPara();
  std::vector<double> masses;
  for (const MatchEnvPtr& env : frag_envs) {
    if (env->getEnvcnnScore() >= para.iso_min_envcnn) {
      masses.push_back(env->getExpEnvPtr()->getMonoNeutralMass());
    }
  }
  std::sort(masses.begin(), masses.end());
  PeakPtrVec sorted_peaks = ms1_peaks;
  std::sort(sorted_peaks.begin(), sorted_peaks.end(), Peak::cmpPosInc);
  std::vector<double> mzs;
  std::vector<double> intes;
  for (const PeakPtr& peak : sorted_peaks) {
    mzs.push_back(peak->getPosition());
    intes.push_back(peak->getIntensity());
  }
  double iso_mass = mass_constant::getIsotopeMass();
  int range = para.iso_range;
  for (size_t p = 0; p < precs.size(); p++) {
    const PrecursorPtr& prec = precs[p];
    double mass = prec->getMonoMass();
    double tol = mass * para.iso_tol_ppm * 1e-6;
    std::vector<int> counts;
    int best = 0;
    for (int j = -range; j <= range; j++) {
      counts.push_back(countPairs(masses, mass + j * iso_mass, tol));
      // the most pairs; a tie goes to the smaller offset, then to -j
      int c = counts.back();
      int b = counts.size() == 1 ? -1 : counts[best + range];
      if (c > b || (c == b && std::abs(j) < std::abs(best))) {
        best = j;
      }
    }
    int at_prec = counts[range];
    int at_best = counts[best + range];
    if (best == 0 || at_best < para.iso_min_pairs || at_best <= at_prec ||
        at_best < para.iso_ratio * at_prec) {
      continue;
    }
    if (para.iso_ms1_check && !mzs.empty()) {
      // the grid of the fit spans at least two isotopes on each side
      int fit_range = std::max(range, 2);
      std::vector<double> fit =
          ms1Fit(mzs, intes, mass, prec->getCharge(), fit_range, para);
      if (fit[fit_range] - fit[best + fit_range] >= para.iso_ms1_margin) {
        continue;
      }
    }
    double mono_mz = prec->getMonoMz() + best * iso_mass / prec->getCharge();
    precs[p] = std::make_shared<Precursor>(prec->getPrecId(),
                                           prec->getFeatureId(), mono_mz,
                                           prec->getCharge(), prec->getInte());
  }
}

}  // namespace prec_pair_refine

}  // namespace toppic
