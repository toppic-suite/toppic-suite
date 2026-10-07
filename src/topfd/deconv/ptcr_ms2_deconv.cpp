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

#include "topfd/deconv/ptcr_ms2_deconv.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/env/avg_mass_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/match_env.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/precursor.hpp"
#include "topfd/deconv/greedy_ms1_deconv.hpp"

namespace toppic {

namespace ptcr_ms2_deconv {

namespace {

// The parameters, set for research with the environment variable
// TOPFD_PTCR_PARA (e.g. "band=0,win_margin=0.5").
struct PtcrPara {
  // the PTCR bound of the search (1) or none (0)
  int band = 1;
  // m/z added to both ends of the isolation window in the bound
  double win_margin = 0;
  // the isotopes of a mass that count in the bound: relative intensity >=
  // iso_min_rel of the apex
  double iso_min_rel = 0.1;
  // envelopes with a lower EnvCNN score are dropped (as greedy MS1's
  // min_envcnn)
  double min_envcnn = 0.2;
  // the signal-to-noise ratio of the seed threshold; 0: the MS/MS ratio (-s)
  double sn_ratio = 0;
  // the mass refined from the matched isotope peaks (1) or the search's (0)
  int refine = 1;
  // fragment PTCR spectra (getFragMode): 1 the MS/MS search with the
  // selected ion as the precursor when no MS1 feature is; 0 as other MS/MS
  // spectra; 2 the charge-series search (deconvFragments)
  int frag = 1;
  // frag 2: fragment envelopes with a lower EnvCNN score are dropped (0:
  // none; the search's -F filters the masses)
  double frag_min_envcnn = 0;
  // frag 2: fragment species within this mass (Da) of the precursor mass,
  // or above it, are not fragments (the surviving precursor)
  double frag_prec_margin = 20;
  // a species matches an MS1 feature within this tolerance (ppm), also one
  // isotope apart
  double feature_ppm = 10;
};

PtcrPara readPara() {
  PtcrPara para;
  const char* str = std::getenv("TOPFD_PTCR_PARA");
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
    if (name == "band") {
      para.band = static_cast<int>(value);
    } else if (name == "win_margin") {
      para.win_margin = value;
    } else if (name == "iso_min_rel") {
      para.iso_min_rel = value;
    } else if (name == "min_envcnn") {
      para.min_envcnn = value;
    } else if (name == "sn_ratio") {
      para.sn_ratio = value;
    } else if (name == "refine") {
      para.refine = static_cast<int>(value);
    } else if (name == "frag") {
      para.frag = static_cast<int>(value);
    } else if (name == "frag_min_envcnn") {
      para.frag_min_envcnn = value;
    } else if (name == "frag_prec_margin") {
      para.frag_prec_margin = value;
    } else if (name == "feature_ppm") {
      para.feature_ppm = value;
    } else {
      LOG_WARN("Unknown TOPFD_PTCR_PARA entry " << name);
    }
  }
  return para;
}

const PtcrPara& getPara() {
  static const PtcrPara para = readPara();
  return para;
}

// The mass offsets (Da) of the first and the last isotope of a mass's
// averagine distribution with a relative intensity >= iso_min_rel, cached
// per 10 Da bin.
class IsoRange {
 public:
  explicit IsoRange(double iso_min_rel) : iso_min_rel_(iso_min_rel) {}

  std::pair<double, double> get(double mono_mass) {
    int key = static_cast<int>(mono_mass / 10);
    auto it = cache_.find(key);
    if (it != cache_.end()) {
      return it->second;
    }
    std::pair<double, double> range(0, 0);
    EnvPtr env = EnvBase::getEnvByMonoMass(mono_mass);
    if (env != nullptr && env->getPeakNum() > 0) {
      double mono = env->getMonoMz();
      double max_inte = env->getReferInte();
      bool found = false;
      for (int i = 0; i < env->getPeakNum(); i++) {
        if (env->getInte(i) >= iso_min_rel_ * max_inte) {
          double offset = env->getMz(i) - mono;
          if (!found) {
            range.first = offset;
            found = true;
          }
          range.second = offset;
        }
      }
    }
    cache_.emplace(key, range);
    return range;
  }

 private:
  double iso_min_rel_;
  std::unordered_map<int, std::pair<double, double>> cache_;
};

// Whether the envelope of a mass at charge z0 overlaps the m/z interval
// [lo, hi]: -1 when it lies below lo, 1 above hi, 0 overlapping.
int compWindow(double mono_mass, int z0, const std::pair<double, double>& iso,
               double lo, double hi) {
  double proton = mass_constant::getProtonMass();
  double first_mz = (mono_mass + iso.first) / z0 + proton;
  double last_mz = (mono_mass + iso.second) / z0 + proton;
  if (last_mz < lo) {
    return -1;
  }
  if (first_mz > hi) {
    return 1;
  }
  return 0;
}

// The charges z0 at which a mass's envelope overlaps [lo, hi]: [min, max],
// or (0, 0). The m/z falls as z0 grows, so they are consecutive.
std::pair<int, int> getIsoCharges(double mono_mass,
                                  const std::pair<double, double>& iso,
                                  double lo, double hi) {
  std::pair<int, int> result(0, 0);
  for (int z0 = 1; z0 <= 1000; z0++) {
    int pos = compWindow(mono_mass, z0, iso, lo, hi);
    if (pos < 0) {
      break;
    }
    if (pos == 0) {
      if (result.first == 0) {
        result.first = z0;
      }
      result.second = z0;
    }
  }
  return result;
}

// The intensity-weighted median of values
double weightedMedian(std::vector<std::pair<double, double>> values) {
  if (values.empty()) {
    return 0;
  }
  std::sort(values.begin(), values.end());
  double total = 0;
  for (const auto& v : values) {
    total += v.second;
  }
  double sum = 0;
  for (const auto& v : values) {
    sum += v.second;
    if (sum >= total / 2) {
      return v.first;
    }
  }
  return values.back().first;
}

// One species from the envelopes of a mass
Species getSpecies(const MatchEnvPtrVec& envs, int refine) {
  Species species;
  species.search_mass = envs[0]->getTheoEnvPtr()->getMonoNeutralMass();
  // the mass errors of the matched isotope peaks: (observed - theoretical)
  // m/z x charge, weighted by the observed intensity
  std::vector<std::pair<double, double>> errors;
  std::vector<std::pair<int, std::pair<double, double>>> charges;
  for (const MatchEnvPtr& env : envs) {
    EnvPtr theo = env->getTheoEnvPtr();
    ExpEnvPtr exp = env->getExpEnvPtr();
    int z = theo->getCharge();
    for (int i = 0; i < theo->getPeakNum() && i < exp->getPeakNum(); i++) {
      if (exp->isExist(i)) {
        errors.emplace_back((exp->getMz(i) - theo->getMz(i)) * z,
                            exp->getInte(i));
      }
    }
    double inte = theo->compInteSum();
    species.inte += inte;
    species.best_envcnn = std::max(species.best_envcnn, env->getEnvcnnScore());
    charges.emplace_back(z, std::make_pair(inte, env->getEnvcnnScore()));
  }
  std::sort(charges.begin(), charges.end());
  for (const auto& c : charges) {
    species.charges.push_back(c.first);
    species.charge_intes.push_back(c.second.first);
    species.charge_envcnns.push_back(c.second.second);
  }
  double shift = weightedMedian(errors);
  species.mono_mass = species.search_mass + (refine != 0 ? shift : 0);
  species.match_peak_num = static_cast<int>(errors.size());
  std::vector<std::pair<double, double>> devs;
  for (const auto& e : errors) {
    devs.emplace_back(std::abs(e.first - shift), e.second);
  }
  species.mass_mad_ppm = weightedMedian(devs) / species.mono_mass * 1e6;
  species.envs = envs;
  return species;
}

// Whether the envelope of a species at charge z covers the m/z mz: the m/z
// range of its isotopes of relative intensity >= iso_min_rel
bool covers(const Species& species, int z, double mz,
            const std::pair<double, double>& iso) {
  double proton = mass_constant::getProtonMass();
  double first_mz = (species.mono_mass + iso.first) / z + proton;
  double last_mz = (species.mono_mass + iso.second) / z + proton;
  return mz >= first_mz && mz <= last_mz;
}

bool hasCharge(const Species& species, int z) {
  return std::find(species.charges.begin(), species.charges.end(), z) !=
         species.charges.end();
}

// The precursor of the spectrum (see the header); marks its species
PrecursorPtr getPrecursor(std::vector<Species>& species_vec,
                          const MsHeaderPtr& header_ptr, IsoRange& iso_range) {
  double sel_mz = header_ptr->getSelectedMz();
  if (sel_mz <= 0) {
    sel_mz = header_ptr->getPrecTargetMz();
  }
  int sel_z = header_ptr->getSelectedCharge();
  // species_vec is sorted by intensity, most intense first
  int best = -1;
  int best_z = 0;
  if (sel_z > 0) {
    for (std::size_t i = 0; i < species_vec.size() && best < 0; i++) {
      const Species& s = species_vec[i];
      if (hasCharge(s, sel_z) &&
          covers(s, sel_z, sel_mz, iso_range.get(s.mono_mass))) {
        best = static_cast<int>(i);
        best_z = sel_z;
      }
    }
    for (std::size_t i = 0; i < species_vec.size() && best < 0; i++) {
      const Species& s = species_vec[i];
      if (s.iso_charge_min <= sel_z && sel_z <= s.iso_charge_max) {
        best = static_cast<int>(i);
        best_z = sel_z;
      }
    }
  }
  if (best < 0 && !species_vec.empty()) {
    best = 0;
    const Species& s = species_vec[0];
    best_z = s.iso_charge_max > 0 ? s.iso_charge_max : s.charges.back();
  }
  if (best < 0) {
    // no species: the selected ion as the file gives it
    int z = sel_z > 0 ? sel_z : 1;
    return std::make_shared<Precursor>(0, -1, sel_mz, z, 0.0);
  }
  Species& s = species_vec[best];
  s.is_prec = true;
  double mono_mz = s.mono_mass / best_z + mass_constant::getProtonMass();
  return std::make_shared<Precursor>(0, -1, mono_mz, best_z, s.inte);
}

std::string joinInts(const std::vector<int>& values) {
  std::ostringstream os;
  for (std::size_t i = 0; i < values.size(); i++) {
    os << (i > 0 ? "," : "") << values[i];
  }
  return os.str();
}

// The charge bound of the search: charge z of mass M only when some z0 in
// [z, z_cap] puts the envelope of M (its isotopes >= iso_min_rel of the
// apex) in [lo, hi]; masses above max_mass are excluded.
greedy_ms1_deconv::ChargeBound makeBound(double lo, double hi, int z_cap,
                                         double max_mass, double iso_min_rel) {
  return [lo, hi, z_cap, max_mass, iso_min_rel](double mono_mass, int z) {
    if (mono_mass > max_mass || z > z_cap) {
      return false;
    }
    // one cache per thread: the search runs in the MS/MS worker threads
    thread_local IsoRange cache(iso_min_rel);
    std::pair<double, double> iso = cache.get(mono_mass);
    for (int z0 = z; z0 <= z_cap; z0++) {
      int pos = compWindow(mono_mass, z0, iso, lo, hi);
      if (pos < 0) {
        return false;
      }
      if (pos == 0) {
        return true;
      }
    }
    return false;
  };
}

// The isolation window of the proton transfer stage, with the margin
bool getPtcrWindow(const MsHeaderPtr& header_ptr, double margin, double& lo,
                   double& hi) {
  double bgn = header_ptr->getPtcrWinBegin();
  double end = header_ptr->getPtcrWinEnd();
  if (bgn <= 0) {
    bgn = header_ptr->getPrecWinBegin();
    end = header_ptr->getPrecWinEnd();
  }
  lo = bgn - margin;
  hi = end + margin;
  return bgn > 0 && hi > lo;
}

std::string getDebugName(const MsHeaderPtr& header_ptr) {
  const char* debug_scan = std::getenv("TOPFD_PTCR_DEBUG_SCAN");
  int scan = header_ptr->getFirstScanNum();
  if (debug_scan != nullptr && std::atoi(debug_scan) == scan) {
    return "ptcr_scan" + std::to_string(scan) + ".txt";
  }
  return "";
}

// The species of the search's masses, with their isolation charges (z0 <=
// z_cap), most intense first
std::vector<Species> getSpeciesVec(const std::vector<MatchEnvPtrVec>& masses,
                                   bool has_window, double lo, double hi,
                                   IsoRange& iso_range, int z_cap) {
  const PtcrPara& para = getPara();
  std::vector<Species> result;
  for (const MatchEnvPtrVec& envs : masses) {
    Species species = getSpecies(envs, para.refine);
    if (has_window) {
      std::pair<int, int> iso_charges = getIsoCharges(
          species.mono_mass, iso_range.get(species.mono_mass), lo, hi);
      species.iso_charge_min = iso_charges.first;
      species.iso_charge_max = std::min(iso_charges.second, z_cap);
      if (species.iso_charge_max < species.iso_charge_min) {
        species.iso_charge_min = 0;
        species.iso_charge_max = 0;
      }
    }
    result.push_back(std::move(species));
  }
  std::stable_sort(
      result.begin(), result.end(),
      [](const Species& a, const Species& b) { return a.inte > b.inte; });
  return result;
}

// The MS1 feature of a mass: the nearest within feature_ppm, else one within
// feature_ppm one isotope apart; nullptr when none
SpecFeaturePtr matchFeature(double mass, const SpecFeaturePtrVec& feat_ptrs) {
  double tol = getPara().feature_ppm * 1e-6 * mass;
  double iso = mass_constant::getIsotopeMass();
  for (int k : {0, -1, 1}) {
    SpecFeaturePtr best = nullptr;
    double best_err = tol;
    for (const SpecFeaturePtr& feat : feat_ptrs) {
      double err = std::abs(feat->getPrecMonoMass() + k * iso - mass);
      if (err <= best_err) {
        best_err = err;
        best = feat;
      }
    }
    if (best != nullptr) {
      return best;
    }
  }
  return nullptr;
}

// The offset (Da) of the averagine apex from the monoisotopic mass
double getApexOffset(double mono_mass) {
  EnvPtr env = EnvBase::getEnvByMonoMass(mono_mass);
  if (env == nullptr || env->getPeakNum() == 0) {
    return 0;
  }
  return env->getReferMz() - env->getMonoMz();
}

}  // namespace

int getFragMode() { return getPara().frag; }

PrecursorPtr getSelectedIonPrec(const MsHeaderPtr& header_ptr,
                                double& max_mass, int& max_charge) {
  int sel_z = header_ptr->getSelectedCharge();
  double sel_mz = header_ptr->getSelectedMz();
  if (sel_mz <= 0) {
    sel_mz = header_ptr->getPrecTargetMz();
  }
  if (sel_z <= 0 || sel_mz <= 0) {
    return nullptr;
  }
  double half_win =
      (header_ptr->getPrecWinEnd() - header_ptr->getPrecWinBegin()) / 2;
  max_mass = (sel_mz - mass_constant::getProtonMass() + half_win) * sel_z;
  max_charge = sel_z;
  return std::make_shared<Precursor>(0, -1, sel_mz, sel_z, 0.0);
}

PtcrSpecPtr deconv(const MzmlMsPtr& ms_ptr, const TopfdParaPtr& topfd_para_ptr,
                   const SpecFeaturePtrVec& feat_ptrs) {
  const PtcrPara& para = getPara();
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  double lo = 0;
  double hi = 0;
  bool has_window = getPtcrWindow(header_ptr, para.win_margin, lo, hi);
  IsoRange iso_range(para.iso_min_rel);
  // the charges up to -c, or the selected ion's when higher
  int max_charge = std::max(topfd_para_ptr->getMaxCharge(),
                            header_ptr->getSelectedCharge());
  double max_mass = topfd_para_ptr->getMaxMass();
  greedy_ms1_deconv::ChargeBound bound;
  if (para.band != 0 && has_window) {
    // charge z of mass M only when some z0 >= z puts the envelope of M in the
    // isolation window
    bound = makeBound(lo, hi, 100000, max_mass, para.iso_min_rel);
  }
  double sn_ratio =
      para.sn_ratio > 0 ? para.sn_ratio : topfd_para_ptr->getMsTwoSnRatio();
  std::vector<MatchEnvPtrVec> masses = greedy_ms1_deconv::deconvSeries(
      ms_ptr->getPeakPtrVec(), topfd_para_ptr, sn_ratio, max_mass, max_charge,
      bound, para.min_envcnn, getDebugName(header_ptr));

  PtcrSpecPtr spec_ptr = std::make_shared<PtcrSpec>();
  spec_ptr->header_ptr = header_ptr;
  spec_ptr->species =
      getSpeciesVec(masses, has_window, lo, hi, iso_range, 100000);
  for (Species& s : spec_ptr->species) {
    SpecFeaturePtr feat = matchFeature(s.mono_mass, feat_ptrs);
    if (feat != nullptr) {
      s.feature_id = feat->getFracFeatureId();
      s.feature_mass = feat->getPrecMonoMass();
    }
  }
  PrecursorPtr prec_ptr = getPrecursor(spec_ptr->species, header_ptr,
                                       iso_range);
  // the precursor species' feature: its id and intensity
  for (const Species& s : spec_ptr->species) {
    if (s.is_prec && s.feature_id >= 0) {
      for (const SpecFeaturePtr& feat : feat_ptrs) {
        if (feat->getFracFeatureId() == s.feature_id) {
          prec_ptr = std::make_shared<Precursor>(
              0, s.feature_id, prec_ptr->getMonoMz(), prec_ptr->getCharge(),
              feat->getPrecInte());
        }
      }
    }
  }
  header_ptr->setPrecPtrVec(PrecursorPtrVec{prec_ptr});
  return spec_ptr;
}

PtcrSpecPtr deconvFragments(const MzmlMsPtr& ms_ptr,
                            const TopfdParaPtr& topfd_para_ptr,
                            double max_mass, int max_charge,
                            const SpecFeaturePtrVec& feat_ptrs) {
  const PtcrPara& para = getPara();
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  double proton = mass_constant::getProtonMass();
  double lo = 0;
  double hi = 0;
  bool has_window = getPtcrWindow(header_ptr, para.win_margin, lo, hi);
  IsoRange iso_range(para.iso_min_rel);
  // the precursor: its charge caps the fragments' charges, and its mass
  // (here the selected ion's, estimated from the m/z of the window centre)
  // their masses
  int sel_z = header_ptr->getSelectedCharge();
  double sel_mz = header_ptr->getSelectedMz();
  if (sel_mz <= 0) {
    sel_mz = header_ptr->getPrecTargetMz();
  }
  double half_win = (header_ptr->getPrecWinEnd() -
                     header_ptr->getPrecWinBegin()) / 2;
  double sel_mass = sel_z > 0 && sel_mz > 0 ? (sel_mz - proton) * sel_z : 0;
  if (sel_mass > 0) {
    max_charge = feat_ptrs.empty() ? sel_z : std::max(max_charge, sel_z);
    double sel_max = sel_mass + half_win * sel_z + 50;
    max_mass = feat_ptrs.empty() ? sel_max : std::max(max_mass, sel_max);
  }
  greedy_ms1_deconv::ChargeBound bound;
  if (para.band != 0 && has_window) {
    bound = makeBound(lo, hi, max_charge, max_mass, para.iso_min_rel);
  } else {
    bound = [max_mass, max_charge](double mono_mass, int z) {
      return mono_mass <= max_mass && z <= max_charge;
    };
  }
  double sn_ratio =
      para.sn_ratio > 0 ? para.sn_ratio : topfd_para_ptr->getMsTwoSnRatio();
  std::vector<MatchEnvPtrVec> masses = greedy_ms1_deconv::deconvSeries(
      ms_ptr->getPeakPtrVec(), topfd_para_ptr, sn_ratio, max_mass, max_charge,
      bound, para.frag_min_envcnn, getDebugName(header_ptr));
  std::vector<Species> all =
      getSpeciesVec(masses, has_window, lo, hi, iso_range, max_charge);

  // the precursors: the MS1 features; without them, the surviving precursor
  // species (its apex at the selected ion's m/z and charge, within the
  // window), else the selected ion
  PrecursorPtrVec prec_ptrs;
  for (std::size_t i = 0; i < feat_ptrs.size(); i++) {
    const SpecFeaturePtr& feat = feat_ptrs[i];
    prec_ptrs.push_back(std::make_shared<Precursor>(
        static_cast<int>(i), feat->getFracFeatureId(), feat->getPrecMonoMz(),
        feat->getPrecCharge(), feat->getPrecInte()));
  }
  int prec_idx = -1;
  if (prec_ptrs.empty() && sel_mass > 0) {
    for (std::size_t i = 0; i < all.size(); i++) {
      double apex = all[i].mono_mass + getApexOffset(all[i].mono_mass);
      if (std::abs(apex - sel_mass) <= half_win * sel_z + 10 &&
          all[i].charges.size() >= 2) {
        prec_idx = static_cast<int>(i);
        break;
      }
    }
    if (prec_idx >= 0) {
      Species& s = all[prec_idx];
      s.is_prec = true;
      prec_ptrs.push_back(std::make_shared<Precursor>(
          0, -1, s.mono_mass / sel_z + proton, sel_z, s.inte));
    } else {
      prec_ptrs.push_back(
          std::make_shared<Precursor>(0, -1, sel_mz, sel_z, 0.0));
    }
  }
  header_ptr->setPrecPtrVec(prec_ptrs);

  // the fragments: below the precursor mass less the margin
  double prec_mass = max_mass;
  if (!prec_ptrs.empty()) {
    prec_mass = 0;
    for (const PrecursorPtr& prec : prec_ptrs) {
      prec_mass = std::max(prec_mass, prec->getMonoMass());
    }
  }
  PtcrSpecPtr spec_ptr = std::make_shared<PtcrSpec>();
  spec_ptr->header_ptr = header_ptr;
  for (Species& s : all) {
    if (s.is_prec || s.mono_mass >= prec_mass - para.frag_prec_margin) {
      continue;
    }
    spec_ptr->species.push_back(std::move(s));
  }
  return spec_ptr;
}

DeconvMsPtr getDeconvMs(const PtcrSpecPtr& spec_ptr) {
  const MsHeaderPtr& header_ptr = spec_ptr->header_ptr;
  int sp_id = header_ptr->getSpecId();
  DeconvPeakPtrVec peaks;
  for (std::size_t i = 0; i < spec_ptr->species.size(); i++) {
    const Species& s = spec_ptr->species[i];
    peaks.push_back(std::make_shared<DeconvPeak>(
        sp_id, static_cast<int>(i), s.mono_mass, s.inte, s.charges.back(),
        s.best_envcnn));
  }
  return std::make_shared<DeconvMs>(header_ptr, peaks);
}

void writeTsv(PtcrSpecPtrVec spec_ptrs, const std::string& tsv_name) {
  std::sort(spec_ptrs.begin(), spec_ptrs.end(),
            [](const PtcrSpecPtr& a, const PtcrSpecPtr& b) {
              return a->header_ptr->getSpecId() < b->header_ptr->getSpecId();
            });
  std::ofstream tsv(tsv_name);
  tsv << "Spectrum ID\tScan\tSpecies\tPrecursor\tMonoisotopic mass"
         "\tSearch mass\tIntensity\tCharges\tHighest charge"
         "\tIsolation charges\tCharge intensities\tCharge EnvCNN scores"
         "\tBest EnvCNN score\tMatched peaks\tMass MAD (ppm)\tFeature ID"
         "\tFeature mass\n";
  for (const PtcrSpecPtr& spec_ptr : spec_ptrs) {
    const MsHeaderPtr& header_ptr = spec_ptr->header_ptr;
    int sp_id = header_ptr->getSpecId();
    for (std::size_t i = 0; i < spec_ptr->species.size(); i++) {
      const Species& s = spec_ptr->species[i];
      std::ostringstream intes;
      std::ostringstream scores;
      for (std::size_t j = 0; j < s.charges.size(); j++) {
        intes << (j > 0 ? "," : "") << std::fixed << std::setprecision(2)
              << s.charge_intes[j];
        scores << (j > 0 ? "," : "") << std::fixed << std::setprecision(3)
               << s.charge_envcnns[j];
      }
      tsv << sp_id << "\t" << header_ptr->getScansString() << "\t" << i
          << "\t" << (s.is_prec ? 1 : 0) << "\t" << std::fixed
          << std::setprecision(6) << s.mono_mass << "\t" << s.search_mass
          << "\t" << std::setprecision(2) << s.inte << "\t"
          << joinInts(s.charges) << "\t" << s.charges.back() << "\t";
      if (s.iso_charge_min > 0) {
        tsv << s.iso_charge_min;
        if (s.iso_charge_max > s.iso_charge_min) {
          tsv << "-" << s.iso_charge_max;
        }
      }
      tsv << "\t" << intes.str() << "\t" << scores.str() << "\t"
          << std::setprecision(3) << s.best_envcnn << "\t" << s.match_peak_num
          << "\t" << std::setprecision(2) << s.mass_mad_ppm << "\t"
          << s.feature_id << "\t";
      if (s.feature_id >= 0) {
        tsv << std::setprecision(6) << s.feature_mass;
      }
      tsv << "\n";
    }
  }
}

void write(PtcrSpecPtrVec spec_ptrs, const std::string& msalign_name,
           const std::string& tsv_name, const std::string& para_str,
           bool avg_mass) {
  std::sort(spec_ptrs.begin(), spec_ptrs.end(),
            [](const PtcrSpecPtr& a, const PtcrSpecPtr& b) {
              return a->header_ptr->getSpecId() < b->header_ptr->getSpecId();
            });
  MsAlignWriter writer(msalign_name);
  writer.writePara(para_str);
  for (const PtcrSpecPtr& spec_ptr : spec_ptrs) {
    DeconvMsPtr ms_ptr = getDeconvMs(spec_ptr);
    if (avg_mass) {
      ms_ptr = avg_mass_util::getAvgMassMs(ms_ptr);
    }
    writer.writeMs(ms_ptr);
  }
  writeTsv(spec_ptrs, tsv_name);
}

}  // namespace ptcr_ms2_deconv

}  // namespace toppic
