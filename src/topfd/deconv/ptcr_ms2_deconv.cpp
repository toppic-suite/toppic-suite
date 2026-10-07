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

}  // namespace

PtcrSpecPtr deconv(const MzmlMsPtr& ms_ptr,
                   const TopfdParaPtr& topfd_para_ptr) {
  const PtcrPara& para = getPara();
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  double lo = header_ptr->getPrecWinBegin() - para.win_margin;
  double hi = header_ptr->getPrecWinEnd() + para.win_margin;
  bool has_window = header_ptr->getPrecWinBegin() > 0 && hi > lo;
  IsoRange iso_range(para.iso_min_rel);
  greedy_ms1_deconv::ChargeBound bound;
  if (para.band != 0 && has_window) {
    // charge z of mass M only when some z0 >= z puts the envelope of M in the
    // isolation window
    bound = [lo, hi, iso_min_rel = para.iso_min_rel](double mono_mass,
                                                      int z) {
      // one cache per thread: the search runs in the MS/MS worker threads
      thread_local IsoRange cache(iso_min_rel);
      std::pair<double, double> iso = cache.get(mono_mass);
      for (int z0 = z; z0 <= z + 1000; z0++) {
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
  double sn_ratio =
      para.sn_ratio > 0 ? para.sn_ratio : topfd_para_ptr->getMsTwoSnRatio();
  std::string debug_name;
  const char* debug_scan = std::getenv("TOPFD_PTCR_DEBUG_SCAN");
  int scan = header_ptr->getFirstScanNum();
  if (debug_scan != nullptr && std::atoi(debug_scan) == scan) {
    debug_name = "ptcr_scan" + std::to_string(scan) + ".txt";
  }
  std::vector<MatchEnvPtrVec> masses = greedy_ms1_deconv::deconvSeries(
      ms_ptr->getPeakPtrVec(), topfd_para_ptr, sn_ratio, bound,
      para.min_envcnn, debug_name);

  PtcrSpecPtr spec_ptr = std::make_shared<PtcrSpec>();
  spec_ptr->header_ptr = header_ptr;
  for (const MatchEnvPtrVec& envs : masses) {
    Species species = getSpecies(envs, para.refine);
    if (has_window) {
      std::pair<int, int> iso_charges = getIsoCharges(
          species.mono_mass, iso_range.get(species.mono_mass), lo, hi);
      species.iso_charge_min = iso_charges.first;
      species.iso_charge_max = iso_charges.second;
    }
    spec_ptr->species.push_back(std::move(species));
  }
  std::stable_sort(
      spec_ptr->species.begin(), spec_ptr->species.end(),
      [](const Species& a, const Species& b) { return a.inte > b.inte; });
  PrecursorPtr prec_ptr = getPrecursor(spec_ptr->species, header_ptr,
                                       iso_range);
  header_ptr->setPrecPtrVec(PrecursorPtrVec{prec_ptr});
  return spec_ptr;
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
  std::ofstream tsv(tsv_name);
  tsv << "Spectrum ID\tScan\tSpecies\tPrecursor\tMonoisotopic mass"
         "\tSearch mass\tIntensity\tCharges\tHighest charge"
         "\tIsolation charges\tCharge intensities\tCharge EnvCNN scores"
         "\tBest EnvCNN score\tMatched peaks\tMass MAD (ppm)\n";
  for (const PtcrSpecPtr& spec_ptr : spec_ptrs) {
    const MsHeaderPtr& header_ptr = spec_ptr->header_ptr;
    int sp_id = header_ptr->getSpecId();
    DeconvPeakPtrVec peaks;
    for (std::size_t i = 0; i < spec_ptr->species.size(); i++) {
      const Species& s = spec_ptr->species[i];
      peaks.push_back(std::make_shared<DeconvPeak>(
          sp_id, static_cast<int>(i), s.mono_mass, s.inte, s.charges.back(),
          s.best_envcnn));
    }
    DeconvMsPtr ms_ptr = std::make_shared<DeconvMs>(header_ptr, peaks);
    if (avg_mass) {
      ms_ptr = avg_mass_util::getAvgMassMs(ms_ptr);
    }
    writer.writeMs(ms_ptr);
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
          << "\t" << std::setprecision(2) << s.mass_mad_ppm << "\n";
    }
  }
}

}  // namespace ptcr_ms2_deconv

}  // namespace toppic
