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

#include "stat/rescore/prsm_rescore.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>

#include "common/base/prot_mod_base.hpp"
#include "common/util/file_util.hpp"
#include "prsm/expected_value.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "seq/alter_type.hpp"
#include "stat/rescore/rescore_model.hpp"

namespace toppic {

namespace prsm_rescore {

namespace {

// common modification masses (Da) for the shiftKnown feature: deamidation
// and +-1 isotope, -2 H, methyl, oxidation, Na, dimethyl, 2 O, acetyl,
// trimethyl, phospho, carbamidomethyl, -NH3, -H2O, H2O, HexNAc, Hex,
// glutathione, K, Ca, the Gly / Glu residues lost and Hex+HexNAc
const double kModShifts[] = {
    0.984016,  -0.984016, 1.003355,   -1.003355,  2.01565,    -2.01565,
    14.01565,  15.994915, 21.981943,  28.0313,    31.989829,  42.010565,
    42.04695,  79.966331, 57.021464,  -17.026549, -18.010565, 18.010565,
    203.079373, 162.052824, 305.068156, 37.955882, 47.984744, -57.021464,
    -129.042593, 365.132196};
const double kModShiftTolerance = 0.02;

const char kAminoAcids[] = "ACDEFGHIKLMNPQRSTVWY";

double log10Floor(double v, double floor) {
  return std::log10(std::max(v, floor));
}

double median(std::vector<double> v) {
  if (v.empty()) {
    return 0.0;
  }
  std::sort(v.begin(), v.end());
  size_t m = v.size() / 2;
  return v.size() % 2 == 1 ? v[m] : (v[m - 1] + v[m]) / 2.0;
}

bool isDecoy(const PrsmPtr& prsm_ptr) {
  return prsm_ptr->getProteoformPtr()->getSeqName().find("DECOY_") == 0;
}

// the amino acid composition, log10 length and log10 mass of the
// proteoform (the retention time model's input)
std::vector<double> compComposition(const PrsmPtr& prsm_ptr) {
  ProteoformPtr proteoform_ptr = prsm_ptr->getProteoformPtr();
  std::string seq = proteoform_ptr->getResSeqPtr()->toAcidString();
  std::vector<double> comp;
  double len = std::max<size_t>(seq.size(), 1);
  for (int k = 0; k < 20; k++) {
    comp.push_back(std::count(seq.begin(), seq.end(), kAminoAcids[k]) / len);
  }
  comp.push_back(std::log10(len));
  comp.push_back(log10Floor(prsm_ptr->getAdjustedPrecMass(), 100.0));
  return comp;
}

std::pair<int, int> prsmKey(const PrsmPtr& prsm_ptr) {
  return std::make_pair(prsm_ptr->getSpectrumId(),
                        prsm_ptr->getPrecursorId());
}

}  // namespace

const std::vector<std::string>& getFeatureNames() {
  static const std::vector<std::string> names = {
      "lgE",        "lgTest",    "lgOneProt",  "frags",      "peaks",
      "peakFrac",   "fragPerRes", "lgMasses",  "shift0",     "shift2",
      "nVar",       "lgAbsShift", "shiftKnown", "nme",       "acet",
      "nIntact",    "cIntact",   "lgLen",      "lgMass",     "charge",
      "lgPrecInte", "precFrac",  "nPrec",      "lgFeatInte", "featScore",
      "intFrac",    "ppmMed",    "absPpmMed",  "ppmMad",     "isoFrac",
      "siteCov",    "longestRun", "compPairs"};
  return names;
}

std::vector<double> compFeatures(const PrsmPtr& prsm_ptr,
                                 const PrsmPtr& filtered_prsm_ptr,
                                 double min_mass) {
  ProteoformPtr proteoform_ptr = prsm_ptr->getProteoformPtr();
  ExpectedValuePtr ev_ptr = prsm_ptr->getExpectedValuePtr();
  double e_value = ev_ptr == nullptr ? 1e300 : ev_ptr->getEValue();
  double test_num = ev_ptr == nullptr ? 1e300 : ev_ptr->getTestNum();
  double one_prot_prob = ev_ptr == nullptr ? 1.0 : ev_ptr->getOneProtProb();
  DeconvMsPtrVec ms_vec = prsm_ptr->getDeconvMsPtrVec();
  double mass_num = 0;
  for (const auto& ms_ptr : ms_vec) {
    mass_num += ms_ptr->size();
  }
  mass_num = std::max(mass_num, 1.0);
  int len = std::max(proteoform_ptr->getLen(), 1);
  int prot_len = static_cast<int>(
      proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec().size());
  MassShiftPtrVec shifts =
      proteoform_ptr->getMassShiftPtrVec(AlterType::UNEXPECTED);
  double abs_shift = 0;
  bool shift_known = false;
  for (const auto& shift_ptr : shifts) {
    double m = shift_ptr->getMassShift();
    abs_shift += std::abs(m);
    for (double s : kModShifts) {
      if (std::abs(m - s) <= kModShiftTolerance) {
        shift_known = true;
      }
    }
  }
  std::string mod_type = proteoform_ptr->getProtModPtr()->getType();
  MsHeaderPtr header_ptr = ms_vec[0]->getMsHeaderPtr();
  const PrecursorPtrVec& precs = header_ptr->getPrecPtrVec();
  double max_inte = 0;
  for (const auto& prec_ptr : precs) {
    max_inte = std::max(max_inte, prec_ptr->getInte());
  }
  // the PrSM's precursor: the spectrum's first, or the one of its
  // precursor id in the scan's list
  PrecursorPtr prec_ptr =
      precs.empty() ? nullptr
                    : header_ptr->getPrecPtrById(prsm_ptr->getPrecursorId());
  double prec_inte = prec_ptr == nullptr ? 0.0 : prec_ptr->getInte();
  int prec_charge = prec_ptr == nullptr ? 0 : prec_ptr->getCharge();

  // the fragment matches on the masses the search used
  double int_frac = 0;
  double ppm_med = 0;
  double ppm_mad = 10.0;
  double iso_frac = 0;
  double site_cov = 0;
  double longest_run = 0;
  double comp_pairs = 0;
  if (filtered_prsm_ptr != nullptr) {
    PeakIonPairPtrVec pairs = peak_ion_pair_util::genePeakIonPairs(
        proteoform_ptr, filtered_prsm_ptr->getRefineMsPtrVec(), min_mass);
    double total_inte = 0;
    for (const auto& ms_ptr : filtered_prsm_ptr->getDeconvMsPtrVec()) {
      for (const auto& peak_ptr : ms_ptr->getPeakPtrVec()) {
        total_inte += peak_ptr->getIntensity();
      }
    }
    std::set<std::pair<int, int>> matched_peaks;
    double matched_inte = 0;
    std::vector<double> errs;
    int iso_num = 0;
    std::set<int> n_sites;
    std::set<int> c_sites;
    for (const auto& pair_ptr : pairs) {
      const ExtendPeakPtr& real_ptr = pair_ptr->getRealPeakPtr();
      DeconvPeakPtr base_ptr = real_ptr->getBasePeakPtr();
      auto key = std::make_pair(pair_ptr->getSpecId(), base_ptr->getPeakId());
      if (matched_peaks.insert(key).second) {
        matched_inte += base_ptr->getIntensity();
      }
      if (std::abs(real_ptr->getMonoMass() - base_ptr->getMonoMass()) > 0.5) {
        iso_num++;
      } else {
        double theo = pair_ptr->getTheoPeakPtr()->getPosition();
        errs.push_back((real_ptr->getPosition() - theo) / theo * 1e6);
      }
      IonPtr ion_ptr = pair_ptr->getTheoPeakPtr()->getIonPtr();
      if (ion_ptr->getIonTypePtr()->isNTerm()) {
        n_sites.insert(ion_ptr->getPos());
      } else {
        c_sites.insert(ion_ptr->getPos());
      }
    }
    if (total_inte > 0) {
      int_frac = matched_inte / total_inte;
    }
    if (!errs.empty()) {
      ppm_med = median(errs);
    }
    if (errs.size() > 1) {
      std::vector<double> dev;
      for (double e : errs) {
        dev.push_back(std::abs(e - ppm_med));
      }
      ppm_mad = median(dev);
    }
    if (!pairs.empty()) {
      iso_frac = static_cast<double>(iso_num) / pairs.size();
    }
    std::set<int> sites(n_sites);
    sites.insert(c_sites.begin(), c_sites.end());
    double cleavages = std::max(len - 1, 1);
    site_cov = sites.size() / cleavages;
    int run = 0;
    int best = 0;
    int prev = -10;
    for (int s : sites) {
      run = (s == prev + 1) ? run + 1 : 1;
      best = std::max(best, run);
      prev = s;
    }
    longest_run = best / cleavages;
    int both = 0;
    for (int s : n_sites) {
      if (c_sites.count(s) > 0) {
        both++;
      }
    }
    comp_pairs = static_cast<double>(both) / len;
  }

  double lg_e = -log10Floor(e_value, 1e-300);
  std::vector<double> f = {
      std::max(-5.0, std::min(lg_e, 100.0)),
      log10Floor(test_num, 1.0),
      one_prot_prob > 0 ? -log10Floor(one_prot_prob, 1e-300) : 0.0,
      prsm_ptr->getMatchTheoFragNum(),
      prsm_ptr->getMatchExpFragNum(),
      prsm_ptr->getMatchExpFragNum() / mass_num,
      prsm_ptr->getMatchTheoFragNum() / len,
      std::log10(mass_num),
      shifts.empty() ? 1.0 : 0.0,
      shifts.size() == 2 ? 1.0 : 0.0,
      static_cast<double>(proteoform_ptr->getVarPtmNum()),
      std::log10(1.0 + abs_shift),
      shift_known ? 1.0 : 0.0,
      (mod_type == ProtModBase::getType_NME() ||
       mod_type == ProtModBase::getType_NME_ACETYLATION())
          ? 1.0
          : 0.0,
      (mod_type == ProtModBase::getType_NME_ACETYLATION() ||
       mod_type == ProtModBase::getType_M_ACETYLATION())
          ? 1.0
          : 0.0,
      proteoform_ptr->getStartPos() <= 1 ? 1.0 : 0.0,
      proteoform_ptr->getEndPos() >= prot_len - 1 ? 1.0 : 0.0,
      std::log10(static_cast<double>(len)),
      log10Floor(prsm_ptr->getAdjustedPrecMass(), 100.0),
      static_cast<double>(prec_charge),
      std::log10(1.0 + std::max(prec_inte, 0.0)),
      max_inte > 0 ? prec_inte / max_inte : 1.0,
      static_cast<double>(std::max<size_t>(precs.size(), 1)),
      std::log10(1.0 + std::max(prsm_ptr->getFracFeatureInte(), 0.0)),
      std::max(prsm_ptr->getFracFeatureScore(), 0.0),
      int_frac,
      ppm_med,
      std::abs(ppm_med),
      ppm_mad,
      iso_frac,
      site_cov,
      longest_run,
      comp_pairs};
  return f;
}

void process(const PrsmParaPtr& prsm_para_ptr,
             const std::string& input_file_ext,
             const std::string& output_file_ext,
             const std::string& second_sp_file_name,
             const EValueFloor& floor) {
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  // every PrSM with the full spectra (as the recount and the tables read
  // them), and the masses the search used for the fragment features
  PrsmPtrVec prsm_ptrs = prsm_reader_util::readPrsmsWithSpectra(
      prsm_para_ptr, input_file_ext, 0.0, second_sp_file_name);
  PrsmPtrVec filtered_ptrs = prsm_reader_util::readPrsmsWithSpectra(
      prsm_para_ptr, input_file_ext, sp_para_ptr->getEnvCnnCutoff(),
      second_sp_file_name);
  std::map<std::pair<int, int>, PrsmPtr> filtered_map;
  for (const auto& p : filtered_ptrs) {
    filtered_map[prsmKey(p)] = p;
  }

  // the searches apart (prsm_fdr counts them apart)
  std::map<std::string, std::vector<size_t>> by_tool;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    by_tool[prsm_ptrs[i]->getTool()].push_back(i);
  }
  rescore_model::Para para;
  for (const auto& tool : by_tool) {
    size_t n = tool.second.size();
    std::vector<bool> decoy(n);
    // the PrSMs with the maximum E-value (EValueFloor) and the rows the
    // model is trained on and scores
    std::vector<bool> floored(n, false);
    std::vector<size_t> model_rows;
    for (size_t k = 0; k < n; k++) {
      const PrsmPtr& p = prsm_ptrs[tool.second[k]];
      decoy[k] = isDecoy(p);
      floored[k] = floor.enabled &&
                   p->getEValue() >= ExpectedValue::getMaxDouble();
      if (!floored[k] || floor.train) {
        model_rows.push_back(k);
      }
    }
    rescore_model::Data data;
    for (size_t k : model_rows) {
      const PrsmPtr& p = prsm_ptrs[tool.second[k]];
      auto it = filtered_map.find(prsmKey(p));
      PrsmPtr fp = it == filtered_map.end() ? nullptr : it->second;
      data.features.push_back(compFeatures(p, fp, sp_para_ptr->getMinMass()));
      data.decoy.push_back(decoy[k]);
      data.spec_ids.push_back(p->getSpectrumId());
      data.comp.push_back(compComposition(p));
      data.rt.push_back(p->getRetentionTime());
    }
    data.init_col = 0;
    bool ok = false;
    std::vector<double> weights;
    std::vector<double> model_scores =
        rescore_model::rescore(data, para, &ok, &weights);
    std::vector<double> scores(n);
    if (!ok) {
      // the -log10 E-value: the FDR then gives the E-value's result
      for (size_t k = 0; k < n; k++) {
        double e = prsm_ptrs[tool.second[k]]->getEValue();
        scores[k] = -log10Floor(e, 1e-300);
      }
    } else {
      for (size_t j = 0; j < model_rows.size(); j++) {
        scores[model_rows[j]] = model_scores[j];
      }
      size_t n_floored = std::count(floored.begin(), floored.end(), true);
      if (n_floored > 0) {
        double value = floor.value;
        if (std::isnan(value)) {
          double low = std::numeric_limits<double>::infinity();
          for (size_t k = 0; k < n; k++) {
            if (!floored[k]) {
              low = std::min(low, scores[k]);
            }
          }
          value = std::isinf(low) ? 0.0 : low - 1.0;
        }
        size_t n_floored_decoy = 0;
        for (size_t k = 0; k < n; k++) {
          if (floored[k]) {
            scores[k] = value;
            n_floored_decoy += decoy[k] ? 1 : 0;
          }
        }
        std::cout << "Rescoring - " << tool.first << ": " << n_floored
                  << " PrSMs with the maximum E-value (" << n_floored_decoy
                  << " decoys) get the rescore " << std::fixed
                  << std::setprecision(2) << value << std::defaultfloat
                  << (floor.train ? "" : " and are left out of the training")
                  << "." << std::endl;
      }
    }
    for (size_t k = 0; k < n; k++) {
      prsm_ptrs[tool.second[k]]->setRescore(scores[k]);
    }
    // the console summary: targets at 1% FDR by E-value and by the score
    std::vector<double> lg_e;
    for (size_t i : tool.second) {
      lg_e.push_back(-log10Floor(prsm_ptrs[i]->getEValue(), 1e-300));
    }
    double thr;
    std::vector<bool> e_pass =
        rescore_model::passFdr(lg_e, decoy, para.train_fdr, &thr);
    std::vector<bool> r_pass =
        rescore_model::passFdr(scores, decoy, para.train_fdr, &thr);
    size_t n_decoy = std::count(decoy.begin(), decoy.end(), true);
    std::cout << "Rescoring - " << tool.first << ": " << tool.second.size()
              << " PrSMs (" << n_decoy << " decoys); target PrSMs at "
              << para.train_fdr * 100 << "% FDR: "
              << std::count(e_pass.begin(), e_pass.end(), true)
              << " by E-value, "
              << std::count(r_pass.begin(), r_pass.end(), true)
              << (ok ? " rescored." : " (no model: the E-value is used).")
              << std::endl;
    if (ok) {
      std::vector<std::string> names = getFeatureNames();
      names.push_back("rtResid");
      std::vector<size_t> order(weights.size());
      for (size_t j = 0; j < order.size(); j++) {
        order[j] = j;
      }
      std::stable_sort(order.begin(), order.end(),
                       [&weights](size_t a, size_t b) {
                         return std::abs(weights[a]) > std::abs(weights[b]);
                       });
      std::ostringstream os;
      os << std::showpos << std::fixed << std::setprecision(2);
      for (size_t j = 0; j < std::min<size_t>(order.size(), 8); j++) {
        os << (j == 0 ? "" : ", ") << names[order[j]] << " "
           << weights[order[j]];
      }
      std::cout << "Rescoring - " << tool.first
                << ": largest weights (standardized features): " << os.str()
                << std::endl;
    }
  }

  std::string output_file_name =
      file_util::basename(prsm_para_ptr->getSpectrumFileName()) + "." +
      output_file_ext;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}

}  // namespace prsm_rescore

}  // namespace toppic
