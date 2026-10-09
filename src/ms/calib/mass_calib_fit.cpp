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

#include "ms/calib/mass_calib_fit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/base/activation.hpp"
#include "common/base/activation_base.hpp"
#include "common/base/amino_acid.hpp"
#include "common/base/amino_acid_base.hpp"
#include "common/base/ion_type.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/ptm.hpp"
#include "common/base/ptm_base.hpp"
#include "common/util/logger.hpp"

namespace toppic {

namespace mass_calib_fit {

namespace {

// the isotope spacing used to round a precursor error to whole isotopes
const double kIsotopeSpacing = 1.00235;
const char* const kAcetylPrefix = "[Acetyl]-";

std::vector<std::string> splitTab(const std::string& line) {
  std::vector<std::string> out;
  size_t start = 0;
  while (true) {
    size_t pos = line.find('\t', start);
    if (pos == std::string::npos) {
      out.push_back(line.substr(start));
      return out;
    }
    out.push_back(line.substr(start, pos - start));
    start = pos + 1;
  }
}

bool toDouble(const std::string& s, double* v) {
  try {
    size_t used = 0;
    *v = std::stod(s, &used);
    return used > 0;
  } catch (const std::exception&) {
    return false;
  }
}

// one PrSM row of the table
struct PrsmRow {
  int scan;
  double rt_min;
  int charge;
  double prec_mass;
  double form_mass;
  double evalue;
  int shift_num;
  int var_ptm_num;
  std::string activation;
  std::string proteoform;
};

std::vector<PrsmRow> readTable(const std::string& file_name, int* row_num) {
  std::ifstream in(file_name);
  if (!in) {
    LOG_ERROR("Cannot open the PrSM table " << file_name << "!");
    std::exit(EXIT_FAILURE);
  }
  std::string line;
  bool found = false;
  while (std::getline(in, line)) {
    if (line.rfind("Data file name", 0) == 0) {
      found = true;
      break;
    }
  }
  if (!found) {
    LOG_ERROR("The PrSM table " << file_name
                                << " has no header line \"Data file name\"!");
    std::exit(EXIT_FAILURE);
  }
  if (!line.empty() && line.back() == '\r') line.pop_back();
  std::vector<std::string> header = splitTab(line);
  std::unordered_map<std::string, int> col;
  for (size_t i = 0; i < header.size(); i++) {
    col[header[i]] = static_cast<int>(i);
  }
  const char* required[] = {"Scan(s)",       "Retention time",
                            "Charge",        "Precursor mass",
                            "Proteoform mass", "Proteoform",
                            "Fragmentation", "# Unexpected modifications",
                            "E-value"};
  for (const char* name : required) {
    if (col.find(name) == col.end()) {
      LOG_ERROR("The PrSM table " << file_name << " has no column \"" << name
                                  << "\"!");
      std::exit(EXIT_FAILURE);
    }
  }
  bool has_var = col.find("# Variable PTMs") != col.end();
  std::vector<PrsmRow> rows;
  *row_num = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    std::vector<std::string> f = splitTab(line);
    if (f.size() < header.size()) continue;
    (*row_num)++;
    PrsmRow r;
    double v;
    std::string scan_str = f[col["Scan(s)"]];
    scan_str = scan_str.substr(0, scan_str.find(' '));
    if (!toDouble(scan_str, &v)) continue;
    r.scan = static_cast<int>(v);
    if (!toDouble(f[col["Retention time"]], &r.rt_min)) continue;
    if (!toDouble(f[col["Charge"]], &v)) continue;
    r.charge = static_cast<int>(v);
    if (!toDouble(f[col["Precursor mass"]], &r.prec_mass)) continue;
    if (!toDouble(f[col["Proteoform mass"]], &r.form_mass)) continue;
    if (!toDouble(f[col["E-value"]], &r.evalue)) continue;
    if (!toDouble(f[col["# Unexpected modifications"]], &v)) continue;
    r.shift_num = static_cast<int>(v);
    r.var_ptm_num = 0;
    if (has_var && toDouble(f[col["# Variable PTMs"]], &v)) {
      r.var_ptm_num = static_cast<int>(v);
    }
    r.activation = f[col["Fragmentation"]];
    r.proteoform = f[col["Proteoform"]];
    rows.push_back(r);
  }
  return rows;
}

// {scan: [(mass, charge)]} of an MS/MS msalign file
std::map<int, std::vector<std::pair<double, int>>> readMsAlignMasses(
    const std::string& file_name) {
  std::ifstream in(file_name);
  if (!in) {
    LOG_ERROR("Cannot open the msalign file " << file_name << "!");
    std::exit(EXIT_FAILURE);
  }
  std::map<int, std::vector<std::pair<double, int>>> out;
  std::vector<std::pair<double, int>> cur;
  int scan = -1;
  bool in_spec = false;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("BEGIN IONS", 0) == 0) {
      in_spec = true;
      cur.clear();
      scan = -1;
    } else if (line.rfind("SCANS=", 0) == 0) {
      scan = std::atoi(line.c_str() + 6);
    } else if (line.rfind("END IONS", 0) == 0) {
      if (in_spec && scan >= 0) out[scan] = cur;
      in_spec = false;
    } else if (in_spec && !line.empty() &&
               std::isdigit(static_cast<unsigned char>(line[0]))) {
      std::vector<std::string> f = splitTab(line);
      double mass;
      if (f.size() >= 3 && toDouble(f[0], &mass)) {
        cur.emplace_back(mass, std::atoi(f[2].c_str()));
      }
    }
  }
  return out;
}

// the residue masses (fixed modifications included) and the N-terminal
// acetylation of a proteoform string without shifts or variable PTMs; false
// when the string holds anything else
bool parseProteoform(const std::string& form, std::vector<double>* res,
                     bool* acetyl) {
  res->clear();
  std::string s = form;
  *acetyl = s.rfind(kAcetylPrefix, 0) == 0;
  if (*acetyl) s = s.substr(std::string(kAcetylPrefix).size());
  size_t i = 0;
  while (i < s.size()) {
    char c = s[i];
    if (c >= 'A' && c <= 'Z') {
      if (!AminoAcidBase::containsOneLetter(std::string(1, c))) return false;
      res->push_back(AminoAcidBase::getAminoAcidPtrByOneLetter(c)->getMonoMass());
      i++;
    } else if (c == '(') {
      size_t j = s.find(')', i);
      if (j == std::string::npos || j + 1 >= s.size() || s[j + 1] != '[') {
        return false;
      }
      size_t k = s.find(']', j);
      if (k == std::string::npos) return false;
      std::string group = s.substr(i + 1, j - i - 1);
      std::string mod = s.substr(j + 2, k - j - 2);
      if (group.size() != 1 || mod.find(';') != std::string::npos ||
          !AminoAcidBase::containsOneLetter(group) ||
          !PtmBase::containsAbbrName(mod)) {
        return false;
      }
      res->push_back(
          AminoAcidBase::getAminoAcidPtrByOneLetter(group[0])->getMonoMass() +
          PtmBase::getPtmPtrByAbbrName(mod)->getMonoMass());
      i = k + 1;
    } else {
      return false;
    }
  }
  return !res->empty();
}

}  // namespace

double median(std::vector<double> v) {
  if (v.empty()) return 0.0;
  size_t n = v.size();
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  double hi = v[n / 2];
  if (n % 2 == 1) return hi;
  double lo = *std::max_element(v.begin(), v.begin() + n / 2);
  return (lo + hi) / 2.0;
}

bool addPrecCalibrant(int scan, double rt_min, int charge, double prec_mass,
                      double form_mass, MassCalibrantSet* set) {
  double d = prec_mass - form_mass;
  if (std::nearbyint(d / kIsotopeSpacing) != 0.0 || charge <= 0) {
    return false;
  }
  MassCalibrant p;
  p.scan = scan;
  p.mass = prec_mass;
  p.charge = charge;
  p.mz = (prec_mass + charge * mass_constant::getProtonMass()) / charge;
  p.rt_min = rt_min;
  p.ppm = d / form_mass * 1e6;
  set->prec.push_back(p);
  return true;
}

int addFragCalibrants(int scan, double rt_min,
                      const std::vector<double>& theo_masses,
                      const std::vector<std::pair<double, int>>& masses,
                      double tol_ppm, MassCalibrantSet* set) {
  double proton = mass_constant::getProtonMass();
  int num = 0;
  for (const auto& mc : masses) {
    double m = mc.first;
    size_t pos = std::lower_bound(theo_masses.begin(), theo_masses.end(), m) -
                 theo_masses.begin();
    // the nearest of the two neighbours, the lighter one on a tie
    bool found = false;
    double best = 0.0;
    size_t cand[2] = {pos - 1, pos};
    for (int t = 0; t < 2; t++) {
      if (t == 0 && pos == 0) continue;
      size_t j = cand[t];
      if (j >= theo_masses.size()) continue;
      double e = (m - theo_masses[j]) / theo_masses[j] * 1e6;
      if (std::fabs(e) <= tol_ppm &&
          (!found || std::fabs(e) < std::fabs(best))) {
        best = e;
        found = true;
      }
    }
    if (found && mc.second > 0) {
      MassCalibrant f;
      f.scan = scan;
      f.mass = m;
      f.charge = mc.second;
      f.mz = (m + mc.second * proton) / mc.second;
      f.rt_min = rt_min;
      f.ppm = best;
      set->frag.push_back(f);
      num++;
    }
  }
  return num;
}

MassCalibrantSet readCalibrants(const std::string& prsm_tsv_file,
                                const std::string& ms2_msalign_file,
                                double max_evalue, double tol_ppm) {
  MassCalibrantSet set;
  std::vector<PrsmRow> rows = readTable(prsm_tsv_file, &set.prsm_num);
  auto masses = readMsAlignMasses(ms2_msalign_file);
  double acetyl_mass = PtmBase::getPtmPtr_Acetylation()->getMonoMass();
  std::vector<double> res;
  for (const PrsmRow& r : rows) {
    if (!(r.evalue <= max_evalue) || r.shift_num != 0 || r.var_ptm_num != 0) {
      continue;
    }
    set.used_prsm_num++;
    // the precursor: only the PrSMs at the monoisotopic peak
    addPrecCalibrant(r.scan, r.rt_min, r.charge, r.prec_mass, r.form_mass,
                     &set);
    // the fragments
    bool acetyl = false;
    if (!parseProteoform(r.proteoform, &res, &acetyl)) continue;
    auto it = masses.find(r.scan);
    if (it == masses.end() || it->second.empty()) continue;
    ActivationPtr act = ActivationBase::getActivationPtrByName(r.activation);
    if (act == nullptr) continue;
    double n_shift = act->getNIonTypePtr()->getShift();
    double c_shift = act->getCIonTypePtr()->getShift();
    std::vector<double> theo;
    double sum = 0.0;
    for (size_t i = 0; i + 1 < res.size(); i++) {
      sum += res[i];
      theo.push_back(sum + (acetyl ? acetyl_mass : 0.0) + n_shift);
    }
    sum = 0.0;
    for (size_t i = res.size() - 1; i >= 1; i--) {
      sum += res[i];
      theo.push_back(sum + c_shift);
    }
    std::sort(theo.begin(), theo.end());
    addFragCalibrants(r.scan, r.rt_min, theo, it->second, tol_ppm, &set);
  }
  return set;
}

std::vector<double> robustFit(const std::vector<double>& x, int cols,
                              const std::vector<double>& y,
                              double* robust_sd) {
  size_t n = y.size();
  std::vector<double> w(n, 1.0);
  std::vector<double> c(cols, 0.0);
  std::vector<double> abs_r(n);
  auto residuals = [&]() {
    for (size_t i = 0; i < n; i++) {
      double fit = 0.0;
      for (int k = 0; k < cols; k++) fit += x[i * cols + k] * c[k];
      abs_r[i] = std::fabs(y[i] - fit);
    }
  };
  for (int iter = 0; iter < 10; iter++) {
    // weighted normal equations: (A^T A) c = A^T b with A = W X, b = W y
    std::vector<double> a(cols * cols, 0.0);
    std::vector<double> b(cols, 0.0);
    for (size_t i = 0; i < n; i++) {
      double w2 = w[i] * w[i];
      for (int p = 0; p < cols; p++) {
        b[p] += w2 * x[i * cols + p] * y[i];
        for (int q = 0; q < cols; q++) {
          a[p * cols + q] += w2 * x[i * cols + p] * x[i * cols + q];
        }
      }
    }
    // Gaussian elimination with partial pivoting
    for (int p = 0; p < cols; p++) {
      int piv = p;
      for (int q = p + 1; q < cols; q++) {
        if (std::fabs(a[q * cols + p]) > std::fabs(a[piv * cols + p])) piv = q;
      }
      if (std::fabs(a[piv * cols + p]) < 1e-300) {
        LOG_ERROR("The calibration fit is singular: too few calibrants?");
        std::exit(EXIT_FAILURE);
      }
      if (piv != p) {
        for (int q = 0; q < cols; q++) std::swap(a[p * cols + q], a[piv * cols + q]);
        std::swap(b[p], b[piv]);
      }
      for (int q = p + 1; q < cols; q++) {
        double f = a[q * cols + p] / a[p * cols + p];
        for (int t = p; t < cols; t++) a[q * cols + t] -= f * a[p * cols + t];
        b[q] -= f * b[p];
      }
    }
    for (int p = cols - 1; p >= 0; p--) {
      double s = b[p];
      for (int q = p + 1; q < cols; q++) s -= a[p * cols + q] * c[q];
      c[p] = s / a[p * cols + p];
    }
    residuals();
    double s = 1.4826 * median(abs_r);
    for (size_t i = 0; i < n; i++) {
      w[i] = std::sqrt(std::min(1.0, 1.5 * s / std::max(abs_r[i], 1e-9)));
    }
  }
  residuals();
  *robust_sd = 1.4826 * median(abs_r);
  return c;
}

MassCalibModel fitModel(const MassCalibrantSet& set) {
  if (set.prec.size() < 3 || set.frag.size() < 4) {
    LOG_ERROR("Too few calibrants for the fit: " << set.prec.size()
              << " precursors, " << set.frag.size() << " fragments.");
    std::exit(EXIT_FAILURE);
  }
  MassCalibModel model;
  std::vector<double> rts;
  for (const auto& p : set.prec) rts.push_back(p.rt_min);
  for (const auto& f : set.frag) rts.push_back(f.rt_min);
  model.rt0_min = median(rts);
  // MS1: c0 + c_rt (rt - rt0) / 60
  std::vector<double> x, y;
  for (const auto& p : set.prec) {
    x.push_back(1.0);
    x.push_back((p.rt_min - model.rt0_min) / 60.0);
    y.push_back(p.ppm);
  }
  std::vector<double> c1 = robustFit(x, 2, y, &model.prec_robust_sd);
  model.ms1.c0 = c1[0];
  model.ms1.c_mz = 0.0;
  model.ms1.c_rt_per_hour = c1[1];
  // MS/MS: c0 + c_mz (m/z - mz0) / 1000 + c_rt (rt - rt0) / 60
  x.clear();
  y.clear();
  std::set<int> scans;
  for (const auto& f : set.frag) {
    x.push_back(1.0);
    x.push_back((f.mz - model.mz0) / 1000.0);
    x.push_back((f.rt_min - model.rt0_min) / 60.0);
    y.push_back(f.ppm);
    scans.insert(f.scan);
  }
  std::vector<double> c2 = robustFit(x, 3, y, &model.frag_robust_sd);
  model.ms2.c0 = c2[0];
  model.ms2.c_mz = c2[1];
  model.ms2.c_rt_per_hour = c2[2];
  model.prec_num = static_cast<int>(set.prec.size());
  model.frag_num = static_cast<int>(set.frag.size());
  model.frag_prsm_num = static_cast<int>(scans.size());
  return model;
}

}  // namespace mass_calib_fit

}  // namespace toppic
