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

#include "toplib/lib_search.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "toplib/lib_builder.hpp"
#include "toplib/lib_numeric.hpp"

namespace toppic {

namespace toplib {

namespace lib_search {

namespace {

const double kIsotopeShift = 1.00235;
// dist_cal keeps a query's best library spectrum if its distance is at
// most 1 (every one without a NaN distance)
const double kDistanceThreshold = 1.0;
// fdr_cal: the cutoffs np.arange(0.1, 1.0, 0.001) and the target FDR
const double kCutoffStart = 0.1;
const double kCutoffStop = 1.0;
const double kCutoffStep = 0.001;
const double kTargetFdrPercent = 1.0;
// the written matches have a similarity (1 - distance) of at least 0.3
const double kSimilarityCutoff = 0.3;

struct QueryPeaks {
  std::vector<double> masses;
  std::vector<double> intensities;
  std::vector<int> charges;
};

// ascending with NaN last, as numpy and pandas sort
bool lessNanLast(double a, double b) {
  if (std::isnan(a)) {
    return false;
  }
  return std::isnan(b) || a < b;
}

// ms_spectrum_query_preprocess: the 50 most intense masses (pandas'
// sort_values, descending), sorted by mass (np.sort, np.argsort), log2
// intensities scaled to unit length. Equal intensities or masses keep their
// file order; numpy's default sort is not stable, so the Python script may
// order them differently.
QueryPeaks queryPeaks(const LibSpectrum& sp) {
  std::vector<size_t> idx(sp.masses.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&sp](size_t a, size_t b) {
    double x = sp.intensities[a];
    double y = sp.intensities[b];
    return !std::isnan(x) && (std::isnan(y) || x > y);
  });
  if (idx.size() > static_cast<size_t>(kLibPeakNum)) {
    idx.resize(kLibPeakNum);
  }
  std::stable_sort(idx.begin(), idx.end(), [&sp](size_t a, size_t b) {
    return lessNanLast(sp.masses[a], sp.masses[b]);
  });
  QueryPeaks q;
  for (size_t i : idx) {
    q.masses.push_back(sp.masses[i]);
    q.intensities.push_back(sp.intensities[i]);
    q.charges.push_back(sp.charges[i]);
  }
  q.intensities = lib_numeric::log2Normalize(q.intensities);
  return q;
}

// precursor_mass_query (ppm: within the tolerance, in ppm of the library
// mass, or one isotope shift apart) and precursor_mass_query_Da
bool precursorMatch(double lib_mass, double query_mass,
                    const LibSearchParams& params) {
  double d = std::fabs(lib_mass - query_mass);
  if (params.tol_type == LibTolType::kDa) {
    return d <= params.precursor_tol;
  }
  double tol = std::max(0.0, params.precursor_tol * lib_mass * 1e-6);
  double th1 = 0.0 + tol;
  double th2 = kIsotopeShift + tol;
  double th3 = kIsotopeShift - tol;
  return d <= th1 || (d <= th2 && d >= th3);
}

// A precursor mass distance beyond which precursorMatch is false, with a
// margin for rounding; infinite for a tolerance of 50% or more.
double windowHalfWidth(double query_mass, const LibSearchParams& params) {
  const double margin = 0.01;
  if (params.tol_type == LibTolType::kDa) {
    return params.precursor_tol + margin;
  }
  double rel = params.precursor_tol * 1e-6;
  if (!(rel < 0.5)) {
    return std::numeric_limits<double>::infinity();
  }
  // a library mass m within the window has m <= q + shift + tol and
  // tol = rel * m <= 2 * rel * (q + shift)
  double tol = std::max(0.0, 2 * rel * (std::fabs(query_mass) + 2));
  return kIsotopeShift + tol + margin;
}

// peak_50_cosine_mass_ppm(library row, query, err = 0, ppm): one minus the
// summed intensity products of the peak pairs of equal charge whose masses
// differ by at most the tolerance (in ppm of the library mass) or by one
// isotope shift; a sum above 1 stops the loop without its last product
double cosineDistance(const double* mz1, const double* int1,
                      const double* ch1, const QueryPeaks& q, double ppm) {
  double cos = 0;
  for (int i = 0; i < kLibPeakNum; i++) {
    double tol = std::max(0.0, ppm * mz1[i] * 1e-6);
    double th1 = 0.0 + tol;
    double th2 = kIsotopeShift + tol;
    double th3 = kIsotopeShift - tol;
    for (size_t j = 0; j < q.masses.size(); j++) {
      if (!(ch1[i] == q.charges[j])) {
        continue;
      }
      double d = std::fabs(q.masses[j] - mz1[i]);
      if (d <= th1 || (d <= th2 && d >= th3)) {
        double p = int1[i] * q.intensities[j];
        cos += p;
        if (cos > 1) {
          return 1.0 - (cos - p);
        }
      }
    }
  }
  return 1.0 - cos;
}

// np.arange(0.1, 1.0, 0.001): the first two values, then start + i * delta
// with delta the difference of the first two (numpy's DOUBLE_fill)
std::vector<double> cutoffs() {
  size_t len = static_cast<size_t>(
      std::ceil((kCutoffStop - kCutoffStart) / kCutoffStep));
  std::vector<double> v(len);
  v[0] = kCutoffStart;
  v[1] = kCutoffStart + kCutoffStep;
  double delta = v[1] - v[0];
  for (size_t i = 2; i < len; i++) {
    v[i] = kCutoffStart + static_cast<double>(i) * delta;
  }
  return v;
}

// fdr_cal: for every cutoff with a match at or below it, the decoy fraction
// of those matches; the first cutoff whose FDR is closest to 1%. Then the
// matches written: the target matches at the cutoff with a similarity of
// at least 0.3.
void selectCutoff(const LibSearchLibrary& lib, LibSearchResult& result) {
  std::vector<double> cuts = cutoffs();
  double best = std::numeric_limits<double>::infinity();
  for (size_t c = 0; c < cuts.size(); c++) {
    int target = 0;
    int decoy = 0;
    for (const LibMatch& m : result.matches) {
      if (m.distance <= cuts[c]) {
        double flag = lib.flags[m.lib_row];
        target += flag == 1 ? 1 : 0;
        decoy += flag == 0 ? 1 : 0;
      }
    }
    if (target == 0 && decoy == 0) {
      continue;
    }
    double fdr =
        static_cast<double>(decoy) / static_cast<double>(target + decoy);
    double gap = std::fabs(fdr * 100 - kTargetFdrPercent);
    if (!result.cutoff_found || gap < best) {
      best = gap;
      result.cutoff_found = true;
      result.cutoff_index = static_cast<int>(c);
      result.cutoff = cuts[c];
      result.cutoff_target_num = target;
      result.cutoff_decoy_num = decoy;
      result.cutoff_fdr = fdr;
    }
  }
  if (!result.cutoff_found) {
    return;
  }
  // The script selects the spectrum ids of the target matches at the
  // cutoff and writes every match of those ids; with one match per
  // spectrum id that is each target match at the cutoff, which is what is
  // kept here, so that spectra of a repeated id (a combined msalign file,
  // on which the script fails) do not take each other's rows.
  for (size_t i = 0; i < result.matches.size(); i++) {
    const LibMatch& m = result.matches[i];
    if (m.distance <= result.cutoff && lib.flags[m.lib_row] == 1 &&
        1 - m.distance >= kSimilarityCutoff) {
      result.reported.push_back(static_cast<int>(i));
    }
  }
}

}  // namespace

LibSearchResult search(const LibSearchLibrary& lib, const LibData& queries,
                       const LibSearchParams& params) {
  // the library rows by precursor mass (a NaN mass never matches)
  std::vector<int> by_mass;
  for (int r = 0; r < lib.row_num; r++) {
    if (!std::isnan(lib.precursor_masses[r])) {
      by_mass.push_back(r);
    }
  }
  std::stable_sort(by_mass.begin(), by_mass.end(), [&lib](int a, int b) {
    return lib.precursor_masses[a] < lib.precursor_masses[b];
  });
  std::vector<double> sorted_masses;
  for (int r : by_mass) {
    sorted_masses.push_back(lib.precursor_masses[r]);
  }

  LibSearchResult result;
  int query_num = static_cast<int>(queries.spectra.size());
  for (int s = 0; s < query_num; s++) {
    const LibSpectrum& sp = queries.spectra[s];
    // ms_query_preprocessing: a precursor intensity and two masses
    if (sp.precursor_intensity == 0 || sp.masses.size() <= 1) {
      continue;
    }
    result.searched_num++;
    // the candidates in library order, as the script lists them
    double w = windowHalfWidth(sp.precursor_mass, params);
    auto lo = sorted_masses.begin();
    auto hi = sorted_masses.end();
    if (std::isfinite(w)) {
      lo = std::lower_bound(sorted_masses.begin(), sorted_masses.end(),
                            sp.precursor_mass - w);
      hi = std::upper_bound(lo, sorted_masses.end(), sp.precursor_mass + w);
    }
    std::vector<int> cands;
    for (auto it = lo; it != hi; ++it) {
      int r = by_mass[it - sorted_masses.begin()];
      if (params.charge_match &&
          !(sp.precursor_charge == lib.precursor_charges[r])) {
        continue;
      }
      if (precursorMatch(lib.precursor_masses[r], sp.precursor_mass, params)) {
        cands.push_back(r);
      }
    }
    if (cands.empty()) {
      continue;
    }
    std::sort(cands.begin(), cands.end());

    QueryPeaks q = queryPeaks(sp);
    std::vector<double> dist(cands.size());
    for (size_t k = 0; k < cands.size(); k++) {
      size_t off = static_cast<size_t>(cands[k]) * kLibPeakNum;
      dist[k] = cosineDistance(&lib.masses[off], &lib.norm_intensities[off],
                               &lib.charges[off], q, params.fragment_ppm);
    }
    // np.argmin: the first NaN, else the first minimum; Python's min(): the
    // first value, replaced by every smaller one
    size_t arg = dist.size();
    for (size_t k = 0; k < dist.size(); k++) {
      if (std::isnan(dist[k])) {
        arg = k;
        break;
      }
    }
    if (arg == dist.size()) {
      arg = 0;
      for (size_t k = 1; k < dist.size(); k++) {
        if (dist[k] < dist[arg]) {
          arg = k;
        }
      }
    }
    double py_min = dist[0];
    for (size_t k = 1; k < dist.size(); k++) {
      if (dist[k] < py_min) {
        py_min = dist[k];
      }
    }
    if (py_min <= kDistanceThreshold) {
      result.matches.push_back(LibMatch{s, cands[arg], dist[arg]});
    }
  }

  std::stable_sort(result.matches.begin(), result.matches.end(),
                   [&queries](const LibMatch& a, const LibMatch& b) {
                     return queries.spectra[a.query].spectrum_id <
                            queries.spectra[b.query].spectrum_id;
                   });
  selectCutoff(lib, result);
  return result;
}

}  // namespace lib_search

}  // namespace toplib

}  // namespace toppic
