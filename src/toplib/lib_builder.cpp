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


#include "toplib/lib_builder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "toplib/lib_cluster.hpp"
#include "toplib/lib_numeric.hpp"

namespace toppic {

namespace toplib {

namespace lib_builder {

namespace {

// ms_library_building.py's constants
const double kPrecursorTolerance = 2.2;
const double kSpectrumDistance = 0.7;
const double kPpm = 10;
const double kIsotopeShift = 1.00235;

// the 19 amino acid residue masses of decoy_mass_gen (no leucine)
const double kAaMasses[19] = {
    71.03711381,  156.10111105, 114.04292747, 115.02694307, 103.00918451,
    129.04259314, 128.05857754, 57.02146373,  137.05891188, 113.08406401,
    128.09496305, 131.04048464, 147.06841395, 97.05276387,  87.03202844,
    101.0476785,  186.07931298, 163.06332857, 99.06841394};

struct Peaks {
  std::vector<double> masses;
  std::vector<double> intensities;
  std::vector<int> charges;
};

Peaks peaksOf(const LibSpectrum& sp) {
  return Peaks{sp.masses, sp.intensities, sp.charges};
}

Peaks select(const Peaks& p, const std::vector<size_t>& idx) {
  Peaks r;
  for (size_t i : idx) {
    r.masses.push_back(p.masses[i]);
    r.intensities.push_back(p.intensities[i]);
    r.charges.push_back(p.charges[i]);
  }
  return r;
}

// np.argsort. numpy's default sort is not stable; a stable one is used
// here, so the order of equal values may differ from the Python script's.
std::vector<size_t> argsort(const std::vector<double>& v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t a, size_t b) { return v[a] < v[b]; });
  return idx;
}

Peaks sortByMass(const Peaks& p) { return select(p, argsort(p.masses)); }

// ms_incorr_remove: np.argsort(inte)[::-1][0:50]
Peaks topPeaksDescending(const Peaks& p) {
  std::vector<size_t> idx = argsort(p.intensities);
  std::reverse(idx.begin(), idx.end());
  if (idx.size() > static_cast<size_t>(kLibPeakNum)) {
    idx.resize(kLibPeakNum);
  }
  return select(p, idx);
}

// cluster_representative_50pks: np.argsort(inte)[-50:], then by mass
Peaks topPeaksByMass(const Peaks& p) {
  std::vector<size_t> idx = argsort(p.intensities);
  if (idx.size() > static_cast<size_t>(kLibPeakNum)) {
    idx.erase(idx.begin(), idx.end() - kLibPeakNum);
  }
  return sortByMass(select(p, idx));
}

// np.log2(inte) / np.linalg.norm(np.log2(inte))
std::vector<double> logNorm(const std::vector<double>& inte) {
  std::vector<double> v(inte.size());
  for (size_t i = 0; i < inte.size(); i++) {
    v[i] = std::log2(inte[i]);
  }
  double norm = lib_numeric::numpyNorm(v);
  for (double& x : v) {
    x = x / norm;
  }
  return v;
}

// distance_calculation.peak_50_cosine_mass_ppm (err = 0, ppm = 10): one
// minus the summed intensity products of the peak pairs of equal charge
// whose masses differ by at most the tolerance or by one isotope shift
double cosineDistance(const Peaks& a, const Peaks& b) {
  double cos = 0;
  for (size_t i = 0; i < a.masses.size(); i++) {
    double tol = std::max(0.0, kPpm * a.masses[i] * 1e-6);
    double th1 = 0.0 + tol;
    double th2 = kIsotopeShift + tol;
    double th3 = kIsotopeShift - tol;
    for (size_t j = 0; j < b.masses.size(); j++) {
      if (a.charges[i] != b.charges[j]) {
        continue;
      }
      double d = std::fabs(b.masses[j] - a.masses[i]);
      if (d <= th1 || (d <= th2 && d >= th3)) {
        cos += a.intensities[i] * b.intensities[j];
        if (cos > 1) {
          return 1.0 - (cos - a.intensities[i] * b.intensities[j]);
        }
      }
    }
  }
  return 1.0 - cos;
}

// linkage + fcluster(criterion='distance') + cluster_index_gen on the
// members of one group, mapped back to their spectrum indices (sorted).
// Pairs for which cannot_link is true get the distance n * n (n the group
// size) instead of dist_func's: a merge of two clusters A and B with such a
// pair then has an average distance of at least n * n / (|A| |B|) >= 4, so
// with average linkage and a cutoff below 4 no flat cluster holds the pair.
template <typename DistFunc, typename CannotLink>
void clusterGroup(const std::vector<int>& members, lib_cluster::Linkage method,
                  double cutoff, DistFunc dist_func, CannotLink cannot_link,
                  std::vector<std::vector<int>>& clusters) {
  int n = static_cast<int>(members.size());
  if (n <= 1) {
    clusters.push_back(members);
    return;
  }
  const double apart = static_cast<double>(n) * n;
  std::vector<double> dist;
  dist.reserve(static_cast<size_t>(n) * (n - 1) / 2);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      dist.push_back(cannot_link(members[i], members[j])
                         ? apart
                         : dist_func(members[i], members[j]));
    }
  }
  std::vector<std::array<double, 4>> z =
      lib_cluster::linkage(std::move(dist), n, method);
  std::vector<int> labels = lib_cluster::fclusterDistance(z, n, cutoff);
  for (const std::vector<int>& local : lib_cluster::clusterMembers(labels)) {
    std::vector<int> global;
    for (int x : local) {
      global.push_back(members[x]);
    }
    std::sort(global.begin(), global.end());
    clusters.push_back(global);
  }
}

// weighted_mass_sum_inte_cal: every peak of a is merged with the peaks of b
// of the same charge within 10 ppm of it (intensity-weighted mass, summed
// intensity); the unmatched peaks of a, then of b, are appended
Peaks mergeTwo(const Peaks& a, const Peaks& b) {
  Peaks r;
  std::vector<bool> used_a(a.masses.size(), false);
  std::vector<bool> used_b(b.masses.size(), false);
  for (size_t i = 0; i < a.masses.size(); i++) {
    double th = kPpm * a.masses[i] * 1e-6;
    th = 0.0 + std::max(0.0, th);
    double int_sum = 0;
    double w_mass_sum = 0;
    for (size_t j = 0; j < b.masses.size(); j++) {
      if (a.charges[i] != b.charges[j] ||
          std::fabs(b.masses[j] - a.masses[i]) > th) {
        continue;
      }
      if (!used_a[i]) {
        w_mass_sum = w_mass_sum + (a.masses[i] * a.intensities[i] +
                                   b.masses[j] * b.intensities[j]);
        int_sum = int_sum + (a.intensities[i] + b.intensities[j]);
      } else {
        w_mass_sum = w_mass_sum + (b.masses[j] * b.intensities[j]);
        int_sum = int_sum + b.intensities[j];
      }
      used_a[i] = true;
      used_b[j] = true;
    }
    if (int_sum != 0) {
      r.masses.push_back(w_mass_sum / int_sum);
      r.intensities.push_back(int_sum);
      r.charges.push_back(a.charges[i]);
    }
  }
  for (size_t i = 0; i < a.masses.size(); i++) {
    if (!used_a[i]) {
      r.masses.push_back(a.masses[i]);
      r.intensities.push_back(a.intensities[i]);
      r.charges.push_back(a.charges[i]);
    }
  }
  for (size_t j = 0; j < b.masses.size(); j++) {
    if (!used_b[j]) {
      r.masses.push_back(b.masses[j]);
      r.intensities.push_back(b.intensities[j]);
      r.charges.push_back(b.charges[j]);
    }
  }
  return r;
}

// precursor_mass_cluster_rep: the members whose integer precursor mass
// differs from the first member's by the most common difference (the
// smallest such difference on a tie); their intensity-weighted mean mass
// and mean intensity
void precursorRep(const std::vector<double>& masses,
                  const std::vector<double>& intes, double& rep_mass,
                  double& rep_inte) {
  std::vector<long long> diff(masses.size());
  for (size_t i = 0; i < masses.size(); i++) {
    diff[i] = static_cast<long long>(masses[i]) -
              static_cast<long long>(masses[0]);
  }
  std::map<long long, int> counts;
  for (long long d : diff) {
    counts[d]++;
  }
  long long best = counts.begin()->first;
  int best_count = 0;
  for (const auto& c : counts) {
    if (c.second > best_count) {
      best = c.first;
      best_count = c.second;
    }
  }
  std::vector<double> prod;
  std::vector<double> sel;
  for (size_t i = 0; i < masses.size(); i++) {
    if (diff[i] == best) {
      prod.push_back(masses[i] * intes[i]);
      sel.push_back(intes[i]);
    }
  }
  double inte_sum = lib_numeric::numpySum(sel);
  rep_mass = lib_numeric::numpySum(prod) / inte_sum;
  rep_inte = inte_sum / static_cast<double>(sel.size());
}

// the E-value of a spectrum, NaN when it has none
double eValueOf(const LibSpectrum& sp, int e_col) {
  if (e_col < 0 || sp.id_values.empty()) {
    return std::nan("");
  }
  const LibValue& v = sp.id_values[e_col];
  if (const double* d = std::get_if<double>(&v)) {
    return *d;
  }
  if (const long long* l = std::get_if<long long>(&v)) {
    return static_cast<double>(*l);
  }
  return std::nan("");
}

// the first member with the lowest E-value (pandas idxmin, or the first row
// of sort_values('e_value')), or -1 if no member has one
int minEValueMember(const std::vector<int>& members,
                    const std::vector<double>& e_values) {
  int best = -1;
  for (int m : members) {
    double e = e_values[m];
    if (!std::isnan(e) && (best < 0 || e < e_values[best])) {
      best = m;
    }
  }
  return best;
}

int columnIndex(const std::vector<LibIdColumn>& cols, const std::string& name) {
  for (size_t i = 0; i < cols.size(); i++) {
    if (cols[i].name == name) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

// the columns of lib_ptm (TopLib's col_names that exist, 'Cluster ID'
// dropped) with the types pandas gives them
std::vector<LibIdColumn> ptmColumns(const LibData& data) {
  auto optionalInt = [](bool all) {
    return all ? LibColType::kInteger : LibColType::kReal;
  };
  std::vector<LibIdColumn> cols = {{"file_id", LibColType::kInteger},
                                   {"file_name", LibColType::kText},
                                   {"spectrum_id", LibColType::kInteger},
                                   {"title", LibColType::kText},
                                   {"scan", LibColType::kInteger},
                                   {"retention_time", LibColType::kReal}};
  if (data.any_ms_level) {
    cols.push_back({"ms_level", optionalInt(data.all_ms_level)});
  }
  if (data.any_ms_one_id) {
    cols.push_back({"ms_one_id", optionalInt(data.all_ms_one_id)});
  }
  if (data.any_ms_one_scan) {
    cols.push_back({"ms_one_scan", optionalInt(data.all_ms_one_scan)});
  }
  if (data.any_window_begin) {
    cols.push_back({"precursor_window_begin", LibColType::kReal});
  }
  if (data.any_window_end) {
    cols.push_back({"precursor_window_end", LibColType::kReal});
  }
  cols.push_back({"activation", LibColType::kText});
  if (data.any_precursor_mz) {
    cols.push_back({"precursor_mz", LibColType::kReal});
  }
  cols.push_back({"precursor_feature_id", LibColType::kInteger});
  for (const LibIdColumn& c : data.id_columns) {
    cols.push_back(c);
  }
  // an integer column that .loc sets to NaN becomes float64
  cols.push_back({"cluster_id", LibColType::kReal});
  return cols;
}

LibValue optionalValue(const std::optional<long long>& v, bool all) {
  if (!v.has_value()) {
    return LibValue();
  }
  if (all) {
    return LibValue(*v);
  }
  return LibValue(static_cast<double>(*v));
}

LibValue optionalValue(const std::optional<double>& v) {
  return v.has_value() ? LibValue(*v) : LibValue();
}

std::vector<LibValue> ptmValues(const LibData& data, const LibSpectrum& sp,
                                int cluster_id) {
  std::vector<LibValue> v = {LibValue(1LL), LibValue(sp.file_name),
                             LibValue(sp.spectrum_id), LibValue(sp.title),
                             LibValue(sp.scan), LibValue(sp.retention_time)};
  if (data.any_ms_level) {
    v.push_back(optionalValue(sp.ms_level, data.all_ms_level));
  }
  if (data.any_ms_one_id) {
    v.push_back(optionalValue(sp.ms_one_id, data.all_ms_one_id));
  }
  if (data.any_ms_one_scan) {
    v.push_back(optionalValue(sp.ms_one_scan, data.all_ms_one_scan));
  }
  if (data.any_window_begin) {
    v.push_back(optionalValue(sp.precursor_window_begin));
  }
  if (data.any_window_end) {
    v.push_back(optionalValue(sp.precursor_window_end));
  }
  v.push_back(LibValue(sp.activation));
  if (data.any_precursor_mz) {
    v.push_back(optionalValue(sp.precursor_mz));
  }
  v.push_back(LibValue(sp.precursor_feature_id));
  if (sp.id_values.empty()) {
    v.resize(v.size() + data.id_columns.size());
  } else {
    v.insert(v.end(), sp.id_values.begin(), sp.id_values.end());
  }
  v.push_back(LibValue(static_cast<double>(cluster_id)));
  return v;
}

}  // namespace

void filterPpir(LibData& data, double cutoff) {
  if (cutoff == 0) {
    return;
  }
  std::vector<LibSpectrum> kept;
  for (LibSpectrum& sp : data.spectra) {
    const std::vector<double>& p = sp.precursor_intensities;
    if (p.empty() || p[0] == 0) {
      continue;
    }
    if (p[0] / lib_numeric::pythonSum(p) >= cutoff) {
      kept.push_back(std::move(sp));
    }
  }
  data.spectra = std::move(kept);
}

LibLibrary build(const LibData& data, LibRepMethod method,
                 bool python_compatible) {
  LibLibrary lib;
  lib.input_spectrum_num = static_cast<int>(data.spectra.size());

  // ms_incorr_remove: the 50 most intense peaks; at least two peaks and a
  // nonzero precursor intensity
  std::vector<const LibSpectrum*> specs;
  std::vector<Peaks> top_peaks;
  for (const LibSpectrum& sp : data.spectra) {
    Peaks top = topPeaksDescending(peaksOf(sp));
    if (top.masses.size() > 1 && sp.precursor_intensity != 0) {
      specs.push_back(&sp);
      top_peaks.push_back(std::move(top));
    }
  }
  int spec_num = static_cast<int>(specs.size());
  lib.valid_spectrum_num = spec_num;

  // ms_spectrum_preprocess(..., 'log'): by mass, log2, unit length
  std::vector<Peaks> log_peaks(spec_num);
  for (int s = 0; s < spec_num; s++) {
    log_peaks[s] = sortByMass(top_peaks[s]);
    log_peaks[s].intensities = logNorm(log_peaks[s].intensities);
  }

  // clustering_precursor_mass_index: by charge (increasing), then complete
  // linkage on the precursor masses
  std::map<int, std::vector<int>> by_charge;
  for (int s = 0; s < spec_num; s++) {
    by_charge[specs[s]->precursor_charge].push_back(s);
  }
  std::vector<std::vector<int>> prec_clusters;
  for (const auto& group : by_charge) {
    clusterGroup(
        group.second, lib_cluster::Linkage::kComplete, kPrecursorTolerance,
        [&specs](int a, int b) {
          return std::fabs(specs[a]->precursor_mass - specs[b]->precursor_mass);
        },
        [](int, int) { return false; }, prec_clusters);
  }
  lib.precursor_cluster_num = static_cast<int>(prec_clusters.size());

  // clustering_index_50peaks: average linkage on the cosine distances.
  // Unlike the Python script, two spectra identified as different
  // proteoforms (different proteoform ids) are never put in one cluster;
  // unidentified spectra may join any cluster.
  int form_id_col = columnIndex(data.id_columns, "proteoform_id");
  std::vector<LibValue> form_ids(spec_num);
  for (int s = 0; s < spec_num; s++) {
    if (form_id_col >= 0 && !specs[s]->id_values.empty()) {
      form_ids[s] = specs[s]->id_values[form_id_col];
    }
  }
  auto different_forms = [&form_ids, python_compatible](int a, int b) {
    return !python_compatible &&
           !std::holds_alternative<std::monostate>(form_ids[a]) &&
           !std::holds_alternative<std::monostate>(form_ids[b]) &&
           form_ids[a] != form_ids[b];
  };
  std::vector<std::vector<int>> clusters;
  for (const std::vector<int>& pc : prec_clusters) {
    clusterGroup(
        pc, lib_cluster::Linkage::kAverage, kSpectrumDistance,
        [&log_peaks](int a, int b) {
          return cosineDistance(log_peaks[a], log_peaks[b]);
        },
        different_forms, clusters);
  }
  int cluster_num = static_cast<int>(clusters.size());
  for (const std::vector<int>& members : clusters) {
    std::vector<int> positions;
    for (int m : members) {
      positions.push_back(static_cast<int>(specs[m] - data.spectra.data()));
    }
    lib.cluster_members.push_back(positions);
  }

  int e_col = columnIndex(data.id_columns, "e_value");
  std::vector<double> e_values(spec_num);
  for (int s = 0; s < spec_num; s++) {
    e_values[s] = eValueOf(*specs[s], e_col);
  }

  // the representatives
  lib.reps.resize(cluster_num);
  for (int k = 0; k < cluster_num; k++) {
    const std::vector<int>& members = clusters[k];
    LibRep& rep = lib.reps[k];
    Peaks peaks;
    if (method == LibRepMethod::kAverage) {
      // cluster_representative_50pks on all peaks, sorted by mass
      peaks = sortByMass(peaksOf(*specs[members[0]]));
      const LibSpectrum& first = *specs[members[0]];
      if (members.size() > 1) {
        for (size_t m = 1; m < members.size(); m++) {
          Peaks next = sortByMass(peaksOf(*specs[members[m]]));
          peaks = sortByMass(mergeTwo(peaks, next));
        }
        std::vector<double> masses;
        std::vector<double> intes;
        for (int m : members) {
          masses.push_back(specs[m]->precursor_mass);
          intes.push_back(specs[m]->precursor_intensity);
        }
        precursorRep(masses, intes, rep.precursor_mass, rep.precursor_intensity);
      } else {
        rep.precursor_mass = first.precursor_mass;
        rep.precursor_intensity = first.precursor_intensity;
      }
      rep.precursor_charge = first.precursor_charge;
      peaks = topPeaksByMass(peaks);
    } else {
      // cluster_representative_min_evalue: the member with the lowest
      // E-value, or the most masses, and its 50 most intense peaks
      int best = minEValueMember(members, e_values);
      if (best < 0) {
        best = members[0];
        for (int m : members) {
          if (specs[m]->masses.size() > specs[best]->masses.size()) {
            best = m;
          }
        }
      }
      peaks = sortByMass(top_peaks[best]);
      rep.precursor_mass = specs[best]->precursor_mass;
      rep.precursor_intensity = specs[best]->precursor_intensity;
      rep.precursor_charge = specs[best]->precursor_charge;
    }
    rep.masses = std::move(peaks.masses);
    rep.intensities = std::move(peaks.intensities);
    rep.charges = std::move(peaks.charges);
    rep.norm_intensities = logNorm(rep.intensities);
  }

  // lib_ptm: the identification of each cluster's lowest-E-value member
  lib.ptm_columns = ptmColumns(data);
  int form_col = columnIndex(data.id_columns, "proteoform");
  std::vector<int> ptm_member(cluster_num);
  std::map<std::pair<int, LibValue>, std::vector<int>> by_form;
  for (int k = 0; k < cluster_num; k++) {
    int best = minEValueMember(clusters[k], e_values);
    if (best < 0) {
      best = clusters[k][0];
    }
    ptm_member[k] = best;
    const LibSpectrum& sp = *specs[best];
    lib.reps[k].ptm_values = ptmValues(data, sp, k);
    bool has_form = !sp.id_values.empty() && form_col >= 0 &&
                    !std::holds_alternative<std::monostate>(sp.id_values[form_col]);
    lib.reps[k].identified = has_form;
    if (!has_form) {
      lib.unidentified_cluster_num++;
      continue;
    }
    // groupby(['precursor_charge', 'proteoform_id']) leaves out a missing id
    if (form_id_col < 0 ||
        std::holds_alternative<std::monostate>(sp.id_values[form_id_col])) {
      continue;
    }
    by_form[{sp.precursor_charge, sp.id_values[form_id_col]}].push_back(k);
  }
  // remove inconsistent identifications: of the clusters with the same
  // proteoform and precursor charge only the one with the lowest E-value
  // keeps it
  for (const auto& group : by_form) {
    const std::vector<int>& ks = group.second;
    if (ks.size() <= 1) {
      continue;
    }
    std::vector<int> order(ks);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
      double ea = e_values[ptm_member[a]];
      double eb = e_values[ptm_member[b]];
      if (std::isnan(ea) || std::isnan(eb)) {
        return !std::isnan(ea) && std::isnan(eb);
      }
      return ea < eb;
    });
    for (size_t i = 1; i < order.size(); i++) {
      lib.reps[order[i]].identified = false;
      lib.inconsistent_cluster_num++;
    }
  }
  // the columns from the proteoform id to the cluster id become NaN
  int first_replaced = static_cast<int>(lib.ptm_columns.size()) -
                       static_cast<int>(data.id_columns.size()) - 1;
  for (LibRep& rep : lib.reps) {
    if (!rep.identified) {
      for (size_t c = first_replaced; c < rep.ptm_values.size(); c++) {
        rep.ptm_values[c] = LibValue();
      }
    }
  }

  // decoy_mass_gen: every mass shifted by the difference of two random
  // residue masses, redrawn while the shifted mass is negative, then while
  // it exceeds the precursor mass
  for (int k = 0; k < cluster_num; k++) {
    LibRep& rep = lib.reps[k];
    lib_numeric::PyRandom rng(static_cast<uint64_t>(k));
    auto shifted = [&rng](double mass) {
      std::vector<int> aa = rng.sample(19, 2);
      return mass + (kAaMasses[aa[0]] - kAaMasses[aa[1]]);
    };
    rep.decoy_masses.resize(rep.masses.size());
    for (size_t i = 0; i < rep.masses.size(); i++) {
      double m = shifted(rep.masses[i]);
      while (m < 0) {
        m = shifted(rep.masses[i]);
      }
      while (m > rep.precursor_mass) {
        m = shifted(rep.masses[i]);
      }
      rep.decoy_masses[i] = m;
    }
  }
  return lib;
}

}  // namespace lib_builder

}  // namespace toplib

}  // namespace toppic
