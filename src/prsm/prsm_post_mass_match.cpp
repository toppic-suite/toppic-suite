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

#include "prsm/prsm_post_mass_match.hpp"

#include <sqlite3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "common/base/activation.hpp"
#include "common/base/ion.hpp"
#include "common/base/ion_type.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/env/env.hpp"
#include "ms/env/env_base.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/deconv_peak.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/spec/theo_peak.hpp"
#include "para/peak_tolerance.hpp"
#include "para/sp_para.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/theo_peak_util.hpp"
#include "sql/sql_util.hpp"

namespace toppic {

namespace prsm_post_mass_match {

namespace {

// Matching parameters of MSPathFinderT (CompositeScorer and
// AbstractFragmentScorer in Informed-Proteomics): theoretical isotopic peaks
// below 10% of the highest one are neither looked for nor scored; an envelope
// is accepted with a Pearson correlation >= 0.7 or a Bhattacharyya distance
// <= 0.03; fragment charges run from 1 to min(15, precursor charge + 2).
constexpr double MIN_RELATIVE_INTE = 0.1;
constexpr double MIN_CORR = 0.7;
constexpr double MAX_DIST = 0.03;
constexpr int MAX_CHARGE = 15;
constexpr int CHARGE_ABOVE_PREC = 2;

struct RawPeak {
  int id;
  double mz;
  double inte;
};

using RawPeakVec = std::vector<RawPeak>;

// One theoretical fragment mass found in the centroided peaks.
struct PostMatch {
  double mono_mass;  // experimental monoisotopic mass
  double inte;       // summed intensity of the observed isotopic peaks
  double apex_inte;  // intensity of the observed most abundant isotopic peak
  double apex_mz;
  int charge;
  double corr;
  int peak_num;  // observed isotopic peaks
  EnvPtr env;    // theoretical envelope shifted and scaled to the apex peak
};

using PostMatchVec = std::vector<PostMatch>;

// Centroided MS/MS peaks of the TopFD SQLite database, read by scan number.
class RawPeakReader {
 public:
  explicit RawPeakReader(sqlite3* db) : db_(db) {
    id_stmt_ = sql_util::prepareSql(
        db_, "SELECT id FROM ms2_spectrum WHERE scan = ?;");
    peak_stmt_ = sql_util::prepareSql(
        db_,
        "SELECT peak_id, mz, intensity FROM ms2_peak WHERE spec_id = ? "
        "ORDER BY mz;");
  }

  ~RawPeakReader() {
    sqlite3_finalize(id_stmt_);
    sqlite3_finalize(peak_stmt_);
  }

  // The database id of the spectrum with the scan number, -1 if absent.
  int getSpecId(int scan) {
    sqlite3_bind_int(id_stmt_, 1, scan);
    int spec_id = -1;
    if (sqlite3_step(id_stmt_) == SQLITE_ROW) {
      spec_id = sqlite3_column_int(id_stmt_, 0);
    }
    sqlite3_clear_bindings(id_stmt_);
    sqlite3_reset(id_stmt_);
    return spec_id;
  }

  RawPeakVec getPeaks(int spec_id) {
    RawPeakVec peaks;
    sqlite3_bind_int(peak_stmt_, 1, spec_id);
    while (sqlite3_step(peak_stmt_) == SQLITE_ROW) {
      peaks.push_back({sqlite3_column_int(peak_stmt_, 0),
                       sqlite3_column_double(peak_stmt_, 1),
                       sqlite3_column_double(peak_stmt_, 2)});
    }
    sqlite3_clear_bindings(peak_stmt_);
    sqlite3_reset(peak_stmt_);
    return peaks;
  }

 private:
  sqlite3* db_;
  sqlite3_stmt* id_stmt_;
  sqlite3_stmt* peak_stmt_;
};

// Index of the most intense peak within mz +- tole, -1 if there is none.
int findPeak(const RawPeakVec& peaks, double mz, double tole) {
  auto it = std::lower_bound(
      peaks.begin(), peaks.end(), mz - tole,
      [](const RawPeak& peak, double value) { return peak.mz < value; });
  int best = -1;
  for (; it != peaks.end() && it->mz <= mz + tole; ++it) {
    int idx = static_cast<int>(it - peaks.begin());
    if (best < 0 || it->inte > peaks[best].inte) {
      best = idx;
    }
  }
  return best;
}

// Bhattacharyya distance and Pearson correlation between the theoretical and
// the observed isotopic intensities (FitScoreCalculator.GetDistanceAndCorrelation
// in Informed-Proteomics): the distance is computed on the two vectors
// normalized to sum 1; the correlation is 0 unless the covariance is positive.
std::pair<double, double> compDistCorr(const std::vector<double>& theo,
                                       const std::vector<double>& obs) {
  size_t n = theo.size();
  double s1 = 0;
  double s2 = 0;
  for (size_t i = 0; i < n; i++) {
    s1 += theo[i];
    s2 += obs[i];
  }
  if (n == 0 || !(s1 > 0) || !(s2 > 0)) {
    return {1.0, 0.0};
  }
  double m1 = s1 / n;
  double m2 = s2 / n;
  double cov = 0;
  double c1 = 0;
  double c2 = 0;
  double bc = 0;
  for (size_t i = 0; i < n; i++) {
    double d1 = theo[i] - m1;
    double d2 = obs[i] - m2;
    cov += d1 * d2;
    c1 += d1 * d1;
    c2 += d2 * d2;
    bc += std::sqrt((theo[i] / s1) * (obs[i] / s2));
  }
  double corr = 0;
  if (c1 > 0 && c2 > 0 && cov > 0) {
    corr = cov / std::sqrt(c1 * c2);
  }
  double dist = -std::log(bc);
  return {dist, corr};
}

// Look for the theoretical envelope of a fragment mass at one charge state in
// the centroided peaks. Returns false if the envelope is not accepted.
bool matchTheoMass(const RawPeakVec& peaks, double mass, int charge,
                   double ppo, int min_peak_num, PostMatch& match) {
  EnvPtr env_ptr = EnvBase::getEnvByMonoMass(mass, charge);
  if (env_ptr == nullptr) {
    return false;
  }
  int refer_idx = env_ptr->getReferIdx();
  double min_inte = MIN_RELATIVE_INTE * env_ptr->getReferInte();
  // the isotopic peaks looked for: contiguous around the most abundant one
  int left = refer_idx;
  while (left > 0 && env_ptr->getInte(left - 1) >= min_inte) {
    left--;
  }
  int right = refer_idx;
  while (right + 1 < env_ptr->getPeakNum() &&
         env_ptr->getInte(right + 1) >= min_inte) {
    right++;
  }
  double refer_mz = env_ptr->getMz(refer_idx);
  int apex_idx = findPeak(peaks, refer_mz, refer_mz * ppo);
  if (apex_idx < 0) {
    return false;
  }
  std::vector<double> theo;
  std::vector<double> obs;
  int peak_num = 0;
  double inte_sum = 0;
  for (int i = left; i <= right; i++) {
    double mz = env_ptr->getMz(i);
    int idx = (i == refer_idx) ? apex_idx : findPeak(peaks, mz, mz * ppo);
    theo.push_back(env_ptr->getInte(i));
    if (idx >= 0) {
      obs.push_back(peaks[idx].inte);
      inte_sum += peaks[idx].inte;
      peak_num++;
    } else {
      obs.push_back(0.0);
    }
  }
  if (peak_num < min_peak_num) {
    return false;
  }
  std::pair<double, double> dist_corr = compDistCorr(theo, obs);
  if (dist_corr.second < MIN_CORR && dist_corr.first > MAX_DIST) {
    return false;
  }
  double apex_mz = peaks[apex_idx].mz;
  // keep the isotopic peaks looked for, move the envelope onto the observed
  // apex peak and scale it to the apex intensity, as TopFD does for its
  // envelopes
  env_ptr = env_ptr->getSubEnv(refer_idx - left, right - refer_idx);
  env_ptr->changeMz(apex_mz - refer_mz);
  env_ptr->changeToAbsInte(peaks[apex_idx].inte);
  match.mono_mass = env_ptr->getMonoNeutralMass();
  match.inte = inte_sum;
  match.apex_inte = peaks[apex_idx].inte;
  match.apex_mz = apex_mz;
  match.charge = charge;
  match.corr = dist_corr.second;
  match.peak_num = peak_num;
  match.env = env_ptr;
  return true;
}

// The theoretical masses of a PrSM that none of its deconvoluted masses
// matched, keyed by ion type (N-terminal or not) and position.
TheoPeakPtrVec getUnmatchedTheoPeaks(const PrsmPtr& prsm_ptr,
                                     const ActivationPtr& activation_ptr,
                                     double min_mass) {
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  TheoPeakPtrVec theo_peaks =
      theo_peak_util::geneProteoformTheoPeak(form_ptr, activation_ptr,
                                             min_mass);
  PeakIonPairPtrVec pairs = peak_ion_pair_util::genePeakIonPairs(
      form_ptr, prsm_ptr->getRefineMsPtrVec(), min_mass);
  std::set<std::pair<bool, int>> matched;
  for (size_t i = 0; i < pairs.size(); i++) {
    IonPtr ion_ptr = pairs[i]->getTheoPeakPtr()->getIonPtr();
    matched.insert({ion_ptr->getIonTypePtr()->isNTerm(), ion_ptr->getPos()});
  }
  TheoPeakPtrVec unmatched;
  for (size_t i = 0; i < theo_peaks.size(); i++) {
    IonPtr ion_ptr = theo_peaks[i]->getIonPtr();
    if (matched.count({ion_ptr->getIonTypePtr()->isNTerm(),
                       ion_ptr->getPos()}) == 0) {
      unmatched.push_back(theo_peaks[i]);
    }
  }
  return unmatched;
}

// Match the unmatched theoretical masses of one PrSM against the centroided
// peaks of its spectrum. Every theoretical mass adds at most one experimental
// mass: the accepted charge state with the most intense apex peak.
PostMatchVec matchPrsm(const PrsmPtr& prsm_ptr, const RawPeakVec& peaks,
                       const SpParaPtr& sp_para_ptr, int min_peak_num) {
  MsHeaderPtr header_ptr =
      prsm_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr();
  TheoPeakPtrVec theo_peaks = getUnmatchedTheoPeaks(
      prsm_ptr, header_ptr->getActivationPtr(), sp_para_ptr->getMinMass());
  int max_charge =
      std::min(MAX_CHARGE, header_ptr->getFirstPrecCharge() + CHARGE_ABOVE_PREC);
  double ppo = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  PostMatchVec matches;
  for (size_t i = 0; i < theo_peaks.size(); i++) {
    double mass = theo_peaks[i]->getModMass();
    bool found = false;
    PostMatch best;
    for (int charge = 1; charge <= max_charge; charge++) {
      PostMatch match;
      if (matchTheoMass(peaks, mass, charge, ppo, min_peak_num, match) &&
          (!found || match.apex_inte > best.apex_inte)) {
        best = match;
        found = true;
      }
    }
    if (found) {
      matches.push_back(best);
    }
  }
  return matches;
}

// Append the matched masses to the PrSM's spectrum and recount its matched
// masses and fragments.
void addMatchesToPrsm(const PrsmPtr& prsm_ptr, const PostMatchVec& matches,
                      const SpParaPtr& sp_para_ptr) {
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  DeconvMsPtr ms_ptr = deconv_ms_ptr_vec[0];
  int sp_id = ms_ptr->getMsHeaderPtr()->getSpecId();
  DeconvPeakPtrVec peaks = ms_ptr->getPeakPtrVec();
  for (size_t i = 0; i < matches.size(); i++) {
    peaks.push_back(std::make_shared<DeconvPeak>(
        sp_id, static_cast<int>(peaks.size()), matches[i].mono_mass,
        matches[i].inte, matches[i].charge, matches[i].corr));
  }
  std::sort(peaks.begin(), peaks.end(), DeconvPeak::cmpPosInc);
  ms_ptr->setPeakPtrVec(peaks);
  ExtendMsPtrVec extend_ms_ptr_vec = extend_ms_factory::geneMsThreePtrVec(
      deconv_ms_ptr_vec, sp_para_ptr, prsm_ptr->getAdjustedPrecMass());
  prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
  prsm_ptr->updateMatchNum(sp_para_ptr);
}

// The parameter block (comment lines) at the top of an msalign file.
std::string readMsalignPara(const std::string& file_name) {
  std::ifstream input(file_name);
  std::string line;
  std::string para_str;
  while (std::getline(input, line)) {
    if (line == "BEGIN IONS") {
      break;
    }
    if (!para_str.empty()) {
      para_str += "\n";
    }
    para_str += line;
  }
  return para_str;
}

// Write the spectra of the input msalign file, with the matched masses
// appended, to the new msalign file.
void writeMsalign(const std::string& sp_file_name,
                  const std::string& post_sp_file_name,
                  const SpParaPtr& sp_para_ptr,
                  const std::map<int, PostMatchVec>& spec_matches) {
  MsAlignWriter writer(post_sp_file_name);
  writer.writePara(readMsalignPara(sp_file_name));
  MsAlignReader reader(sp_file_name, 1, sp_para_ptr->getActivationPtr());
  DeconvMsPtr ms_ptr = reader.getNextMsPtr();
  while (ms_ptr != nullptr) {
    int sp_id = ms_ptr->getMsHeaderPtr()->getSpecId();
    auto it = spec_matches.find(sp_id);
    if (it != spec_matches.end()) {
      DeconvPeakPtrVec peaks = ms_ptr->getPeakPtrVec();
      for (size_t i = 0; i < it->second.size(); i++) {
        const PostMatch& match = it->second[i];
        peaks.push_back(std::make_shared<DeconvPeak>(
            sp_id, static_cast<int>(peaks.size()), match.mono_mass,
            match.inte, match.charge, match.corr));
      }
      ms_ptr->setPeakPtrVec(peaks);
    }
    writer.writeMs(ms_ptr);
    ms_ptr = reader.getNextMsPtr();
  }
}

// Add the matched masses as envelopes to the ms2_env and ms2_env_peak tables,
// continuing the envelope ids of each spectrum.
void writeSqlEnvs(sqlite3* db,
                  const std::map<int, PostMatchVec>& sql_spec_matches) {
  sqlite3_stmt* max_stmt = sql_util::prepareSql(
      db, "SELECT COALESCE(MAX(env_id) + 1, 0) FROM ms2_env WHERE spec_id = ?;");
  sqlite3_stmt* env_stmt = sql_util::prepareSql(
      db,
      "INSERT INTO ms2_env(spec_id, env_id, mono_mass, ref_mass, charge, "
      "intensity, envcnn_score, peak_num) VALUES (?, ?, ?, ?, ?, ?, ?, ?);");
  sqlite3_stmt* peak_stmt = sql_util::prepareSql(
      db,
      "INSERT INTO ms2_env_peak(spec_id, env_id, peak_id, mz, intensity) "
      "VALUES (?, ?, ?, ?, ?);");
  sql_util::execSql(db, "BEGIN TRANSACTION;");
  for (const auto& entry : sql_spec_matches) {
    int spec_id = entry.first;
    sqlite3_bind_int(max_stmt, 1, spec_id);
    int env_id = 0;
    if (sqlite3_step(max_stmt) == SQLITE_ROW) {
      env_id = sqlite3_column_int(max_stmt, 0);
    }
    sqlite3_clear_bindings(max_stmt);
    sqlite3_reset(max_stmt);
    for (const PostMatch& match : entry.second) {
      EnvPtr env_ptr = match.env;
      int peak_num = env_ptr->getPeakNum();
      sqlite3_bind_int(env_stmt, 1, spec_id);
      sqlite3_bind_int(env_stmt, 2, env_id);
      sqlite3_bind_double(env_stmt, 3, match.mono_mass);
      sqlite3_bind_double(env_stmt, 4, env_ptr->getReferNeutralMass());
      sqlite3_bind_int(env_stmt, 5, match.charge);
      sqlite3_bind_double(env_stmt, 6, env_ptr->compInteSum());
      sqlite3_bind_double(env_stmt, 7, match.corr);
      sqlite3_bind_int(env_stmt, 8, peak_num);
      sql_util::stepAndReset(db, env_stmt);
      for (int k = 0; k < peak_num; k++) {
        sqlite3_bind_int(peak_stmt, 1, spec_id);
        sqlite3_bind_int(peak_stmt, 2, env_id);
        sqlite3_bind_int(peak_stmt, 3, k);
        sqlite3_bind_double(peak_stmt, 4, env_ptr->getMz(k));
        sqlite3_bind_double(peak_stmt, 5, env_ptr->getInte(k));
        sql_util::stepAndReset(db, peak_stmt);
      }
      env_id++;
    }
  }
  sql_util::execSql(db, "END TRANSACTION;");
  sqlite3_finalize(max_stmt);
  sqlite3_finalize(env_stmt);
  sqlite3_finalize(peak_stmt);
}

}  // namespace

std::string process(const PrsmParaPtr& prsm_para_ptr,
                    const std::string& input_file_ext,
                    const std::string& output_file_ext, int min_peak_num) {
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string base_name = file_util::basename(sp_file_name);
  std::string sql_base = base_name;
  if (str_util::endsWith(sql_base, "_ms2")) {
    sql_base = sql_base.substr(0, sql_base.size() - 4);
  }
  std::string post_base_name = sql_base + "_post_ms2";
  std::string post_sp_file_name = post_base_name + ".msalign";
  std::string sql_file_name = sql_base + ".sqlite";
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();

  PrsmPtrVec prsm_ptrs =
      prsm_reader_util::readPrsmsWithSpectra(prsm_para_ptr, input_file_ext);

  // spectrum id -> matched masses (for the msalign file); database spectrum
  // id -> matched masses (for the SQLite envelope tables)
  std::map<int, PostMatchVec> spec_matches;
  std::map<int, PostMatchVec> sql_spec_matches;
  int match_num = 0;
  int prsm_num = 0;
  if (!std::filesystem::exists(sql_file_name)) {
    std::cout << "SQLite database " << sql_file_name
              << " not found (run topfd without -N): no centroided peaks, "
                 "post mass matching skipped."
              << std::endl;
  } else {
    sqlite3* db = nullptr;
    if (sqlite3_open(sql_file_name.c_str(), &db) != SQLITE_OK) {
      LOG_ERROR("Cannot open the database " << sql_file_name << ": "
                                            << sqlite3_errmsg(db));
      exit(EXIT_FAILURE);
    }
    {
      RawPeakReader reader(db);
      for (size_t i = 0; i < prsm_ptrs.size(); i++) {
        PrsmPtr prsm_ptr = prsm_ptrs[i];
        MsHeaderPtr header_ptr =
            prsm_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr();
        int sql_spec_id = reader.getSpecId(header_ptr->getFirstScanNum());
        if (sql_spec_id < 0) {
          LOG_WARN("Scan " << header_ptr->getFirstScanNum()
                           << " is not in the SQLite database.");
          continue;
        }
        RawPeakVec peaks = reader.getPeaks(sql_spec_id);
        PostMatchVec matches =
            matchPrsm(prsm_ptr, peaks, sp_para_ptr, min_peak_num);
        if (matches.empty()) {
          continue;
        }
        addMatchesToPrsm(prsm_ptr, matches, sp_para_ptr);
        int sp_id = header_ptr->getSpecId();
        PostMatchVec& spec_vec = spec_matches[sp_id];
        spec_vec.insert(spec_vec.end(), matches.begin(), matches.end());
        PostMatchVec& sql_vec = sql_spec_matches[sql_spec_id];
        sql_vec.insert(sql_vec.end(), matches.begin(), matches.end());
        match_num += static_cast<int>(matches.size());
        prsm_num++;
      }
      writeSqlEnvs(db, sql_spec_matches);
    }
    sqlite3_close(db);
  }
  std::cout << "Post mass matching: " << match_num << " masses added to "
            << prsm_num << " of " << prsm_ptrs.size() << " PrSMs." << std::endl;

  // PrSMs with the recounted matched masses and fragments, under the name of
  // the new spectrum file
  std::string post_sp_short_name =
      file_util::filenameFromEntirePath(post_sp_file_name);
  PrsmXmlWriter prsm_writer(post_base_name + "." + output_file_ext);
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    prsm_ptrs[i]->setFileName(post_sp_short_name);
    prsm_writer.write(prsm_ptrs[i]);
  }
  prsm_writer.close();

  writeMsalign(sp_file_name, post_sp_file_name, sp_para_ptr, spec_matches);

  // the downstream steps look for the TopFD feature file under the new name
  std::string feature_file_name = base_name + ".feature";
  if (std::filesystem::exists(feature_file_name)) {
    std::filesystem::copy_file(
        feature_file_name, post_base_name + ".feature",
        std::filesystem::copy_options::overwrite_existing);
  }
  return post_sp_file_name;
}

}  // namespace prsm_post_mass_match

}  // namespace toppic
