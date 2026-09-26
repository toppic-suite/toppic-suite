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

#include "topfd/ecscore/env_coll/env_coll_detect.hpp"

#include <limits>

#include "topfd/ecscore/score/comp_env_cnn_score.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

#include <algorithm>
#include <ctime>
#include <iostream>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/feature/frac_feature_writer.hpp"
#include "ms/feature/spec_feature_writer.hpp"
#include "ms/msmap/ms_map.hpp"
#include "ms/mzml/mzml_ms_group_reader.hpp"
#include "ms/spec/deconv_ms.hpp"
#include <cstdlib>

#include "ms/spec/ms_header.hpp"
#include "ms/spec/msalign_reader_util.hpp"
#include "topfd/common/topfd_para.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"
#include "topfd/ecscore/env_coll/env_coll_assign.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"
#include "topfd/ecscore/score/ecscore.hpp"
#include "topfd/ecscore/score/ecscore_writer.hpp"

namespace toppic {

namespace env_coll_detect {


// Windowed variant of processMs1 (enabled with TOPFD_MS1_WINDOW=1): the raw
// peaks are kept only as compact m/z + intensity arrays, and the peak matrix
// holds materialized rows only for a sliding retention-time window. Seeds are
// processed window by window (intensity order inside a window) instead of in
// one global intensity order, so results can differ slightly from the
// whole-run implementation near window boundaries.
void processMs1Windowed(const TopfdParaPtr& topfd_para_ptr) {
  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(
      topfd_para_ptr->getFracId(), topfd_para_ptr->getMzmlFileName(),
      topfd_para_ptr->getMaxCharge(), topfd_para_ptr->getMs1MinScanNum());
  score_para_ptr->seed_env_inte_corr_tole_cutoff_ =
      topfd_para_ptr->getSeedEnvInteCorrToleCutoff();
  // read deconvoluted MS1 peaks
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();
  std::string ms1_file_name = output_base_name + "_ms1.msalign";
  DeconvMsPtrVec deconv_ms1_ptr_vec;
  msalign_reader_util::readAllSpectra(ms1_file_name, deconv_ms1_ptr_vec);

  // read ms1 raw peaks incrementally into compact arrays and ms2 headers
  std::vector<std::vector<double>> mz_2d;
  std::vector<std::vector<double>> inte_2d;
  MsHeaderPtr2D ms2_header_ptr_2d;
  {
    MzmlMsGroupReaderPtr mzml_reader_ptr = std::make_shared<MzmlMsGroupReader>(
        topfd_para_ptr->getMzmlFileName(), topfd_para_ptr->getPrecWindowWidth(),
        topfd_para_ptr->getActivation(), topfd_para_ptr->getFracId(),
        topfd_para_ptr->isFaims(), topfd_para_ptr->getFaimsVoltage(),
        topfd_para_ptr->isMissingLevelOne());
    while (true) {
      MzmlMsGroupPtr group_ptr = mzml_reader_ptr->getNextMsGroupPtr();
      if (group_ptr == nullptr) {
        break;
      }
      const PeakPtrVec& peak_list = group_ptr->getMsOnePtr()->getPeakPtrVec();
      std::vector<double> mz_list;
      std::vector<double> inte_list;
      mz_list.reserve(peak_list.size());
      inte_list.reserve(peak_list.size());
      for (auto& peak : peak_list) {
        mz_list.push_back(peak->getPosition());
        inte_list.push_back(peak->getIntensity());
      }
      mz_2d.push_back(std::move(mz_list));
      inte_2d.push_back(std::move(inte_list));
      const MzmlMsPtrVec& ms2_ptr_vec = group_ptr->getMsTwoPtrVec();
      MsHeaderPtrVec header_vec;
      for (size_t i = 0; i < ms2_ptr_vec.size(); i++) {
        header_vec.push_back(ms2_ptr_vec[i]->getMsHeaderPtr());
      }
      ms2_header_ptr_2d.push_back(header_vec);
    }
  }
  size_t ms1_raw_spec_num = mz_2d.size();

  // Prepare seed envelopes
  SeedEnvPtrVec seed_ptrs;
  SeedEnvPtr2D seed_ptr_2d;
  for (auto& ms1_data : deconv_ms1_ptr_vec) {
    DeconvPeakPtrVec peaks = ms1_data->getPeakPtrVec();
    SeedEnvPtrVec one_spec_seed_ptrs;
    for (auto& peak : peaks) {
      SeedEnvPtr seed_ptr = std::make_shared<SeedEnv>(peak);
      seed_ptrs.push_back(seed_ptr);
      one_spec_seed_ptrs.push_back(seed_ptr);
    }
    seed_ptr_2d.push_back(one_spec_seed_ptrs);
  }
  std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnv::cmpSeedInteDec);

  double sn_ratio = topfd_para_ptr->getMsOneSnRatio();
  bool single_scan_noise = topfd_para_ptr->isUseSingleScanNoiseLevel();
  MsHeaderPtrVec ms1_header_ptr_vec;
  for (auto& ms1_spectrum : deconv_ms1_ptr_vec) {
    ms1_header_ptr_vec.push_back(ms1_spectrum->getMsHeaderPtr());
  }
  // Peak matrix shell: global statistics, hollow rows
  MsMapPtr matrix_ptr = std::make_shared<MsMap>(
      mz_2d, inte_2d, ms1_header_ptr_vec, score_para_ptr->bin_size_);
  double neighbor_tole = -1;
  if (score_para_ptr->min_scan_num_ >= 2) {
    neighbor_tole = score_para_ptr->neighbor_mz_tole_;
  }

  const int win_rows = 150;
  const int win_margin = 150;
  int row_num = matrix_ptr->getRowNum();
  int num_wins = (row_num + win_rows - 1) / win_rows;
  int max_live_rows = 0;
  int filled_lo = 0;
  int filled_hi = 0;
  double fill_sn_ratio = sn_ratio;
  double fill_neighbor_tole = neighbor_tole;
  auto slide = [&](int lo, int hi) {
    while (filled_lo < lo && filled_lo < filled_hi) {
      matrix_ptr->evictRow(filled_lo);
      filled_lo++;
    }
    if (filled_lo >= filled_hi) {
      filled_lo = lo;
      filled_hi = lo;
    }
    while (filled_hi < hi) {
      matrix_ptr->fillRow(filled_hi, mz_2d[filled_hi], inte_2d[filled_hi],
                          fill_sn_ratio, single_scan_noise,
                          fill_neighbor_tole);
      filled_hi++;
    }
    if (filled_hi - filled_lo > max_live_rows) {
      max_live_rows = filled_hi - filled_lo;
    }
  };

  // bucket seeds by window of their apex spectrum; the global sort keeps
  // intensity order inside each bucket
  std::vector<SeedEnvPtrVec> win_seed_ptrs(num_wins);
  for (auto& seed_ptr : seed_ptrs) {
    int w = seed_ptr->getSpecId() / win_rows;
    if (w >= 0 && w < num_wins) {
      win_seed_ptrs[w].push_back(seed_ptr);
    }
  }

  /// Extract Features
  env_coll_util::EnvCollRtCache rt_cache;
  LOG_DEBUG("Number of seed envelopes: " << seed_ptrs.size());
  int seed_num = seed_ptrs.size();
  int count = 0;
  EnvCollPtrVec env_coll_list;
  ECScorePtrVec ecscore_list;
  for (int w = 0; w < num_wins; w++) {
    slide(std::max(0, w * win_rows - win_margin),
          std::min(row_num, (w + 1) * win_rows + win_margin));
    for (auto& sorted_seed_ptr : win_seed_ptrs[w]) {
      count = count + 1;
      if (count % 100 == 0 || count == seed_num) {
        double perc = static_cast<int>(count * 100 / seed_num);
        std::cout << "\r" << "Processing feature " << count << " ...       "
                  << perc << "\% finished." << std::flush;
      }
      SeedEnvPtr seed_ptr = seed_env_util::preprocessSeedEnvPtr(
          sorted_seed_ptr, matrix_ptr, score_para_ptr, sn_ratio);
      if (seed_ptr == nullptr) continue;
      EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(
          matrix_ptr, seed_ptr, score_para_ptr, sn_ratio,
          topfd_para_ptr->getSplitIntensityRatio());
      if (env_coll_ptr == nullptr) continue;
      if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                               env_coll_list, score_para_ptr,
                                               sn_ratio, &rt_cache)) {
        env_coll_ptr->removePeakData(matrix_ptr);
        continue;
      }
      env_coll_ptr->refineMonoMass();
      ECScorePtr ecscore_ptr =
          std::make_shared<ECScore>(env_coll_ptr, matrix_ptr, sn_ratio);
      if (ecscore_ptr->getScore() < topfd_para_ptr->getMs1EcscoreCutoff()) {
        continue;
      }
      env_coll_ptr->setEcscore(ecscore_ptr->getScore());
      env_coll_ptr->removePeakData(matrix_ptr);
      env_coll_list.push_back(env_coll_ptr);
      ecscore_list.push_back(ecscore_ptr);
    }
  }
  std::cout << std::endl;

  if (topfd_para_ptr->isSearchPrecWindow()) {
    // add ms1 feature based on precursor windows
    // set min match envelope to 1 to accept single scan features
    score_para_ptr->min_scan_num_ = 1;
    score_para_ptr->min_match_peak_ = 1;
    sn_ratio = 0;
    rt_cache.clear();
    // rematerialize rows with sn 0 and no neighbor filter, window by window
    // (the whole-run implementation reconstructs the full matrix here; a
    // refilled row restores the original peak intensities, while reconstruct
    // keeps the reductions made by removePeakData)
    while (filled_lo < filled_hi) {
      matrix_ptr->evictRow(filled_lo);
      filled_lo++;
    }
    filled_lo = 0;
    filled_hi = 0;
    fill_sn_ratio = 0;
    fill_neighbor_tole = -1;
    int ms1_spec_num = deconv_ms1_ptr_vec.size();
    for (std::size_t ms1_idx = 0; ms1_idx < deconv_ms1_ptr_vec.size();
         ms1_idx++) {
      slide(std::max(0, static_cast<int>(ms1_idx) - win_margin),
            std::min(row_num, static_cast<int>(ms1_idx) + win_margin + 1));
      double perc = static_cast<int>((ms1_idx + 1) * 100 / ms1_spec_num);
      int scan =
          deconv_ms1_ptr_vec[ms1_idx]->getMsHeaderPtr()->getFirstScanNum();
      std::cout << "\r" << "Additional feature search MS1 spectrum scan "
                << scan << " ...       " << perc << "\% finished."
                << std::flush;
      for (std::size_t i = 0; i < ms2_header_ptr_2d[ms1_idx].size(); i++) {
        MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_idx][i];
        if (env_coll_assign::checkEnvColl(header_ptr, env_coll_list)) {
          continue;
        }
        double prec_win_begin = header_ptr->getPrecWinBegin();
        double prec_win_end = header_ptr->getPrecWinEnd();
        SeedEnvPtrVec seed_ptr_list = seed_ptr_2d[ms1_idx];
        SeedEnvPtrVec selected_seed_list;
        for (std::size_t j = 0; j < seed_ptr_list.size(); j++) {
          double ref_mz = seed_ptr_list[j]->getReferMz();
          if (ref_mz > prec_win_begin && ref_mz < prec_win_end) {
            selected_seed_list.push_back(seed_ptr_list[j]);
          }
        }
        SeedEnvPtr seed_ptr;
        if (selected_seed_list.size() > 0) {
          // choose the highest intensity one
          std::sort(selected_seed_list.begin(), selected_seed_list.end(),
                    SeedEnv::cmpSeedInteDec);
          seed_ptr = selected_seed_list[0];
          seed_ptr = seed_env_util::relaxProcessSeedEnvPtr(
              seed_ptr, matrix_ptr, score_para_ptr, sn_ratio);
        }
        if (seed_ptr == nullptr) {
          continue;
        }
        EnvCollPtr env_coll_ptr = env_coll_util::findEnvCollWithSingleEnv(
            matrix_ptr, seed_ptr, score_para_ptr, sn_ratio,
            topfd_para_ptr->getSplitIntensityRatio());
        if (env_coll_ptr == nullptr) continue;
        if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                                 env_coll_list, score_para_ptr,
                                                 sn_ratio, &rt_cache)) {
          env_coll_ptr->removePeakData(matrix_ptr);
          continue;
        }
        env_coll_ptr->refineMonoMass();
        ECScorePtr ecscore_ptr =
            std::make_shared<ECScore>(env_coll_ptr, matrix_ptr, sn_ratio);
        if (ecscore_ptr->getScore() < 0 ||
            std::isnan(ecscore_ptr->getScore())) {
          continue;
        }
        env_coll_ptr->setEcscore(ecscore_ptr->getScore());
        env_coll_ptr->removePeakData(matrix_ptr);
        env_coll_list.push_back(env_coll_ptr);
        ecscore_list.push_back(ecscore_ptr);
      }
    }
    std::cout << std::endl;
  }

  std::cout << "Sliding window: max resident rows " << max_live_rows << " of "
            << row_num << std::endl;

  FracFeaturePtrVec frac_features;
  for (std::size_t i = 0; i < env_coll_list.size(); i++) {
    EnvCollPtr env_coll_ptr = env_coll_list[i];
    FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(
        i, deconv_ms1_ptr_vec, score_para_ptr->frac_id_,
        score_para_ptr->file_name_, env_coll_ptr, matrix_ptr, sn_ratio);
    frac_features.push_back(frac_feat_ptr);
  }

  // map MS2 features
  SpecFeaturePtrVec ms2_features;
  env_coll_assign::assignEnvColls(frac_features, env_coll_list,
                                  ms2_header_ptr_2d, ms2_features,
                                  topfd_para_ptr->getMs1EcscoreCutoff());
  std::cout << "Number of proteoform features: " << env_coll_list.size()
            << std::endl;
  /// output files
  if (topfd_para_ptr->isOutputCsvFeatureFile()) {
    std::string feat_file_name = output_base_name + "_ms1.csv";
    ecscore_writer::writeScores(feat_file_name, ecscore_list);
    std::string batmass_file_name =
        output_base_name + "_" + "frac_ms1.mzrt.csv";
    int num_spec = ms1_raw_spec_num;
    frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features,
                                              num_spec);
  }

  std::string frac_feat_xml_file_name = output_base_name + "_feature.xml";
  frac_feature_writer::writeXmlFeatures(frac_feat_xml_file_name, frac_features);
  std::string frac_feat_file_name = output_base_name + "_" + "ms1.feature";
  frac_feature_writer::writeFeatures(frac_feat_file_name, frac_features);

  std::string ms2_feat_file_name = output_base_name + "_" + "ms2.feature";
  spec_feature_writer::writeFeatures(ms2_feat_file_name, ms2_features);
}


namespace {

// Batched EnvCNN scoring for the MS1 feature detection (issue #586).
//
// Scoring one candidate per EnvCNN call streams the model's 83 MB of weights
// once per candidate: 2.85 ms each, 246 s of a 610 s run on 480_F3_02.
// Scoring 64 at a time costs 0.162 ms per candidate. The loop cannot simply
// score in batches, because an accepted feature removes its peaks from the
// matrix that later candidates are built and scored from, so a candidate is
// parked with the region it reads and is rebuilt and rescored on its own when
// a removal made after it was parked touches that region. On 480_F3_02 that
// happens to 0.5% of the candidates (toppic_test experiment 62).
//
// Two candidates that checkExistingFeatures would merge have the same mass
// and overlapping retention times, so their regions intersect and the second
// one is rebuilt; that is why a parked candidate needs no second existing
// feature check.
constexpr size_t kEnvcnnBatchSize = 64;

// the m/z range an envelope set covers, widened by the 0.1 Dalton the EnvCNN
// input reaches for its noise peaks, and the scans it spans
struct ReadRegion {
  double mz_lo, mz_hi;
  int spec_lo, spec_hi;
};

inline bool regionsIntersect(const ReadRegion& a, const ReadRegion& b) {
  return a.mz_lo <= b.mz_hi && b.mz_lo <= a.mz_hi && a.spec_lo <= b.spec_hi &&
         b.spec_lo <= a.spec_hi;
}

std::vector<ReadRegion> getReadRegions(const EnvCollPtr& coll_ptr) {
  std::vector<ReadRegion> regions;
  for (auto& env_set_ptr : coll_ptr->getEnvSetList()) {
    if (env_set_ptr == nullptr) continue;
    EnvPeakPtrVec peaks = env_set_ptr->getSeedPtr()->getPeakPtrList();
    double lo = std::numeric_limits<double>::max();
    double hi = std::numeric_limits<double>::lowest();
    for (auto& peak_ptr : peaks) {
      if (peak_ptr == nullptr) continue;
      lo = std::min(lo, peak_ptr->getPosition());
      hi = std::max(hi, peak_ptr->getPosition());
    }
    if (hi < lo) continue;
    regions.push_back({lo - 0.1, hi + 0.1, env_set_ptr->getStartSpecId(),
                       env_set_ptr->getEndSpecId()});
  }
  return regions;
}

// a candidate waiting for its batch to be scored
struct ParkedCand {
  SeedEnvPtr orig_seed_ptr;            // to rebuild it if it goes stale
  EnvCollPtr coll_ptr;
  std::vector<float> tensor;           // the EnvCNN input as of parking
  std::vector<ReadRegion> regions;
  size_t removal_num;                  // removals made before it was parked
};

}  // namespace

void processMs1(const TopfdParaPtr& topfd_para_ptr) {
  if (topfd_para_ptr->isMissingLevelOne()) {
    return;
  }
  if (std::getenv("TOPFD_MS1_WINDOW") != nullptr) {
    processMs1Windowed(topfd_para_ptr);
    return;
  }

  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(
      topfd_para_ptr->getFracId(), topfd_para_ptr->getMzmlFileName(),
      topfd_para_ptr->getMaxCharge(), topfd_para_ptr->getMs1MinScanNum());
  score_para_ptr->seed_env_inte_corr_tole_cutoff_ =
      topfd_para_ptr->getSeedEnvInteCorrToleCutoff();
  // read deconvoluted MS1 peaks
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();
  std::string ms1_file_name = output_base_name + "_ms1.msalign";
  DeconvMsPtrVec deconv_ms1_ptr_vec;
  msalign_reader_util::readAllSpectra(ms1_file_name, deconv_ms1_ptr_vec);

  // read ms1 raw peaks and ms2_headers
  PeakPtrVec2D ms1_mzml_peaks;
  MsHeaderPtr2D ms2_header_ptr_2d;
  MzmlMsGroupReaderPtr mzml_reader_ptr = std::make_shared<MzmlMsGroupReader>(
      topfd_para_ptr->getMzmlFileName(), topfd_para_ptr->getPrecWindowWidth(),
      topfd_para_ptr->getActivation(), topfd_para_ptr->getFracId(),
      topfd_para_ptr->isFaims(), topfd_para_ptr->getFaimsVoltage(),
      topfd_para_ptr->isMissingLevelOne());
  mzml_reader_ptr->getMs1Map(ms1_mzml_peaks, ms2_header_ptr_2d);
  mzml_reader_ptr = nullptr;

  // Prepare seed envelopes
  SeedEnvPtrVec seed_ptrs;
  SeedEnvPtr2D seed_ptr_2d;
  for (auto& ms1_data : deconv_ms1_ptr_vec) {
    DeconvPeakPtrVec peaks = ms1_data->getPeakPtrVec();
    SeedEnvPtrVec one_spec_seed_ptrs;
    for (auto& peak : peaks) {
      // SeedEnv objects are immutable after construction; the global
      // intensity-sorted list and the per-spectrum list share them.
      SeedEnvPtr seed_ptr = std::make_shared<SeedEnv>(peak);
      seed_ptrs.push_back(seed_ptr);
      one_spec_seed_ptrs.push_back(seed_ptr);
    }
    seed_ptr_2d.push_back(one_spec_seed_ptrs);
  }

  std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnv::cmpSeedInteDec);
  // write_out_files::write_seed_envelopes(seed_envs, "envs.csv");

  double sn_ratio = topfd_para_ptr->getMsOneSnRatio();
  bool single_scan_noise = topfd_para_ptr->isUseSingleScanNoiseLevel();
  MsHeaderPtrVec ms1_header_ptr_vec;
  for (auto& ms1_spectrum : deconv_ms1_ptr_vec) {
    ms1_header_ptr_vec.push_back(ms1_spectrum->getMsHeaderPtr());
  }
  /// Prepare data -- Peak Matrix
  MsMapPtr matrix_ptr = std::make_shared<MsMap>(
      ms1_mzml_peaks, ms1_header_ptr_vec, score_para_ptr->bin_size_, sn_ratio,
      single_scan_noise);
  // The raw peaks were copied into the matrix; only the spectrum count is
  // needed below. Release the whole-run raw peak buffer.
  size_t ms1_raw_spec_num = ms1_mzml_peaks.size();
  PeakPtrVec2D().swap(ms1_mzml_peaks);

  if (score_para_ptr->min_scan_num_ >= 2) {
    matrix_ptr->removeNonNeighbors(score_para_ptr->neighbor_mz_tole_);
  }

  /// Extract Features
  LOG_DEBUG("Number of seed envelopes: " << seed_ptrs.size());
  int seed_num = seed_ptrs.size();
  EnvCollPtrVec env_coll_list;
  ECScorePtrVec ecscore_list;
  env_coll_util::EnvCollRtCache rt_cache;
  // the candidates waiting to be scored, and the peak removals made since
  // the oldest of them was parked (issue #586)
  std::vector<ParkedCand> parked;
  std::vector<ReadRegion> removals;
  long batch_call_num = 0, cand_num = 0, rebuilt_num = 0;

  auto removePeaks = [&](const EnvCollPtr& coll_ptr) {
    std::vector<ReadRegion> regions = getReadRegions(coll_ptr);
    removals.insert(removals.end(), regions.begin(), regions.end());
    coll_ptr->removePeakData(matrix_ptr);
  };

  // score the parked candidates in one call and apply them in parking order
  auto flushBatch = [&]() {
    if (parked.empty()) return;
    std::vector<float> tensor;
    tensor.reserve(parked.size() * parked[0].tensor.size());
    for (auto& cand : parked)
      tensor.insert(tensor.end(), cand.tensor.begin(), cand.tensor.end());
    std::vector<double> batch_scores =
        onnx_env_cnn::predict(static_cast<int>(parked.size()), tensor);
    batch_call_num++;
    for (size_t i = 0; i < parked.size(); i++) {
      ParkedCand& cand = parked[i];
      bool stale = false;
      for (size_t r = cand.removal_num; r < removals.size() && !stale; r++)
        for (const auto& region : cand.regions)
          if (regionsIntersect(removals[r], region)) {
            stale = true;
            break;
          }
      EnvCollPtr env_coll_ptr = cand.coll_ptr;
      double envcnn_score = batch_scores[i];
      if (stale) {
        // a removal touched what this candidate reads: redo it against the
        // matrix as it stands now, exactly as the sequential loop would
        rebuilt_num++;
        SeedEnvPtr seed_ptr = seed_env_util::preprocessSeedEnvPtr(
            cand.orig_seed_ptr, matrix_ptr, score_para_ptr, sn_ratio);
        if (seed_ptr == nullptr) continue;
        env_coll_ptr = env_coll_util::findEnvColl(
            matrix_ptr, seed_ptr, score_para_ptr, sn_ratio,
            topfd_para_ptr->getSplitIntensityRatio());
        if (env_coll_ptr == nullptr) continue;
        if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                                 env_coll_list, score_para_ptr,
                                                 sn_ratio, &rt_cache)) {
          removePeaks(env_coll_ptr);
          continue;
        }
        env_coll_ptr->refineMonoMass();
        envcnn_score =
            comp_env_cnn_score::compEnvcnnScore(matrix_ptr, env_coll_ptr);
      }
      ECScorePtr ecscore_ptr = std::make_shared<ECScore>(
          env_coll_ptr, matrix_ptr, sn_ratio, envcnn_score);
      if (ecscore_ptr->getScore() < topfd_para_ptr->getMs1EcscoreCutoff()) {
        continue;
      }
      env_coll_ptr->setEcscore(ecscore_ptr->getScore());
      removePeaks(env_coll_ptr);
      env_coll_list.push_back(env_coll_ptr);
      ecscore_list.push_back(ecscore_ptr);
    }
    parked.clear();
    removals.clear();
  };

  for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
    int count = seed_env_idx + 1;
    if (count % 100 == 0 || count == seed_num) {
      double perc = static_cast<int>(count * 100 / seed_num);
      std::cout << "\r" << "Processing feature " << count << " ...       "
                << perc << "\% finished." << std::flush;
    }
    SeedEnvPtr orig_seed_ptr = seed_ptrs[seed_env_idx];
    SeedEnvPtr seed_ptr = seed_env_util::preprocessSeedEnvPtr(
        orig_seed_ptr, matrix_ptr, score_para_ptr, sn_ratio);
    if (seed_ptr == nullptr) continue;
    EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(
        matrix_ptr, seed_ptr, score_para_ptr, sn_ratio,
        topfd_para_ptr->getSplitIntensityRatio());
    if (env_coll_ptr == nullptr) continue;
    if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                             env_coll_list, score_para_ptr,
                                             sn_ratio, &rt_cache)) {
      removePeaks(env_coll_ptr);
      continue;
    }
    env_coll_ptr->refineMonoMass();
    cand_num++;
    parked.push_back({orig_seed_ptr, env_coll_ptr,
                      comp_env_cnn_score::getEnvcnnTensor(matrix_ptr,
                                                          env_coll_ptr),
                      getReadRegions(env_coll_ptr), removals.size()});
    if (parked.size() >= kEnvcnnBatchSize) flushBatch();
  }
  flushBatch();
  std::cout << std::endl;
  std::cout << "EnvCNN scoring: " << cand_num << " candidates in "
            << batch_call_num << " batched calls, " << rebuilt_num
            << " rebuilt" << std::endl;

  if (topfd_para_ptr->isSearchPrecWindow()) {
    // add ms1 feature based on precursor windows
    // set min match envelope to 1 to accept single scan features
    score_para_ptr->min_scan_num_ = 1;
    score_para_ptr->min_match_peak_ = 1;
    sn_ratio = 0;
    matrix_ptr->reconstruct(sn_ratio, single_scan_noise);
    rt_cache.clear();
    int ms1_spec_num = deconv_ms1_ptr_vec.size();
    for (std::size_t ms1_idx = 0; ms1_idx < deconv_ms1_ptr_vec.size();
         ms1_idx++) {
      double perc = static_cast<int>((ms1_idx + 1) * 100 / ms1_spec_num);
      int scan =
          deconv_ms1_ptr_vec[ms1_idx]->getMsHeaderPtr()->getFirstScanNum();
      std::cout << "\r" << "Additional feature search MS1 spectrum scan "
                << scan << " ...       " << perc << "\% finished."
                << std::flush;
      for (std::size_t i = 0; i < ms2_header_ptr_2d[ms1_idx].size(); i++) {
        MsHeaderPtr header_ptr = ms2_header_ptr_2d[ms1_idx][i];
        if (env_coll_assign::checkEnvColl(header_ptr, env_coll_list)) {
          continue;
        }
        double prec_win_begin = header_ptr->getPrecWinBegin();
        double prec_win_end = header_ptr->getPrecWinEnd();
        SeedEnvPtrVec seed_ptr_list = seed_ptr_2d[ms1_idx];
        SeedEnvPtrVec selected_seed_list;
        LOG_DEBUG("ms1 id  " << ms1_idx << " seed number "
                             << seed_ptr_list.size());
        for (std::size_t i = 0; i < seed_ptr_list.size(); i++) {
          double ref_mz = seed_ptr_list[i]->getReferMz();
          if (ref_mz > prec_win_begin && ref_mz < prec_win_end) {
            selected_seed_list.push_back(seed_ptr_list[i]);
          }
        }
        SeedEnvPtr seed_ptr;
        if (selected_seed_list.size() > 0) {
          // choose the highest intensity one
          std::sort(selected_seed_list.begin(), selected_seed_list.end(),
                    SeedEnv::cmpSeedInteDec);
          seed_ptr = selected_seed_list[0];
          seed_ptr = seed_env_util::relaxProcessSeedEnvPtr(
              seed_ptr, matrix_ptr, score_para_ptr, sn_ratio);
        }
        if (seed_ptr == nullptr) {
          continue;
        }

        EnvCollPtr env_coll_ptr = env_coll_util::findEnvCollWithSingleEnv(
            matrix_ptr, seed_ptr, score_para_ptr, sn_ratio,
            topfd_para_ptr->getSplitIntensityRatio());
        if (env_coll_ptr == nullptr) continue;
        if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                                 env_coll_list, score_para_ptr,
                                                 sn_ratio, &rt_cache)) {
          env_coll_ptr->removePeakData(matrix_ptr);
          continue;
        }
        env_coll_ptr->refineMonoMass();
        ECScorePtr ecscore_ptr =
            std::make_shared<ECScore>(env_coll_ptr, matrix_ptr, sn_ratio);
        if (ecscore_ptr->getScore() < 0 ||
            std::isnan(ecscore_ptr->getScore())) {
          continue;
        }
        env_coll_ptr->setEcscore(ecscore_ptr->getScore());
        env_coll_ptr->removePeakData(matrix_ptr);
        env_coll_list.push_back(env_coll_ptr);
        ecscore_list.push_back(ecscore_ptr);
      }
    }
    std::cout << std::endl;
  }

  FracFeaturePtrVec frac_features;
  for (std::size_t i = 0; i < env_coll_list.size(); i++) {
    EnvCollPtr env_coll_ptr = env_coll_list[i];
    FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(
        i, deconv_ms1_ptr_vec, score_para_ptr->frac_id_,
        score_para_ptr->file_name_, env_coll_ptr, matrix_ptr, sn_ratio);
    frac_features.push_back(frac_feat_ptr);
  }

  // map MS2 features
  SpecFeaturePtrVec ms2_features;
  env_coll_assign::assignEnvColls(frac_features, env_coll_list,
                                  ms2_header_ptr_2d, ms2_features,
                                  topfd_para_ptr->getMs1EcscoreCutoff());
  std::cout << "Number of proteoform features: " << env_coll_list.size()
            << std::endl;
  /// output files
  if (topfd_para_ptr->isOutputCsvFeatureFile()) {
    std::string feat_file_name = output_base_name + "_ms1.csv";
    ecscore_writer::writeScores(feat_file_name, ecscore_list);
    std::string batmass_file_name =
        output_base_name + "_" + "frac_ms1.mzrt.csv";
    int num_spec = ms1_raw_spec_num;
    frac_feature_writer::writeBatMassFeatures(batmass_file_name, frac_features,
                                              num_spec);
  }

  std::string frac_feat_xml_file_name = output_base_name + "_feature.xml";
  frac_feature_writer::writeXmlFeatures(frac_feat_xml_file_name, frac_features);
  std::string frac_feat_file_name = output_base_name + "_" + "ms1.feature";
  frac_feature_writer::writeFeatures(frac_feat_file_name, frac_features);

  // SampleFeaturePtrVec sample_features;
  // feature_util::getSampleFeatures(sample_features, frac_features,
  // ms2_features); std::string sample_feature_file_name = output_base_name +
  // "_"  + "ms1.feature";
  // sample_feature_writer::writeFeatures(sample_feature_file_name,
  // sample_features);
  std::string ms2_feat_file_name = output_base_name + "_" + "ms2.feature";
  spec_feature_writer::writeFeatures(ms2_feat_file_name, ms2_features);
}

void processMs2(const TopfdParaPtr& topfd_para_ptr) {
  EcscoreParaPtr score_para_ptr = std::make_shared<EcscorePara>(
      topfd_para_ptr->getFracId(), topfd_para_ptr->getMzmlFileName(),
      topfd_para_ptr->getMaxCharge(), topfd_para_ptr->getMs2MinScanNum());
  score_para_ptr->seed_env_inte_corr_tole_cutoff_ =
      topfd_para_ptr->getSeedEnvInteCorrToleCutoff();

  // read deconvoluted MS2 peaks
  std::string output_base_name = topfd_para_ptr->getOutputBaseName();
  std::string ms2_file_name = output_base_name + "_ms2_raw.msalign";
  DeconvMsPtrVec deconv_ms2_ptr_vec;
  msalign_reader_util::readAllSpectra(ms2_file_name, deconv_ms2_ptr_vec);

  // get isolation window
  std::set<std::pair<double, double>> win_set;
  for (auto& ms2_data : deconv_ms2_ptr_vec) {
    if (ms2_data->getMsHeaderPtr()->getMsLevel() == 1) continue;
    std::pair<double, double> cur_win(
        ms2_data->getMsHeaderPtr()->getPrecWinBegin(),
        ms2_data->getMsHeaderPtr()->getPrecWinEnd());
    win_set.insert(cur_win);
  }
  std::vector<std::pair<double, double>> win_list(win_set.begin(),
                                                  win_set.end());

  for (std::size_t i = 0; i < win_list.size(); i++) {
    std::pair<double, double> cur_win = win_list[i];
    std::cout << std::fixed << std::setprecision(6)
              << "Processing isolation window: " << (i + 1) << " of "
              << win_list.size() << " [" << cur_win.first << ","
              << cur_win.second << "]" << std::endl;
    std::cout << std::setprecision(2);
    // read ms1 raw peaks and ms2_headers
    PeakPtrVec2D ms2_mzml_peaks;

    MzmlMsGroupReaderPtr mzml_reader_ptr = std::make_shared<MzmlMsGroupReader>(
        topfd_para_ptr->getMzmlFileName(), topfd_para_ptr->getPrecWindowWidth(),
        topfd_para_ptr->getActivation(), topfd_para_ptr->getFracId(),
        topfd_para_ptr->isFaims(), topfd_para_ptr->getFaimsVoltage(),
        topfd_para_ptr->isMissingLevelOne());
    int peak_num = mzml_reader_ptr->getMs2Map(ms2_mzml_peaks, cur_win.first);
    mzml_reader_ptr = nullptr;

    FracFeaturePtrVec frac_features;
    ECScorePtrVec ecscore_list;
    // Get envelope collections and features
    if (peak_num == 0) {
      std::cout << "No MS2 peaks in this isolation window." << std::endl;
    } else {
      // Prepare seed envelopes
      int spec_id = 0;
      DeconvMsPtrVec deconv_ms2_ptr_shortlisted_vec;
      SeedEnvPtrVec seed_ptrs;
      for (auto& ms2_data : deconv_ms2_ptr_vec) {
        double ms_win_begin = ms2_data->getMsHeaderPtr()->getPrecWinBegin();
        if (ms_win_begin != cur_win.first) continue;
        deconv_ms2_ptr_shortlisted_vec.push_back(ms2_data);
        DeconvPeakPtrVec peaks = ms2_data->getPeakPtrVec();
        for (auto& peak : peaks) {
          peak->setSpId(spec_id);
          SeedEnvPtr seed_ptr_1 = std::make_shared<SeedEnv>(peak);
          seed_ptrs.push_back(seed_ptr_1);
        }
        spec_id = spec_id + 1;
      }

      std::sort(seed_ptrs.begin(), seed_ptrs.end(), SeedEnv::cmpSeedInteDec);

      double sn_ratio = topfd_para_ptr->getMsTwoSnRatio();
      bool single_scan_noise = topfd_para_ptr->isUseSingleScanNoiseLevel();
      MsHeaderPtrVec ms2_header_ptr_vec;
      for (auto& ms2_spectrum : deconv_ms2_ptr_shortlisted_vec) {
        ms2_header_ptr_vec.push_back(ms2_spectrum->getMsHeaderPtr());
      }
      /// Prepare data -- Peak Matrix
      MsMapPtr matrix_ptr = std::make_shared<MsMap>(
          ms2_mzml_peaks, ms2_header_ptr_vec, score_para_ptr->bin_size_,
          topfd_para_ptr->getMsTwoSnRatio(), single_scan_noise);

      if (score_para_ptr->min_scan_num_ > 2) {
        matrix_ptr->removeNonNeighbors(score_para_ptr->neighbor_mz_tole_);
      }

      /// Extract Features
      LOG_DEBUG("Number of seed envelopes: " << seed_ptrs.size());
      int seed_num = static_cast<int>(seed_ptrs.size());
      EnvCollPtrVec env_coll_list;
      env_coll_util::EnvCollRtCache rt_cache;
      for (int seed_env_idx = 0; seed_env_idx < seed_num; seed_env_idx++) {
        int count = seed_env_idx + 1;
        if (count % 100 == 0 || count == seed_num) {
          double perc = static_cast<int>(count * 100 / seed_num);
          std::cout << "\r" << "Processing feature " << count << " ...       "
                    << perc << "\% finished." << std::flush;
        }
        SeedEnvPtr seed_ptr = seed_ptrs[seed_env_idx];
        seed_ptr = seed_env_util::preprocessSeedEnvPtr(
            seed_ptr, matrix_ptr, score_para_ptr, sn_ratio);
        if (seed_ptr == nullptr) continue;
        EnvCollPtr env_coll_ptr = env_coll_util::findEnvColl(
            matrix_ptr, seed_ptr, score_para_ptr, sn_ratio,
            topfd_para_ptr->getSplitIntensityRatio());
        if (env_coll_ptr == nullptr) continue;
        if (env_coll_util::checkExistingFeatures(matrix_ptr, env_coll_ptr,
                                                 env_coll_list, score_para_ptr,
                                                 sn_ratio, &rt_cache)) {
          env_coll_ptr->removePeakData(matrix_ptr);
          continue;
        }
        env_coll_ptr->refineMonoMass();
        ECScorePtr ecscore_ptr = std::make_shared<ECScore>(
            env_coll_ptr, matrix_ptr, topfd_para_ptr->getMsTwoSnRatio());
        if (ecscore_ptr->getScore() < topfd_para_ptr->getMs2EcscoreCutoff()) {
          continue;
        }
        ecscore_list.push_back(ecscore_ptr);
        env_coll_ptr->setEcscore(ecscore_ptr->getScore());
        env_coll_ptr->removePeakData(matrix_ptr);
        env_coll_list.push_back(env_coll_ptr);
      }
      std::cout << std::endl;

      std::cout << "Number of fragment features: " << env_coll_list.size()
                << std::endl;

      for (std::size_t i = 0; i < env_coll_list.size(); i++) {
        EnvCollPtr env_coll_ptr = env_coll_list[i];
        FracFeaturePtr frac_feat_ptr = env_coll_util::getFracFeature(
            i, deconv_ms2_ptr_shortlisted_vec, score_para_ptr->frac_id_,
            score_para_ptr->file_name_, env_coll_ptr, matrix_ptr, sn_ratio);
        frac_features.push_back(frac_feat_ptr);
      }
    }

    std::string feat_file_name =
        output_base_name + "_" + std::to_string(cur_win.first) + "_ms2.csv";
    ecscore_writer::writeScores(feat_file_name, ecscore_list);
    /// output files
    if (topfd_para_ptr->isOutputCsvFeatureFile()) {
      std::string batmass_file_name;
      batmass_file_name = output_base_name + "_" +
                          std::to_string(cur_win.first) + "_frac_ms2.mzrt.csv";
      int num_spec = ms2_mzml_peaks.size();
      frac_feature_writer::writeBatMassFeatures(batmass_file_name,
                                                frac_features, num_spec);
    }
  }
}

}  // namespace env_coll_detect

}  // namespace toppic
