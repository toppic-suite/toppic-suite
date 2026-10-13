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

#ifndef TOPPIC_TOPFD_COMMON_TOPFD_PARA_HPP_
#define TOPPIC_TOPFD_COMMON_TOPFD_PARA_HPP_

#include <sqlite3.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace toppic {

class TopfdPara;

using TopfdParaPtr = std::shared_ptr<TopfdPara>;

class TopfdPara {
 public:
  TopfdPara() = default;

  ~TopfdPara();

  // msalign_ext ("ms1.msalign", "ms2.msalign", ...) names the msalign file
  // the block is written into: it adds an "Msalign file name:" line
  // (<output base name>_<ext>, so FAIMS runs show their per-voltage file).
  // Empty for the console printout, which belongs to no file.
  std::string getTopfdParaStr(const std::string& prefix,
                              const std::string& sep,
                              const std::string& msalign_ext = "") const;

  std::string getParaStr(const std::string& prefix, const std::string& sep,
                         const std::string& msalign_ext = "") const;

  void createSqlDb(const std::string& db_name);

  // Compacts (VACUUM) and closes the SQLite database once it is complete;
  // nothing when none is open. TopFD and TopDIA call it at the end of each
  // file or FAIMS fraction (not topfd -T, whose database holds one
  // spectrum: createSqlDb and the destructor close a database still open
  // without compacting it).
  void closeSqlDb();

  const std::string& getExeDir() const { return exe_dir_; }
  const std::string& getResourceDir() const { return resource_dir_; }
  bool isMissingLevelOne() const { return missing_level_one_; }
  int getMaxCharge() const { return max_charge_; }

  double getMassGroupBoundary() const { return mass_group_boundary_; }
  double getMaxPeakSharingRatio() const { return max_peak_sharing_ratio_; }
  double getMaxMass() const { return max_mass_; }
  double getMzError() const { return mz_error_; }
  int getMaxMissPeakNum() const { return max_miss_peak_num_; }
  double getMsOneRankPeakDistance() const { return ms1_rank_peak_distance_; }
  double getMsTwoRankPeakDistance() const { return ms2_rank_peak_distance_; }
  bool isMsTwoHarmonicCheck() const { return ms2_harmonic_check_; }
  bool isRunFilterByMz() const { return run_filter_by_mz_; }

  bool isEstimateMinInte() const { return estimate_min_inte_; }
  double getMsOneSnRatio() const { return ms_one_sn_ratio_; }
  double getMsTwoSnRatio() const { return ms_two_sn_ratio_; }
  double getSplitIntensityRatio() const { return split_intensity_ratio_; }
  double getPrecWindowWidth() const { return prec_window_; }
  // the score that ranks the envelopes (--env-rank of topfd and topdia):
  // "envcnn" or "msdeconv"
  static const std::vector<std::string>& getDeconvScores();
  const std::string& getDeconvScore() const { return deconv_score_; }
  bool isSortUseMsDeconv() const { return deconv_score_ == "msdeconv"; }
  bool isAANumBasedFilter() const { return aa_num_based_filter_; }
  double getMs1EnvCnnScoreCutoff() const { return ms1_env_cnn_score_cutoff_; }
  double getMs2EnvCnnScoreCutoff() const { return ms2_env_cnn_score_cutoff_; }
  const std::string& getActivation() const { return activation_; }
  bool isGeneSql() const { return gene_sql_; }
  bool isKeepUnusedPeaks() const { return keep_unused_peaks_; }
  bool isOutputAvgMass() const { return output_avg_mass_; }
  bool isOutputCsvFeatureFile() const { return output_csv_feature_file_; }
  bool isOutputDpEnvs() const { return output_dp_envs_; }
  // the deconvolution methods, named as the --deconv-ms1 and --deconv-ms2
  // options of topfd and topdia name them
  static const std::vector<std::string>& getMs1DeconvMethods();
  static const std::vector<std::string>& getMs2DeconvMethods();
  const std::string& getMs1DeconvMethod() const { return ms1_deconv_method_; }
  const std::string& getMs2DeconvMethod() const { return ms2_deconv_method_; }
  bool isMs1GreedyDeconv() const { return ms1_deconv_method_ == "greedy"; }
  bool isMs2GreedyDeconv() const { return ms2_deconv_method_ == "greedy"; }
  // the greedy searches' limit of the peaks a new envelope may share with
  // each envelope of its charge already reported (-1: no limit)
  int getMs1MaxSharedPeakNum() const { return ms1_max_shared_peak_num_; }
  int getMs2MaxSharedPeakNum() const { return ms2_max_shared_peak_num_; }
  bool isUseFlashDeconv() const { return ms2_deconv_method_ == "flashdeconv"; }
  bool isMs2PrecPairRefine() const { return ms2_prec_pair_refine_; }

  // topfd --ptcr: every MS/MS spectrum is a PTCR spectrum (the activation
  // is then PTCR; see ptcr_ms2_deconv)
  bool isPtcrForced() const { return ptcr_forced_; }
  // the MS3 spectral deconvolution (DeconvMs3Process), on by default;
  // hidden --disable-ms3 turns it off
  bool isMs3Deconv() const { return ms3_deconv_; }
  int getMs3ScanNum() const { return ms_3_scan_num_; }
  int getThreadNum() const { return thread_num_; }
  double getMs1EcscoreCutoff() const { return ms1_ecscore_cutoff_; }
  double getMs2EcscoreCutoff() const { return ms2_ecscore_cutoff_; }
  bool isSearchPrecWindow() const { return search_prec_window_; }
  // True when the input mzML carries MS/MS precursor windows (so prec_window_,
  // the default width, is unused). Set from MzmlProfile::hasPrecWindow().
  bool isFilePrecWindow() const { return file_prec_window_; }
  bool isUseSingleScanNoiseLevel() const {
    return use_single_scan_noise_level_;
  }
  bool isTextPeakList() const { return text_peak_list_; }
  double getPrecInteCutoffRatio() const { return prec_inte_cutoff_ratio_; }

  const std::string& getMzmlFileName() const { return mzml_file_name_; }
  const std::string& getOutputBaseName() const { return output_base_name_; }

  int getFracId() const { return frac_id_; }
  bool isFaims() const { return is_faims_; }
  double getFaimsVoltage() const { return faims_volt_; }
  int getMs1ScanNum() const { return ms_1_scan_num_; }
  int getMs2ScanNum() const { return ms_2_scan_num_; }
  int getMs1MinScanNum() const { return ms1_min_scan_num_; }
  double getSeedEnvInteCorrToleCutoff() const {
    return seed_env_inte_corr_tole_cutoff_;
  }
  int getMs2MinScanNum() const { return ms2_min_scan_num_; }
  sqlite3* getSqlDb() const { return sql_db_; }

  void setExeDir(std::string dir) { exe_dir_ = std::move(dir); }
  void setResourceDir(std::string dir) { resource_dir_ = std::move(dir); }
  void setMissingLevelOne(bool missing) { missing_level_one_ = missing; }
  void setMaxCharge(int charge) { max_charge_ = charge; }

  void setMassGroupBoundary(double mass) {
    mass_group_boundary_ = mass;
  }
  void setMaxPeakSharingRatio(double ratio) {
    max_peak_sharing_ratio_ = ratio;
  }
  void setMs1MaxSharedPeakNum(int num) { ms1_max_shared_peak_num_ = num; }
  void setMs2MaxSharedPeakNum(int num) { ms2_max_shared_peak_num_ = num; }
  void setMaxMass(double mass) { max_mass_ = mass; }
  void setMzError(double error) { mz_error_ = error; }
  void setMsOneSnRatio(double ratio) { ms_one_sn_ratio_ = ratio; }
  void setMsTwoSnRatio(double ratio) { ms_two_sn_ratio_ = ratio; }
  void setSplitIntensityRatio(double ratio) { split_intensity_ratio_ = ratio; }
  void setPrecWindowWidth(double window) { prec_window_ = window; }
  // false, and the score unchanged, for a name not in getDeconvScores()
  bool setDeconvScore(const std::string& score);
  void setAANumBasedFilter(bool filter) { aa_num_based_filter_ = filter; }
  void setMs1EnvCnnScoreCutoff(double cutoff) {
    ms1_env_cnn_score_cutoff_ = cutoff;
  }
  void setMs2EnvCnnScoreCutoff(double cutoff) {
    ms2_env_cnn_score_cutoff_ = cutoff;
  }
  void setActivation(std::string activation) {
    activation_ = std::move(activation);
  }
  void setKeepUnusedPeaks(bool keep) { keep_unused_peaks_ = keep; }
  void setOutputAvgMass(bool output) { output_avg_mass_ = output; }
  void setOutputCsvFeatureFile(bool output) {
    output_csv_feature_file_ = output;
  }
  void setGeneSql(bool gene_sql) { gene_sql_ = gene_sql; }
  void setOutputDpEnvs(bool output) { output_dp_envs_ = output; }
  // false, and the method unchanged, for a name not in getMs1DeconvMethods()
  // (getMs2DeconvMethods())
  bool setMs1DeconvMethod(const std::string& method);
  bool setMs2DeconvMethod(const std::string& method);
  void setMs2PrecPairRefine(bool refine) { ms2_prec_pair_refine_ = refine; }

  void setPtcrForced(bool ptcr_forced) { ptcr_forced_ = ptcr_forced; }
  void setMs3Deconv(bool ms3_deconv) { ms3_deconv_ = ms3_deconv; }
  void setMs3ScanNumber(int ms3_scan_num) { ms_3_scan_num_ = ms3_scan_num; }
  void setMaxMissPeakNum(int num) { max_miss_peak_num_ = num; }
  void setMsOneRankPeakDistance(double d) { ms1_rank_peak_distance_ = d; }
  void setMsTwoRankPeakDistance(double d) { ms2_rank_peak_distance_ = d; }
  void setMsTwoHarmonicCheck(bool check) { ms2_harmonic_check_ = check; }
  void setRunFilterByMz(bool run) { run_filter_by_mz_ = run; }
  void setThreadNum(int num) { thread_num_ = num; }
  void setSearchPrecWindow(bool search) { search_prec_window_ = search; }
  void setFilePrecWindow(bool file_prec_window) {
    file_prec_window_ = file_prec_window;
  }
  void setUseSingleScanNoiseLevel(bool single_scan_noise) {
    use_single_scan_noise_level_ = single_scan_noise;
  }
  void setMs1EcscoreCutoff(double cutoff) { ms1_ecscore_cutoff_ = cutoff; }
  void setMs2EcscoreCutoff(double cutoff) { ms2_ecscore_cutoff_ = cutoff; }
  void setMs1MinScanNum(int min_scan_num) { ms1_min_scan_num_ = min_scan_num; }
  void setSeedEnvInteCorrToleCutoff(double cutoff) {
    seed_env_inte_corr_tole_cutoff_ = cutoff;
  }
  void setMs2MinScanNum(int min_scan_num) { ms2_min_scan_num_ = min_scan_num; }

  void setFracId(int frac_id) { frac_id_ = frac_id; }
  void setMzmlFileNameAndFaims(const std::string& mzml_file_name, bool is_faims,
                               double voltage);
  void setMs1ScanNumber(int ms1_scan_num) { ms_1_scan_num_ = ms1_scan_num; }
  void setMs2ScanNumber(int ms2_scan_num) { ms_2_scan_num_ = ms2_scan_num; }

  void setTextPeakList(bool text_peak_list) {
    text_peak_list_ = text_peak_list;
  }
  void setOutputMatchEnv(bool output_match_env) {
    output_match_env_ = output_match_env;
  }

 private:
  // Width of the label column in the parameter printout: every value starts at
  // the same column so the report is aligned across all sections. It is wide
  // enough for the longest label ("Filtering fragments using estimated fragment
  // number:").
  static constexpr int para_label_width_ = 53;

  // Total width of the "### <title> ###" section banners.
  static constexpr int para_banner_width_ = 55;

  // A banner line with the title centered and padded with '#' to a fixed width,
  // e.g. "############### Parameters ###############".
  static std::string banner(const std::string& prefix,
                            const std::string& title);

  std::string exe_dir_;
  std::string resource_dir_;

  // parameters for deconcovolution
  int max_charge_ = 30;

  // hidden option --mass-group-boundary: EnvPara::mass_group_boundary_[2],
  // above which envelope validation needs three matched and three
  // consecutive isotopic peaks instead of two
  // keep in step with EnvPara::DEFAULT_MASS_GROUP_BOUNDARY; this header
  // deliberately has no project includes, so the value is repeated
  double mass_group_boundary_ = 4250;
  // hidden option --max-peak-sharing-ratio: DpPara::max_peak_sharing_ratio_,
  // above which two envelopes sharing peaks do not coexist in the DP
  double max_peak_sharing_ratio_ = 1.0;
  double max_mass_ = 50000;
  // precursor window is used only when the mzML file does not
  // contain the precursor window information
  double prec_window_ = 3.0;
  bool missing_level_one_ = false;
  double mz_error_ = 0.02;
  // Maximum number of missing peaks allowed in a matched envelope; used to
  // initialize EnvPara::max_miss_peak_num_.
  int max_miss_peak_num_ = 1;
  // When false, the "filtering by mz" step is skipped; used to initialize
  // EnvPara::run_filter_by_mz_.
  bool run_filter_by_mz_ = true;
  // The window of the filtering-by-mz step, in Daltons, one value per
  // deconvolution level; used to set EnvPara::rank_peak_distance_, which is
  // shared by the two levels. Keep in step with
  // EnvPara::DEFAULT_MS1_RANK_PEAK_DISTANCE and
  // DEFAULT_MS2_RANK_PEAK_DISTANCE; this header deliberately has no project
  // includes, so the values are repeated.
  double ms1_rank_peak_distance_ = 12;
  double ms2_rank_peak_distance_ = 2;
  // When true, an MS/MS envelope selected by the dynamic programming is
  // replaced by a better scoring harmonic parent. On by default.
  bool ms2_harmonic_check_ = true;
  double ms_one_sn_ratio_ = 3.0;
  double ms_two_sn_ratio_ = 1.0;
  int thread_num_ = 1;
  std::string activation_ = "FILE";
  // keep unused peaks in dynamic programming
  bool keep_unused_peaks_ = false;
  // the score that ranks the envelopes (--env-rank of topfd and topdia):
  // "envcnn" (the EnvCNN score, the default) or "msdeconv"
  // (the MS-Deconv score)
  std::string deconv_score_ = "envcnn";
  // the EnvCNN score cutoffs of the MS1 and MS/MS envelopes (topfd's
  // --envcnn-cutoff-ms1 and --envcnn-cutoff-ms2) of the greedy and MS-Deconv
  // deconvolutions; not the envelopes added for precursor windows, not
  // FLASHDeconv. MS1 0.2: the greedy MS1 search's own minimum until then
  // (its envelopes below 0.2 are a third of them, seed the features that the
  // feature detection would otherwise spend 50% more time on, and gave no
  // identification; toppic_test experiment 104, run g6)
  double ms1_env_cnn_score_cutoff_ = 0.2;
  double ms2_env_cnn_score_cutoff_ = 0.0;
  bool aa_num_based_filter_ = false;
  bool output_csv_feature_file_ = false;
  // the database holds the MS1 peaks with their 3D levels (toppic_test
  // experiment 169)
  bool gene_sql_ = true;
  // set per input file from MzmlProfile::hasPrecWindow()
  bool file_prec_window_ = false;

  // parameters for feature identification
  double split_intensity_ratio_ = 2.5;
  bool search_prec_window_ = true;
  bool use_single_scan_noise_level_ = false;
  double ms1_ecscore_cutoff_ = 0.1;
  double ms2_ecscore_cutoff_ = 0;
  int ms1_min_scan_num_ = 1;
  // Minimum correlation between a seed envelope's intensities and the
  // theoretical ones in feature detection (EcscorePara). TopDIA sets it per
  // stage from its -y (MS1) / -Y (MS2) options; TopFD uses the default.
  double seed_env_inte_corr_tole_cutoff_ = 0.5;
  int ms2_min_scan_num_ = 1;

  // For an MS/MS spectrum, the precursor is not reported if its intensity
  // is less than the cutoff ratio * the intensity of the first precusor
  double prec_inte_cutoff_ratio_ = 0.1;

  //** Fixed parameter setting **
  // estimate min intensity using the method in Thrash.
  bool estimate_min_inte_ = true;
  // report average masses (not monoisotopic masses) of precursor and
  // fragment ions in the MS/MS msalign file
  bool output_avg_mass_ = false;
  bool output_match_env_ = false;
  // When true, dump the windowed candidate envelopes and the DP-selected
  // envelopes of every deconvoluted spectrum to win_envs.txt / dp_envs.txt
  // (debugging aid).
  bool output_dp_envs_ = false;
  // the MS1 deconvolution method (option --deconv-ms1): "greedy"
  // (greedy_ms1_deconv, deconv_ms1_process::deconvMsOneGreedy; TopFD's
  // default) or "msdeconv" (MS-Deconv's precursor-window step +
  // DeconvSingleSp; TopDIA's default, which its argument parser sets)
  std::string ms1_deconv_method_ = "greedy";
  // the MS/MS deconvolution method (option --deconv-ms2): "greedy"
  // (greedy_ms2_deconv, in deconv_ms2_process::deconvMsTwo and
  // topfd_single_process; TopFD's default), "msdeconv" (MS-Deconv's
  // DeconvSingleSp; TopDIA's default) or "flashdeconv" (FlashDeconvMs2Process
  // instead of DeconvMs2Process)
  std::string ms2_deconv_method_ = "greedy";
  // hidden options --ms1-max-shared-peak-num and --ms2-max-shared-peak-num:
  // the greedy MS1 and MS/MS searches do not report an envelope that shares
  // more than this number of peaks with an envelope of the same charge
  // already reported in the spectrum (toppic_test experiment 157;
  // experiment 158 limited the test to the same charge). The default -1 is
  // no limit (the test is off, the searches before experiment 157): no k
  // raised the identifications on 480_F3_02 (reports 0157 and 0158).
  int ms1_max_shared_peak_num_ = -1;
  int ms2_max_shared_peak_num_ = -1;
  // the precursor masses of an MS/MS spectrum corrected by whole isotopes
  // with the complementary pairs of its fragment masses and the MS1 spectrum
  // (prec_pair_refine::correctIsotopes; toppic_test experiment 116): on by
  // default, off with --skip-ms2-prec-pair-refine and in TopDIA. It runs in
  // DeconvMs2Process (not FLASHDeconv) and only for MS/MS spectra with an
  // MS1 spectrum.
  bool ms2_prec_pair_refine_ = true;
  // topfd --ptcr (isPtcrForced)
  bool ptcr_forced_ = false;
  // the MS3 spectral deconvolution (isMs3Deconv) and the MS3 scans of the
  // file without a proton transfer stage
  bool ms3_deconv_ = true;
  int ms_3_scan_num_ = 0;

  //** information for each run **
  int frac_id_ = -1;
  std::string mzml_file_name_ = "";
  bool is_faims_ = false;
  double faims_volt_ = -1;
  std::string output_base_name_ = "deconv";
  std::string sql_file_name_ = "";
  sqlite3* sql_db_ = nullptr;

  int ms_1_scan_num_ = -1;
  int ms_2_scan_num_ = -1;

  // call function for processing a simple text file
  // containing a mass list
  bool text_peak_list_ = false;
};

}  // namespace toppic

#endif
