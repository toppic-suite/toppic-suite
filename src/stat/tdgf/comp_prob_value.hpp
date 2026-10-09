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

#ifndef TOPPIC_STAT_TDGF_COMP_PROB_VALUE_HPP_
#define TOPPIC_STAT_TDGF_COMP_PROB_VALUE_HPP_

#include <vector>

#include "common/base/residue_freq.hpp"
#include "ms/spec/prm_peak.hpp"
#include "prsm/prsm.hpp"
#include "stat/tdgf/prob_peak.hpp"

namespace toppic {

class CompProbValue {
 public:
  // strict_complement: the REVERSED (complementary) PRM peaks get their
  // strict tolerance instead of the relaxed one (ProbPeak)
  CompProbValue(double convert_ratio, const ResFreqPtrVec& residue_ptrs,
                double residue_avg_len, int max_layer_num, int max_table_height,
                double max_sp_prec_mass, bool strict_complement = false);

  ~CompProbValue();

  void compute(const ResFreqPtrVec& n_term_residue_ptrs,
               const PrmPeakPtrVec2D& peak_ptr_2d, int thresh, int shift_num,
               bool strict, double prob_prec_mass,
               const PeakTolerancePtr& tole_ptr);
  // the residues of the computations that follow, in place of those of the
  // constructor: the tables are allocated once and serve both residue lists
  // of CompPValueArray (toppic_test experiment 129)
  void setResidues(const ResFreqPtrVec& residue_ptrs, int residue_avg_len);

  // main function to get probabilities
  double getCondProb(int shift, int thresh);
  double getCondProbOneValue(int shift, int value);

 private:
  static int const ORI_BLOCK_LEN = 50;
  static double K1() { return 0.1; }
  static double K2() { return 0.01; }

  // double to integer convert ratio
  double convert_ratio_;

  // the strict tolerance for the REVERSED PRM peaks (ProbPeak)
  bool strict_complement_;

  //**********************************************************
  // Amino acids
  //**********************************************************
  std::vector<int> n_term_acid_masses_;
  std::vector<double> n_term_acid_frequencies_;
  std::vector<int> residue_masses_;
  std::vector<double> residue_frequencies_;
  int residue_avg_len_ = 0;

  //**********************************************************
  // DP Table
  //**********************************************************
  // Number of unexpected mutations
  int max_layer_num_ = 0;
  // Maximum score
  int max_table_height_ = 0;

  // Block length in integer
  int block_len_ = 0;

  // Spectrum
  int max_sp_len_ = 0;
  // the extent of pos_scores_/tmp_pos_scores_ used by this computation
  int pos_len_ = 0;

  std::vector<ProbPeak> prob_peaks_;
  int sp_len_;

  // table height
  int height_;

  std::vector<double> factors_;  // normalization factors;

  int shift_num_;

  std::vector<std::vector<std::vector<double>>>
      results_;  // nLayer peak number, height;
  double shift_prob_;

  // the prob that a randam protein has the same precursor to the spectrum
  std::vector<double> prec_probs_;

  short* pos_scores_;
  short* tmp_pos_scores_;

  void setFactors();

  // computation
  void clearVar();
  void setMassErr(const PrmPeakPtrVec2D& peak_ptr_2d, bool strict,
                  double prec_mass, const PeakTolerancePtr& tole_ptr);

  void updatePosScores(const std::vector<ProbPeak>& prob_peaks,
                       int spectrum_id);

  void setPosScores(const std::vector<ProbPeak>& prob_peaks,
                    int group_spec_num);

  void setVars(int thresh);

  void setPeakBgnEnd(const std::vector<int>& peak_masses,
                     const std::vector<int>& peak_tolerances);

  double getShiftProb();
  void comp();
  void compPrecProbs();

  // The sliced table (toppic_test experiment 150). A layer's table is
  // swept in slices of kChunkRows score rows: the recursion along the mass
  // axis stays within a row (a residue moves the mass, not the score) and
  // only the position score moves a column's values up, by at most the
  // position score, so a slice needs the slice below it only at the
  // columns with a position score (saved_rows_). A slice's working set,
  // the residue look-back window of kChunkRows rows, stays in a core's
  // cache, where the whole-height sweep (one page of 5,000 Da x the table
  // height, 138 MB per thread) streamed the table through memory and the
  // E-value stage was bound by the memory bandwidth of the machine. The
  // per-cell operations and their order are those of the whole-height
  // sweep. Every slice of every layer starts from a cleared page, so the
  // columns before column 0 hold zeros; the whole-height sweep cleared the
  // tail of its page once per computation, and for precursors whose table
  // wrapped the page (above about 4.7 kDa) a shift layer's first blocks
  // read the columns the previous layer left there, which raised the
  // one-shift and two-shift E-values (toppic_test experiment 152).
  static const int kChunkRows = 8;
  int chunk_block_cols_ = 0;   // columns per block (at most the smallest
                               // residue mass, so a block never reads
                               // itself)
  int chunk_page_blocks_ = 0;  // blocks in the circular page
  std::vector<double> chunk_table_;
  std::vector<int> score_cols_;     // the columns with a position score
  int max_pos_score_ = 0;
  std::vector<double> saved_rows_;  // per score column, the slice below's
                                    // last max_pos_score_ rows before the
                                    // column update
  std::vector<double> shift_probs_;  // layer 0's sums over the prior
                                     // columns of the precursor, per row
  void setChunkVars();
  void compOneLayer(std::vector<std::vector<double>>& prev_results,
                    bool is_first_layer,
                    std::vector<std::vector<double>>& cur_results);
  void sweepChunk(std::vector<std::vector<double>>& prev_results,
                  bool is_first_layer,
                  std::vector<std::vector<double>>& cur_results, int s0,
                  int s1);
};

using CompProbValuePtr = std::shared_ptr<CompProbValue>;

}  // namespace toppic
#endif
