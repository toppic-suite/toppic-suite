// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_COMP_PROB_VALUE_HPP_
#define PROT_COMP_PROB_VALUE_HPP_

#include <vector>

#include "base/residue_freq.hpp"
#include "spec/prm_peak.hpp"
#include "spec/base_peak_type.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class ProbPeak {
 public:
  ProbPeak(PrmPeakPtr peak_ptr, int spectrum_id, int height, 
           bool strict, double convert_ratio);
  int mass_;
  int tolerance_;
  BasePeakTypePtr base_type_ptr_;
  int spectrum_id_;
  int mass_bgn_;
  int mass_end_;
  int table_bgn_;
  int table_end_;
};

class CompProbValue;
typedef std::shared_ptr<CompProbValue> CompProbValuePtr;

class CompProbValue {
 public:
  CompProbValue(double convert_ratio, const ResFreqPtrVec &residue_ptrs, 
                double residue_avg_len, int max_layer_num, int max_table_height, 
                double max_sp_prec_mass);

  ~CompProbValue();

  void compute(const ResFreqPtrVec &n_term_residue_ptrs, 
               const PrmPeakPtrVec2D &peak_ptr_2d, 
               int thresh, int shift_num, bool strict);
  // main function to get probabilities
  double getCondProb(int shift, int thresh);
  double getCondProbOneValue(int shift, int value);

  static void compProbArray(CompProbValuePtr comp_prob_ptr, const ResFreqPtrVec &n_term_residue_ptrs, 
                            const PrmPeakPtrVec2D &peak_ptr_2d, const PrsmPtrVec &prsm_ptrs, bool strict, 
                            std::vector<double> &results);

 private:
  static int const ORI_PAGE_LEN = 5000;
  static int const ORI_BLOCK_LEN = 50;
  static double K1() {return 0.1;}
  static double K2() {return 0.1;}

  /* double to integer convert ratio */
  double convert_ratio_;

  /**********************************************************
   * Amino acids
   **********************************************************/
  std::vector<int> n_term_acid_masses_;
  std::vector<double> n_term_acid_frequencies_;
  std::vector<int> residue_masses_;
  std::vector<double> residue_frequencies_;
  int residue_avg_len_ = 0;

  /**********************************************************
   * DP Table
   **********************************************************/
  /** number of unexpected mutations */
  int max_layer_num_ = 0;
  /** maximum score */
  int max_table_height_ = 0;

  /** page length in integer */
  int page_len_ = 0;
  int block_len_ = 0;

  /** spectrum */
  int max_sp_len_ = 0;

  std::vector<ProbPeak> prob_peaks_;
  int sp_len_;

  /** table height */
  int height_;

  int sp_table_size_;      // spLen * height
  int page_table_size_;    // pageLen * height
  int block_table_size_;   // blockLen * height

  std::vector<int> acid_dists_;   // acidMass * height;
  std::vector<double> factors_; //normalization factors;

  int shift_num_;

  std::vector<std::vector<std::vector<double>>> results_; // nLayer peak number, height;
  double shift_prob_;  

  // the prob that a randam protein has the same precursor to the spectrum 
  std::vector<double> prec_probs_; 

  double* page_table_;
  short* pos_scores_;
  short* tmp_pos_scores_;

  void setFactors();

  // computation
  void clearVar();
  void setMassErr(const PrmPeakPtrVec2D &peak_ptr_2d, bool strict);

  void updatePosScores(const std::vector<ProbPeak> &prob_peaks, 
                       int spectrum_id);

  void setPosScores(const std::vector<ProbPeak> &prob_peaks, 
                    int group_spec_num);

  void setVars(int thresh);

  void setPeakBgnEnd(const std::vector<int> &peak_masses, 
                     const std::vector<int> &peak_tolerances);

  double getShiftProb();
  void comp();
  void compPrecProbs();
  void runClear(int page_pos);
  void runFirstLayerInit(int win_table_bgn, int win_table_end);
  void runInit(std::vector<std::vector<double>> &results, 
               int win_table_bgn, int win_table_end);
  void runAddProb(int page_pos, int prev_pos, int size, double f);
  void updateCol(int col_end, int scr);
  void runUpdate(int page_end, int scr_bgn, int scr_end);
  void compOneLayer(std::vector<std::vector<double>> &prev_results, 
                    bool is_first_layer,  
                    std::vector<std::vector<double>> &cur_results);

};

}
#endif
