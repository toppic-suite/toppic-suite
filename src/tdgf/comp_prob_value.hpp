//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef TOPPIC_TDGF_COMP_PROB_VALUE_HPP_
#define TOPPIC_TDGF_COMP_PROB_VALUE_HPP_

#include <vector>

#include "common/base/residue_freq.hpp"
#include "spec/prm_peak.hpp"
#include "spec/base_peak_type.hpp"
#include "prsm/prsm.hpp"

namespace toppic {

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
               int thresh, int shift_num, bool strict, 
               double prob_prec_mass, PeakTolerancePtr tole_ptr);
  // main function to get probabilities
  double getCondProb(int shift, int thresh);
  double getCondProbOneValue(int shift, int value);

  static void compProbArray(CompProbValuePtr comp_prob_ptr, const ResFreqPtrVec &n_term_residue_ptrs, 
                            const PrmPeakPtrVec2D &peak_ptr_2d, const PrsmPtrVec &prsm_ptrs, bool strict,
                            double prob_prec_mass, PeakTolerancePtr tole_ptr, std::vector<double> &results);

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
  void setMassErr(const PrmPeakPtrVec2D &peak_ptr_2d, bool strict, 
                  double prec_mass, PeakTolerancePtr tole_ptr);

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
