#ifndef PROT_COMP_PROB_VALUE_HPP_
#define PROT_COMP_PROB_VALUE_HPP_

#include <vector>

#include "base/residue_freq.hpp"
#include "spec/prm_peak.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class CompProbValue {
 public:
  CompProbValue(double convert_ratio, const ResFreqPtrVec &residue_ptrs, 
                int max_layer_num, int max_table_height, 
                double max_sp_prec_mass);

  ~CompProbValue();

  void compute(const ResFreqPtrVec &n_term_residue_ptrs, 
               const PrmPeakPtrVec &peak_ptrs, 
               int thresh, int shift_num, bool strict);
  // main function to get probabilities
  double getCondProb(int shift, int thresh);
  double getCondProbOneValue(int shift, int value);

 private:
  static int const ORI_PAGE_LEN = 5000;
  static int const ORI_BLOCK_LEN = 50;
  static double K() {return 0.55;}

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

  std::vector<int> peak_masses_;
  std::vector<int> peak_tolerances_;
  std::vector<int> base_types_;
  int sp_len_;

  /** table height */
  int height_;

  int sp_table_size_;      // spLen * height
  int page_table_size_;    // pageLen * height
  int block_table_size_;   // blockLen * height

  std::vector<int> acid_dists_;   // acidMass * height;
  std::vector<double> factors_; //normalization factors;

  std::vector<int> peak_mass_bgns_;
  std::vector<int> peak_mass_ends_;
  std::vector<int> peak_table_bgns_;
  std::vector<int> peak_table_ends_;

  int shift_num_;

  std::vector<std::vector<std::vector<double>>> results_; // nLayer peak number, height;
  std::vector<std::vector<std::vector<double>>> priors_;  // peak number, height

  // the prob that a randam protein has the same precursor to the spectrum 
  std::vector<double> prec_probs_; 

  double* page_table_;
  short* pos_scores_;

  void setFactors();

  // computation
  void clearVar();
  void setMassErr(const PrmPeakPtrVec &peak_ptrs, bool strict);

  void setPosScores(const std::vector<int> &peak_masses, 
                    const std::vector<int> &peak_tolerances,
                    const std::vector<int> &base_types);

  void setHeight(int thresh, int max_peak_mass);
  void setPeakBgnEnd(const std::vector<int> &peak_masses, 
                     const std::vector<int> &peak_tolerances);
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
                    std::vector<std::vector<double>> &cur_results, 
                    std::vector<std::vector<double>> &cur_priors);
};

typedef std::shared_ptr<CompProbValue> CompProbValuePtr;

int computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio);

void compProbArray(CompProbValuePtr comp_prob_ptr, const ResFreqPtrVec &n_term_residue_ptrs, 
                   const PrmPeakPtrVec &peak_ptrs, const PrsmPtrVec &prsm_ptrs, bool strict, 
                   std::vector<double> &results);

}
#endif
