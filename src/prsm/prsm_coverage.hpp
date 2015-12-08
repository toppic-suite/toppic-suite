#ifndef PROT_PRSM_PRSM_COVERAGE_HPP_
#define PROT_PRSM_PRSM_COVERAGE_HPP_

#include "base/string_util.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class PrsmCoverage {
 public:
  PrsmCoverage(PrsmParaPtr prsm_para_ptr, const std::string &input_file_ext,
               const std::string &output_file_ext);
  void processSingleCoverage();
  void processCombineCoverage();

 private:
  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;

  void printTitle(std::ofstream &file);
  void printTwoTitle(std::ofstream &file);
  void computeCoverage(std::ofstream &file, PrsmPtr prsm_ptr, 
                       PeakIonPairPtrVec &pair_ptrs, PrsmParaPtr prsm_para_ptr);
  void compOneCoverage(std::ofstream &file, PrsmPtr prsm, 
                       PeakIonPairPtrVec &pair_ptrs, PrsmParaPtr prsm_para_ptr);
  void compTwoCoverage(std::ofstream &file, PrsmPtr prsm_ptr, 
                    PeakIonPairPtrVec &pair_ptrs_1, PeakIonPairPtrVec &pair_ptrs_2,
                    PeakIonPairPtrVec &pair_ptrs_3, PrsmParaPtr prsm_para_ptr);
  void processOnePrsm(std::ofstream &file, PrsmPtr prsm_ptr, 
                      PrsmParaPtr prsm_para_ptr);
  void processTwoPrsms(std::ofstream &file, PrsmPtr prsm_ptr_1, 
                       PrsmPtr prsm_ptr_2, PrsmParaPtr prsm_para_ptr);

};

typedef std::shared_ptr<PrsmCoverage> PrsmCoveragePtr;

} /* namespace prot */

#endif /* PROT_PRSM_COVERAGE_HPP_ */
