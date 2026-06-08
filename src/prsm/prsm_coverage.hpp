//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_PRSM_PRSM_COVERAGE_HPP_
#define TOPPIC_PRSM_PRSM_COVERAGE_HPP_

#include <string>

#include "prsm/peak_ion_pair.hpp"
#include "prsm/prsm.hpp"
#include "para/prsm_para.hpp"

namespace toppic {

class PrsmCoverage {
 public:
  PrsmCoverage(const PrsmParaPtr &prsm_para_ptr, 
               const std::string &input_file_ext,
               const std::string &output_file_ext);

  void processSingleCoverage();
  void processCombineCoverage();

 private:
  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;

  void printTitle(std::ofstream &file);
  void printTwoTitle(std::ofstream &file);
  void computeCoverage(std::ofstream &file, const PrsmPtr &prsm_ptr,
                       PeakIonPairPtrVec &pair_ptrs, const PrsmParaPtr &prsm_para_ptr);
  void compOneCoverage(std::ofstream &file, const PrsmPtr &prsm,
                       PeakIonPairPtrVec &pair_ptrs, const PrsmParaPtr &prsm_para_ptr);
  void compTwoCoverage(std::ofstream &file, const PrsmPtr &prsm_ptr,
                       PeakIonPairPtrVec &pair_ptrs_1, PeakIonPairPtrVec &pair_ptrs_2,
                       PeakIonPairPtrVec &pair_ptrs_3, const PrsmParaPtr &prsm_para_ptr);
  void processOnePrsm(std::ofstream &file, const PrsmPtr &prsm_ptr,
                      const PrsmParaPtr &prsm_para_ptr);
  void processTwoPrsms(std::ofstream &file, const PrsmPtr &prsm_ptr_1,
                       const PrsmPtr &prsm_ptr_2, const PrsmParaPtr &prsm_para_ptr);

  void outputMatchPeaks(const PrsmPtr &prsm_ptr, const PrsmParaPtr &prsm_para_ptr);
};

using PrsmCoveragePtr = std::shared_ptr<PrsmCoverage>;

}  // namespace toppic

#endif

