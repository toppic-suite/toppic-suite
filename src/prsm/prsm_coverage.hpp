// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
