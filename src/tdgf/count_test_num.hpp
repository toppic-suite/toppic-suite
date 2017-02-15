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


#ifndef PROT_COUNT_TEST_NUM_HPP_
#define PROT_COUNT_TEST_NUM_HPP_

#include "base/proteoform.hpp"
#include "base/residue_freq.hpp"
#include "tdgf/tdgf_mng.hpp"

namespace prot {

class CountTestNum {
 public:
  /*
  CountTestNum(const ProteoformPtrVec &raw_proteo_ptrs, 
               const ProteoformPtrVec &mod_proteo_ptrs,
               const ResFreqPtrVec &residue_ptrs, 
               double convert_ratio, double max_prec_mass,
               double max_ptm_mass);
               */
  CountTestNum(TdgfMngPtr mng_ptr);

  ~CountTestNum();

  double compCandNum(AlignTypePtr type_ptr, int index, 
                     double ori_mass, double ori_tolerance);

  ResFreqPtrVec getResFreqPtrVec() {return residue_ptrs_;}
  ResFreqPtrVec getNTermResFreqPtrVec() {return prot_n_term_residue_ptrs_;}
  double getResidueAvgLen() {return residue_avg_len_;}

 private:
  static double PREFIX_SUFFIX_ADJUST() {return 0.693;}
  static double INTERNAL_ADJUST() {return 0.508;}

  //ProteoformPtrVec raw_proteo_ptrs_;
  //ProteoformPtrVec mod_proteo_ptrs_;

  double *comp_mass_cnts_;
  double *pref_mass_cnts_;
  double *suff_mass_cnts_;
  double *internal_mass_cnts_;

  ResFreqPtrVec residue_ptrs_; 
  ResFreqPtrVec prot_n_term_residue_ptrs_;

  std::vector<int> raw_proteo_lens_;
  std::vector<int> mod_proteo_lens_;

  double convert_ratio_;
  int max_sp_len_;
  int residue_avg_len_;
  double max_ptm_mass_;
  double norm_factor_;
  

  int convertMass(double m);
  void init(PrsmParaPtr para_ptr);
  void initCompMassCnt(const ProteoformPtrVec &prot_mod_forms);
  void initPrefMassCnt(const ProteoformPtrVec &prot_mod_forms);
  void initSuffMassCnt(const ProteoformPtrVec &raw_forms);
  void initInternalMassCnt();

  double compNonPtmCandNum(AlignTypePtr type_ptr, 
                           double ori_mass, double ori_tolerance);
  double compPtmCandNum (AlignTypePtr type_ptr, double ori_mass);
  double compPtmRestrictCandNum (AlignTypePtr type_ptr, int shift_num, double ori_mass);
  double compSeqNum(AlignTypePtr type_ptr, int low, int high);
  double compMassNum(double *cnts, int low, int high);
};

typedef std::shared_ptr<CountTestNum> CountTestNumPtr;


}

#endif
