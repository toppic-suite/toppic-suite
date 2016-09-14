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


#ifndef PROT_LOCAL_UTIL_HPP_
#define PROT_LOCAL_UTIL_HPP_

#include "base/ptm.hpp"
#include "base/algorithm.hpp"
#include "base/acid_base.hpp"
#include "base/change.hpp"
#include "prsm/prsm.hpp"
#include "local_mng.hpp"

namespace prot {

class LocalUtil {
 public:
  static void init(LocalMngPtr mng_ptr);

  static std::vector<double> normalize(const std::vector<double> & scr);

  static void scr_filter(std::vector<double> & scr, int & bgn, int & end,
                         double & conf, double thread);

  // compute the number of supporting peaks at the left and right of the change
  static void compSupPeakNum(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                             ChangePtr change, double min_mass, int & left, int & right);

  static PtmPtrVec getPtmPtrVecByMass(double mass, double err);

  static PtmPairVec getPtmPairVecByMass(double mass1, double mass2, double err);

  static bool modifiable(ProteoformPtr proteoform_ptr, int i, PtmPtr ptm_ptr);

  static void compOnePtmScr(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec, 
                            std::vector<double> &scr_vec, double & raw_scr, PtmPtrVec & ptm_vec);

  static void compTwoPtmScr(ProteoformPtr proteoform, int num_match, 
                            const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                            double & raw_scr, PtmPairVec & ptm_pair_vec);

  static double dpTwoPtmScr(ProteoformPtr proteoform, int num_match, 
                            const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                            double mass1, double mass2, PtmPtr p1, PtmPtr p2);

  static void onePtmTermAdjust(ProteoformPtr & proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                               double & mass, double err);

  static void twoPtmTermAdjust(ProteoformPtr & proteoform, int num_match, 
                               const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                               double & mass1, double & mass2);

  static void compSplitPoint(ProteoformPtr & proteoform, int num_match, const ExtendMsPtrVec & extend_ms_ptr_vec,
                             double prec_mass);

  static int compNumPeakIonPairs(const ProteoformPtr &proteoform_ptr, const ExtendMsPtrVec &ms_ptr_vec){
    return PeakIonPairFactory::genePeakIonPairs(proteoform_ptr, ms_ptr_vec, mng_ptr_->min_mass_).size();
  }

 private:
  static void readPtmTxt(const std::string &file_name);

  static void getNtermTruncRange(ProteoformPtr proteoform, const ExtendMsPtrVec &ms_ptr_vec, 
                                 int & min, int & max);

  static void getCtermTruncRange(ProteoformPtr proteoform, const ExtendMsPtrVec &ms_ptr_vec,
                                 int & min, int & max);

  static LocalMngPtr mng_ptr_;
  static PtmPtrVec var_ptm_list_;
  static PtmPairVec ptm_pair_vec_;
  static ModPtrVec mod_list_N_;
  static ModPtrVec mod_list_C_;
  static ModPtrVec mod_list_any_;

  static double ppm_, p1_, p2_;

};

inline ChangePtr geneUnexpectedChange(ChangePtr change, double mass) {
  return std::make_shared<Change>(change->getLeftBpPos(), 
                                  change->getRightBpPos(), 
                                  ChangeType::UNEXPECTED, mass, 
                                  std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(), 
                                                        ResidueBase::getEmptyResiduePtr()));
}

inline double getPeptideMass(const std::string & seq) {
  double m = 0;
  for (size_t i = 0; i < seq.length(); i++) {
    m += AcidBase::getAcidPtrByOneLetter(seq.substr(i, 1))->getMonoMass();
  }
  return m;
}

inline ChangePtrVec getExpectedChangeVec(ProteoformPtr proteoform) {
  ChangePtrVec res;
  for (size_t i = 0; i < proteoform->getChangePtrVec().size(); i++) {
    if (proteoform->getChangePtrVec()[i]->getChangeTypePtr() != ChangeType::UNEXPECTED) {
      ChangePtr tmp = proteoform->getChangePtrVec()[i];
      res.push_back(std::make_shared<Change>(tmp->getLeftBpPos(), 
                                             tmp->getRightBpPos(), 
                                             tmp->getChangeTypePtr(), tmp->getMassShift(), 
                                             tmp->getModPtr()));
    }
  }
  return res;
}

} // namespace prot

#endif
