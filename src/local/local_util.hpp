//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_LOCAL_UTIL_HPP_
#define PROT_LOCAL_UTIL_HPP_

#include "base/ptm.hpp"
#include "base/acid_base.hpp"
#include "base/change.hpp"
#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair_util.hpp"
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

  static int compNumPeakIonPairs(const ProteoformPtr &proteoform_ptr, const ExtendMsPtrVec &ms_ptr_vec) {
    return peak_ion_pair_util::genePeakIonPairs(proteoform_ptr, ms_ptr_vec, mng_ptr_->min_mass_).size();
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
