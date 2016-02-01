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
      res.push_back(proteoform->getChangePtrVec()[i]);
    }
  }
  return res;
}

} // namespace prot

#endif
