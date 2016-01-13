#ifndef PROT_LOCAL_UTIL_HPP_
#define PROT_LOCAL_UTIL_HPP_

#include "base/ptm.hpp"
#include "base/algorithm.hpp"
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
  static void compSupPeakNum(ProteoformPtr proteoform,
                             ExtendMsPtrVec extend_ms_ptr_vec,
                             ChangePtr change,
                             double min_mass, int & left, int & right);

  static PtmPtrVec getPtmPtrVecByMass(double mass, double err);

  static PtmPairVec getPtmPairVecByMass(double mass1, double mass2, double err);

  static bool modifiable(ProteoformPtr proteoform_ptr, int i, PtmPtr ptm_ptr);

  static void compOnePtmScr(ProteoformPtr proteoform, ExtendMsPtrVec extend_ms_ptr_vec, 
                            std::vector<double> &scr_vec, double & raw_scr, PtmPtrVec & ptm_vec);

  static void compTwoPtmScr(ProteoformPtr proteoform, int num_match, 
                            ExtendMsPtrVec extend_ms_ptr_vec, double prec_mass,
                            double & raw_scr, PtmPairVec ptm_pair_vec);

  static double dpTwoPtmScr(ProteoformPtr proteoform, int num_match, 
                            ExtendMsPtrVec extend_ms_ptr_vec, double prec_mass,
                            double mass1, double mass2, PtmPtr p1, PtmPtr p2);

  static void onePtmTermAdjust(ProteoformPtr proteoform, ExtendMsPtrVec extend_ms_ptr_vec,
                               double & mass, double err);

  static void twoPtmTermAdjust(ProteoformPtr proteoform, int num_match, 
                               ExtendMsPtrVec extend_ms_ptr_vec, double prec_mass,
                               double & mass1, double & mass2);

  static void getNtermTruncRange(ProteoformPtr proteoform, int & min, int & max, double max_mass);

  static void getCtermTruncRange(ProteoformPtr proteoform, int & min, int & max, double max_mass);

  static void compSplitPoint(ProteoformPtr proteoform, int num_match, ExtendMsPtrVec extend_ms_ptr_vec,
                             double prec_mass);

  static int compNumPeakIonPairs(const ProteoformPtr &proteoform_ptr, 
                                 const ExtendMsPtrVec &ms_ptr_vec);

 private:
  static void readPtmTxt(const std::string &file_name);

  static LocalMngPtr mng_ptr_;
  static PtmPtrVec var_ptm_list_;
  static PtmPairVec ptm_pair_vec_;
  static ModPtrVec mod_list_N_;
  static ModPtrVec mod_list_C_;
  static ModPtrVec mod_list_any_;

  static double ppm_, p1_, p2_;

};
} // namespace prot

#endif
