#ifndef PROT_COUNT_TEST_NUM_HPP_
#define PROT_COUNT_TEST_NUM_HPP_

#include "base/proteoform.hpp"
#include "base/residue_freq.hpp"
#include "tdgf/tdgf_mng.hpp"

namespace prot {

class CountTestNum {
 public:
  CountTestNum(ProteoformPtrVec &raw_forms, ProteoformPtrVec &prot_mod_forms,
               ResFreqPtrVec &residues, double convert_ratio, double max_prec_mass,
               double max_ptm_mass);

  ~CountTestNum();

  double compCandNum(SemiAlignTypePtr type, int shift_num, 
                     double ori_mass, double ori_tolerance);

 private:
  static double PREFIX_SUFFIX_ADJUST() {return 0.693;}
  static double INTERNAL_ADJUST() {return 0.508;}


  ProteoformPtrVec raw_forms_;
  ProteoformPtrVec prot_mod_forms_;

  double *comp_mass_cnts_;
  double *pref_mass_cnts_;
  double *suff_mass_cnts_;
  double *internal_mass_cnts_;

  double convert_ratio_;
  int max_sp_len_;
  int residue_avg_len_;
  double max_ptm_mass_;
  double norm_factor_;

  int convertMass(double m);
  void initCompMassCnt(ProteoformPtrVec &prot_mod_forms);
  void initPrefMassCnt(ProteoformPtrVec &prot_mod_forms);
  void initSuffMassCnt(ProteoformPtrVec &raw_forms);
  void initInternalMassCnt();
  double compNonPtmCandNum(SemiAlignTypePtr type, int shift_num, 
                           double ori_mass, double ori_tolerance);
  double compPtmCandNum (SemiAlignTypePtr type, int shift_num, double ori_mass);
  double compPtmRestrictCandNum (SemiAlignTypePtr type, int shift_num, double ori_mass);
  double compSeqNum(SemiAlignTypePtr type, int low, int high);
  double compMassNum(double *cnts, int low, int high);
};

typedef std::shared_ptr<CountTestNum> CountTestNumPtr;

}

#endif
