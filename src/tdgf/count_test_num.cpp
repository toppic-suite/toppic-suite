#include <cmath>

#include "base/logger.hpp"
#include "tdgf/count_test_num.hpp"
#include "tdgf/comp_prob_value.hpp"

namespace prot {

CountTestNum::CountTestNum(ProteoformPtrVec &raw_forms, 
                           ProteoformPtrVec &prot_mod_forms,
                           ResFreqPtrVec &residues,
                           double convert_ratio,
                           double max_prec_mass,
                           double max_ptm_mass) {
  raw_forms_ = raw_forms;
  prot_mod_forms_ = prot_mod_forms;
  convert_ratio_ = convert_ratio;
  max_ptm_mass_ = max_ptm_mass;
  max_sp_len_ = (int)std::round(max_prec_mass * convert_ratio_);
  residue_avg_len_ = computeAvgLength(residues, convert_ratio_);
  LOG_DEBUG("get residue average length");
  initCompMassCnt(prot_mod_forms);
  LOG_DEBUG("complete mass count initialized");
  initPrefMassCnt(prot_mod_forms);
  LOG_DEBUG("prefix mass count initialized");
  initSuffMassCnt(raw_forms);
  LOG_DEBUG("suffix mass count initialized");
  initInternalMassCnt();
  LOG_DEBUG("internal mass count initialized");
}

CountTestNum::~CountTestNum() {
  if (comp_mass_cnts_ != nullptr) {
    delete comp_mass_cnts_;
  }
  if (pref_mass_cnts_ != nullptr) {
    delete pref_mass_cnts_;
  }
  if (suff_mass_cnts_ != nullptr) {
    delete suff_mass_cnts_;
  }
  if (internal_mass_cnts_ != nullptr) {
    delete internal_mass_cnts_;
  }
}

int CountTestNum::convertMass(double m) {
  int n = (int) std::round(m * convert_ratio_);
  if (n < 0) {
    LOG_WARN("Negative mass value: " << m);
    return 0;
  }
  if (n >= max_sp_len_) {
    n = max_sp_len_ - 1;
  }
  return n;
}

/** initialize the four tables for mass counts */
void CountTestNum::initCompMassCnt(ProteoformPtrVec &prot_mod_forms) {
  comp_mass_cnts_ = new double[max_sp_len_]; 
  for (unsigned int i = 0; i < prot_mod_forms.size(); i++) {
    double m = prot_mod_forms[i]->getResSeqPtr()->getResMassSum();
    //System.out.println("mass " + m);
    comp_mass_cnts_[convertMass(m)] += 1.0;
  }
}

/** initialize the four tables for mass counts */
void CountTestNum::initPrefMassCnt(ProteoformPtrVec &prot_mod_forms) {
  pref_mass_cnts_ = new double[max_sp_len_];
  for (unsigned int i = 0; i < prot_mod_forms.size(); i++) {
    std::vector<double> prm_masses = prot_mod_forms[i]->getBpSpecPtr()->getPrmMasses();
    // prefix
    for (unsigned int j = 1; j < prm_masses.size() - 1; j++) {
      pref_mass_cnts_[convertMass(prm_masses[j])] += 1.0;
    }
  }
}

/** initialize the four tables for mass counts */
void CountTestNum::initSuffMassCnt(ProteoformPtrVec &raw_forms) {
  // sequence mass 
  suff_mass_cnts_ = new double[max_sp_len_];
  for (unsigned int i = 0; i < raw_forms.size(); i++) {
    BreakPointPtrVec break_points = raw_forms[i]->getBpSpecPtr()->getBreakPointPtrVec();
    // suffix
    for (unsigned int j = 1; j < break_points.size() - 1; j++) {
      suff_mass_cnts_[convertMass(break_points[j]->getSrm())] += 1.0;
    }
  }
}

void CountTestNum::initInternalMassCnt() {
  internal_mass_cnts_ = new double[max_sp_len_];
  // middle
  double norm_count = 0;
  // use approxiation to speed up
  for (int i = max_sp_len_ - 1; i >= 0; i--) {
    norm_count += suff_mass_cnts_[i];
    internal_mass_cnts_[i] = norm_count/ residue_avg_len_;
  }
}

double CountTestNum::compCandNum(SemiAlignTypePtr type, int shift_num, double ori_mass, 
                                 double ori_tolerance) {
  double cand_num = 0;
  if (shift_num == 0) {
    cand_num = compNonPtmCandNum(type, shift_num, ori_mass, ori_tolerance);
  }
  // with shifts 
  else if (shift_num >= 1){
    if (max_ptm_mass_ >=10000) {
      // max shift mass is larger than 10k, we treat it as no limitation 
      cand_num = compPtmCandNum(type, shift_num, ori_mass);
    }
    else {
      cand_num = compPtmRestrictCandNum(type, shift_num, ori_mass);
    }
    // multiple adjustment 
    if (type == SemiAlignTypeFactory::getPrefixPtr() 
        || type == SemiAlignTypeFactory::getSuffixPtr()) {
      cand_num = cand_num * PREFIX_SUFFIX_ADJUST();
    }
    else if (type == SemiAlignTypeFactory::getInternalPtr()) {
      cand_num = cand_num * INTERNAL_ADJUST();
    }
  }

  if (cand_num == 0.0) {
    LOG_WARN("candidate number is ZERO");
  }
  return cand_num;
}

double CountTestNum::compNonPtmCandNum(SemiAlignTypePtr type, int shift_num, 
                                       double ori_mass, double ori_tolerance) {
  int low = std::floor((ori_mass - ori_tolerance) * convert_ratio_);
  int high = std::ceil((ori_mass + ori_tolerance) * convert_ratio_);
  double cand_num = compSeqNum(type, low, high);
  return cand_num;
}

double CountTestNum::compPtmCandNum (SemiAlignTypePtr type, 
                                     int shift_num, double ori_mass) {
  double cand_num = 0;
  if (type == SemiAlignTypeFactory::getCompletePtr()) {
    cand_num = prot_mod_forms_.size();
  } else if (type == SemiAlignTypeFactory::getPrefixPtr()) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += prot_mod_forms_[i]->getResSeqPtr()->getLen();
    }
  } else if (type == SemiAlignTypeFactory::getSuffixPtr()) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += raw_forms_[i]->getResSeqPtr()->getLen();
    }
  } else if (type == SemiAlignTypeFactory::getInternalPtr()) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num = cand_num + raw_forms_[i]->getResSeqPtr()->getLen() 
          * raw_forms_[i]->getResSeqPtr()->getLen();
    }
  }
  return cand_num;
}

double CountTestNum::compPtmRestrictCandNum (SemiAlignTypePtr type, 
                                             int shift_num, double ori_mass) {
  double shift = max_ptm_mass_ * shift_num;
  int low = std::floor((ori_mass - shift) * convert_ratio_);
  int high = std::ceil((ori_mass + shift) * convert_ratio_);
  double cand_num = compSeqNum(type, low, high);
  return cand_num;
}

double CountTestNum::compSeqNum(SemiAlignTypePtr type, int low, int high) {
  double candNum = 0;
  if (type == SemiAlignTypeFactory::getCompletePtr()) {
    candNum = compMassNum(comp_mass_cnts_, low, high);
  } else if (type == SemiAlignTypeFactory::getPrefixPtr()) {
    candNum = compMassNum(pref_mass_cnts_, low, high);
  } else if (type == SemiAlignTypeFactory::getSuffixPtr()) {
    candNum = compMassNum(suff_mass_cnts_, low, high);
  } else if (type == SemiAlignTypeFactory::getInternalPtr()) {
    candNum = compMassNum(internal_mass_cnts_, low, high);
  }
  return candNum;
}

double CountTestNum::compMassNum(double *cnts, int low, int high) {
  double cnt = 0;
  if (high >= max_sp_len_) {
    high = max_sp_len_ - 1;
  }
  if (low < 0) {
    low = 0;
  }
  if (low > high) {
    low = high;
  }
  for (int i = low; i <= high; i++) {
    cnt += cnts[i];
  }
  return cnt;
}

}
