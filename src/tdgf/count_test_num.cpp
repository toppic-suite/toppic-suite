#include <cmath>

#include "base/logger.hpp"
#include "tdgf/count_test_num.hpp"
#include "tdgf/comp_prob_value.hpp"

namespace prot {

CountTestNum::CountTestNum(ProteoformPtrVec &raw_forms, 
                           ProteoformPtrVec &prot_mod_forms,
                           ResFreqPtrVec &n_term_residues,
                           ResFreqPtrVec &residues,
                           TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  raw_forms_ = raw_forms;
  prot_mod_forms_ = prot_mod_forms;
  convert_ratio_ = mng_ptr_->double_to_int_constant_;
  max_sp_len_ = (int)std::round(mng_ptr_->max_sp_prec_mass_ * convert_ratio_);
  residue_avg_len_ = computeAvgLength(residues, convert_ratio_);
  norm_factor_ = 0;
  for (unsigned int i = 0; i < n_term_residues.size(); i++) {
    norm_factor_ += n_term_residues[i]->getFreq();
  }
  initCompMassCnt(prot_mod_forms);
  initPrefMassCnt(prot_mod_forms);
  initSuffMassCnt(raw_forms);
  initInternalMassCnt();
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
      suff_mass_cnts_[convertMass(break_points[i]->getSrm())] += 1.0;
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

double CountTestNum::compCandNum(int type, int shift_num, double ori_mass, 
                                 double ori_tolerance) {
  double cand_num = 0;
  if (shift_num == 0) {
    cand_num = compNormNonPtmCandNum(type, shift_num, ori_mass, ori_tolerance);
  }
  // with shifts 
  else if (shift_num == 1){
    cand_num = compOnePtmCandNum(type, shift_num, ori_mass);
  }
  else {
    cand_num = compMultiplePtmCandNum(type, shift_num, ori_mass);
  }

  if (cand_num == 0.0) {
    LOG_WARN("candidate number is ZERO");
  }
  if (type == SEMI_ALIGN_TYPE_PREFIX || type == SEMI_ALIGN_TYPE_SUFFIX) {
    cand_num = cand_num * PREFIX_SUFFIX_ADJUST();
  }
  else if (type == SEMI_ALIGN_TYPE_INTERNAL) {
    cand_num = cand_num * INTERNAL_ADJUST();
  }
  return cand_num;
}

double CountTestNum::compNormNonPtmCandNum(int type, int shift_num, 
                                           double ori_mass, double ori_tolerance) {
  int low = std::floor((ori_mass - ori_tolerance) * convert_ratio_);
  int high = std::ceil((ori_mass + ori_tolerance) * convert_ratio_);
  double cand_num = compSeqNum(type, low, high);
  // normalization: the reason is that we a residue list with frequency sum > 1 in CompProbValue 
  if (type == SEMI_ALIGN_TYPE_COMPLETE || type == SEMI_ALIGN_TYPE_PREFIX) {
    cand_num = cand_num / norm_factor_;
    //System.out.println("nCandidate " + nCandidates + " normFactor " + normFactor);
  } 
  return cand_num;
}

double CountTestNum::compOnePtmCandNum (int type, int shift_num, double ori_mass) {
  double cand_num = 0;
  if (type == SEMI_ALIGN_TYPE_COMPLETE) {
    cand_num = raw_forms_.size();
  } else if (type == SEMI_ALIGN_TYPE_PREFIX) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += raw_forms_[i]->getResSeqPtr()->getLen();
    }
  } else if (type == SEMI_ALIGN_TYPE_SUFFIX) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += raw_forms_[i]->getResSeqPtr()->getLen();
    }
  } else if (type == SEMI_ALIGN_TYPE_INTERNAL) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num = cand_num + raw_forms_[i]->getResSeqPtr()->getLen() 
          * raw_forms_[i]->getResSeqPtr()->getLen();
    }
  }
  return cand_num;
}

double CountTestNum::compSeqNum(int type, int low, int high) {
  double candNum = 0;
  if (type == SEMI_ALIGN_TYPE_COMPLETE) {
    candNum = compMassNum(comp_mass_cnts_, low, high);
  } else if (type == SEMI_ALIGN_TYPE_PREFIX) {
    candNum = compMassNum(pref_mass_cnts_, low, high);
  } else if (type == SEMI_ALIGN_TYPE_SUFFIX) {
    candNum = compMassNum(suff_mass_cnts_, low, high);
  } else if (type == SEMI_ALIGN_TYPE_INTERNAL) {
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

double CountTestNum::compMultiplePtmCandNum (int type, int shift_num, 
                                             double ori_mass) {
  double cand_num = 0;
  if (type == SEMI_ALIGN_TYPE_COMPLETE) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += std::pow(raw_forms_[i]->getResSeqPtr()->getLen(), shift_num - 1);
    }
  } else if (type == SEMI_ALIGN_TYPE_PREFIX) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += std::pow(raw_forms_[i]->getResSeqPtr()->getLen(), shift_num);
    }
  } else if (type == SEMI_ALIGN_TYPE_SUFFIX) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += std::pow(raw_forms_[i]->getResSeqPtr()->getLen(), shift_num);
    }
  } else if (type == SEMI_ALIGN_TYPE_INTERNAL) {
    for (unsigned int i = 0; i < raw_forms_.size(); i++) {
      cand_num += std::pow(raw_forms_[i]->getResSeqPtr()->getLen(), shift_num + 1);
    }
  }
  // use average proportion to estimate the number of matched substring. 
  cand_num = cand_num * getAvgProportion(ori_mass, 
                                         mng_ptr_->sp_para_ptr_->getPeakTolerance()->getPpo(),
                                         convert_ratio_, residue_avg_len_);
  return cand_num;
}

double getAvgProportion(double mass, double ppo, 
                        double convert_ratio, double residue_avg_len) {
  double proportion = mass * ppo * 2  * convert_ratio/ residue_avg_len;
  if (proportion == 0) {
    LOG_ERROR("Error in computing proportion = 0");
  }
  return proportion;
}

}
