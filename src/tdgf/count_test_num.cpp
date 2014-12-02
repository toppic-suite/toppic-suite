#include <cmath>

#include "base/logger.hpp"
#include "base/fasta_reader.hpp"
#include "tdgf/count_test_num.hpp"
#include "tdgf/comp_prob_value.hpp"

namespace prot {
/*
CountTestNum::CountTestNum(const ProteoformPtrVec &raw_proteo_ptrs, 
                           const ProteoformPtrVec &mod_proteo_ptrs,
                           const ResFreqPtrVec &residue_ptrs,
                           double convert_ratio,
                           double max_prec_mass,
                           double max_ptm_mass) {
  raw_proteo_ptrs_ = raw_proteo_ptrs;
  mod_proteo_ptrs_ = mod_proteo_ptrs;
  convert_ratio_ = convert_ratio;
  max_ptm_mass_ = max_ptm_mass;
  max_sp_len_ = (int)std::round(max_prec_mass * convert_ratio_);
  residue_avg_len_ = computeAvgLength(residue_ptrs, convert_ratio_);
  LOG_DEBUG("get residue average length");
  initCompMassCnt(mod_proteo_ptrs);
  LOG_DEBUG("complete mass count initialized");
  initPrefMassCnt(mod_proteo_ptrs);
  LOG_DEBUG("prefix mass count initialized");
  initSuffMassCnt(raw_proteo_ptrs);
  LOG_DEBUG("suffix mass count initialized");
  initInternalMassCnt();
  LOG_DEBUG("internal mass count initialized");
}
*/

CountTestNum::CountTestNum(TdgfMngPtr mng_ptr) {
  convert_ratio_ = mng_ptr->convert_ratio_;
  max_ptm_mass_ = mng_ptr->max_ptm_mass_;
  max_sp_len_ = (int)std::round(mng_ptr->max_prec_mass_ * convert_ratio_);
  init(mng_ptr->prsm_para_ptr_);
  LOG_DEBUG("count numbers initialized");
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

inline int CountTestNum::convertMass(double m) {
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

void updateResidueCounts(const ResiduePtrVec &residue_list, 
                         std::vector<double> &counts,
                         ProteoformPtr prot_ptr) {
  ResSeqPtr seq_ptr = prot_ptr->getResSeqPtr();    
  for (int i = 0; i < seq_ptr->getLen(); i++) {
    ResiduePtr res_ptr = seq_ptr->getResiduePtr(i);
    int pos = findResidue(residue_list, res_ptr);
    if (pos >= 0) {
      // found 
      counts[pos] = counts[pos]+1;
    }
  }
}

ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list, 
                              const std::vector<double> &counts) {
  double sum = 0;
  for (size_t i = 0; i < counts.size(); i++) {
    sum = sum + counts[i];
  }
  ResFreqPtrVec res_freq_list;
  for (size_t i = 0; i < residue_list.size(); i++) {
    ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(), 
                                            residue_list[i]->getPtmPtr(),
                                            counts[i]/sum));
    res_freq_list.push_back(res_freq_ptr);
  }
  return res_freq_list;
}

void updateNTermResidueCounts(ResiduePtrVec &residue_list, std::vector<double> &counts,
                              const ProteoformPtrVec &mod_proteo_ptrs) {
  for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
    ResSeqPtr seq_ptr = mod_proteo_ptrs[i]->getResSeqPtr();    
    if (seq_ptr->getLen() >= 1) {
      ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
      int pos = findResidue(residue_list, res_ptr);
      if (pos >= 0) {
        // found 
        counts[pos] = counts[pos]+1;
      }
      else {
        residue_list.push_back(res_ptr);
        counts.push_back(1);
      }
    }
  }
}

void CountTestNum::init(PrsmParaPtr para_ptr) {
  std::string db_file_name = para_ptr->getSearchDbFileName();
  comp_mass_cnts_ = new double[max_sp_len_](); 
  pref_mass_cnts_ = new double[max_sp_len_]();
  suff_mass_cnts_ = new double[max_sp_len_]();

  ResiduePtrVec residue_list = para_ptr->getFixModResiduePtrVec();
  
  std::vector<double> residue_counts(residue_list.size(), 0.0);

  ResiduePtrVec n_term_residue_list;
  std::vector<double> n_term_residue_counts;

  FastaReader reader(db_file_name);
  ProtModPtrVec prot_mods = para_ptr->getAllowProtModPtrVec();
  ProteoformPtr proteo_ptr = reader.getNextProteoformPtr(residue_list);
  
  ProteoformPtrVec raw_proteo_ptrs 
      = readFastaToProteoform(para_ptr->getSearchDbFileName(), 
                              para_ptr->getFixModResiduePtrVec());

  ResFreqPtrVec residue_freqs 
      = compResidueFreq(para_ptr->getFixModResiduePtrVec(), raw_proteo_ptrs); 
  
  residue_avg_len_ = computeAvgLength(residue_freqs, convert_ratio_);
  
  while (proteo_ptr != nullptr) {
    ProteoformPtrVec mod_proteo_ptrs = generateProtModProteoform(proteo_ptr, prot_mods);
    for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
      // complete
      double m = mod_proteo_ptrs[i]->getResSeqPtr()->getResMassSum();
      comp_mass_cnts_[convertMass(m)] += 1.0;
      // prefix
      std::vector<double> prm_masses = mod_proteo_ptrs[i]->getBpSpecPtr()->getPrmMasses();
      for (size_t j = 1; j < prm_masses.size() - 1; j++) {
        pref_mass_cnts_[convertMass(prm_masses[j])] += 1.0;
      }
    }
    //suffix
    BreakPointPtrVec break_points = proteo_ptr->getBpSpecPtr()->getBreakPointPtrVec();
    for (size_t i = 1; i < break_points.size() - 1; i++) {
      suff_mass_cnts_[convertMass(break_points[i]->getSrm())] += 1.0;
    }

    // length
    raw_proteo_lens_.push_back(proteo_ptr->getLen());
    for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
      mod_proteo_lens_.push_back(mod_proteo_ptrs[i]->getLen());
    }

    // update residue counts 
    updateResidueCounts(residue_list,residue_counts, proteo_ptr);

    // update n terminal residue counts
    updateNTermResidueCounts(n_term_residue_list, n_term_residue_counts, mod_proteo_ptrs);

    // next protein 
    proteo_ptr = reader.getNextProteoformPtr(residue_list);
  }
  // compute residue freq;
  residue_ptrs_ =  compResidueFreq(residue_list, residue_counts);

  // compute residue average length
  residue_avg_len_ = computeAvgLength(residue_ptrs_, convert_ratio_);

  // compute n term residue freq;
  prot_n_term_residue_ptrs_ =  compResidueFreq(n_term_residue_list, n_term_residue_counts);

  // internal 
  initInternalMassCnt();
}

/** initialize the four tables for mass counts */
/*
inline void CountTestNum::initCompMassCnt(const ProteoformPtrVec &mod_proteo_ptrs) {
  // init to 0
  comp_mass_cnts_ = new double[max_sp_len_](); 
  for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
    double m = mod_proteo_ptrs[i]->getResSeqPtr()->getResMassSum();
    //System.out.println("mass " + m);
    comp_mass_cnts_[convertMass(m)] += 1.0;
  }
}
*/

/** initialize the four tables for mass counts */
/*
inline void CountTestNum::initPrefMassCnt(const ProteoformPtrVec &mod_proteo_ptrs) {
  // init to 0
  pref_mass_cnts_ = new double[max_sp_len_]();
  for (size_t i = 0; i < mod_proteo_ptrs.size(); i++) {
    std::vector<double> prm_masses = mod_proteo_ptrs[i]->getBpSpecPtr()->getPrmMasses();
    // prefix
    for (size_t j = 1; j < prm_masses.size() - 1; j++) {
      pref_mass_cnts_[convertMass(prm_masses[j])] += 1.0;
    }
  }
}
*/

/** initialize the four tables for mass counts */
/*
inline void CountTestNum::initSuffMassCnt(const ProteoformPtrVec &raw_proteo_ptrs) {
  // sequence mass 
  // init to 0
  suff_mass_cnts_ = new double[max_sp_len_]();
  for (size_t i = 0; i < raw_proteo_ptrs.size(); i++) {
    BreakPointPtrVec break_points = raw_proteo_ptrs[i]->getBpSpecPtr()->getBreakPointPtrVec();
    // suffix
    for (size_t j = 1; j < break_points.size() - 1; j++) {
      suff_mass_cnts_[convertMass(break_points[j]->getSrm())] += 1.0;
    }
  }
}
*/

inline void CountTestNum::initInternalMassCnt() {
  // init to 0
  internal_mass_cnts_ = new double[max_sp_len_]();
  // middle
  double norm_count = 0;
  // use approxiation to speed up
  LOG_DEBUG("residue_avg_len_ " << residue_avg_len_);
  for (int i = max_sp_len_ - 1; i >= 0; i--) {
    norm_count += suff_mass_cnts_[i];
    internal_mass_cnts_[i] = norm_count/ residue_avg_len_;
  }
}

double CountTestNum::compCandNum(SemiAlignTypePtr type_ptr, int shift_num, double ori_mass, 
                                 double ori_tolerance) {
  double cand_num = 0;
  if (shift_num == 0) {
    cand_num = compNonPtmCandNum(type_ptr, shift_num, ori_mass, ori_tolerance);
  }
  // with shifts 
  else if (shift_num >= 1){
    if (max_ptm_mass_ >=10000) {
      // max shift mass is larger than 10k, we treat it as no limitation 
      cand_num = compPtmCandNum(type_ptr, shift_num, ori_mass);
    }
    else {
      cand_num = compPtmRestrictCandNum(type_ptr, shift_num, ori_mass);
    }
    // multiple adjustment 
    if (type_ptr == SemiAlignTypeFactory::getPrefixPtr() 
        || type_ptr == SemiAlignTypeFactory::getSuffixPtr()) {
      cand_num = cand_num * PREFIX_SUFFIX_ADJUST();
    }
    else if (type_ptr == SemiAlignTypeFactory::getInternalPtr()) {
      cand_num = cand_num * INTERNAL_ADJUST();
    }
  }

  if (cand_num == 0.0) {
    LOG_WARN("candidate number is ZERO");
  }
  return cand_num;
}

double CountTestNum::compNonPtmCandNum(SemiAlignTypePtr type_ptr, int shift_num, 
                                       double ori_mass, double ori_tolerance) {
  int low = std::floor((ori_mass - ori_tolerance) * convert_ratio_);
  int high = std::ceil((ori_mass + ori_tolerance) * convert_ratio_);
  double cand_num = compSeqNum(type_ptr, low, high);
  
  //if (type_ptr == SemiAlignTypeFactory::getCompletePtr()) {
    LOG_DEBUG("low " << low << " high " << high << " cand num " << cand_num);
  //}
  
  return cand_num;
}

double CountTestNum::compPtmCandNum (SemiAlignTypePtr type_ptr, 
                                     int shift_num, double ori_mass) {
  double cand_num = 0;
  if (type_ptr == SemiAlignTypeFactory::getCompletePtr()) {
    cand_num = mod_proteo_lens_.size();
  } else if (type_ptr == SemiAlignTypeFactory::getPrefixPtr()) {
    for (size_t i = 0; i < raw_proteo_lens_.size(); i++) {
      cand_num += mod_proteo_lens_[i];
    }
  } else if (type_ptr == SemiAlignTypeFactory::getSuffixPtr()) {
    for (size_t i = 0; i < raw_proteo_lens_.size(); i++) {
      cand_num += raw_proteo_lens_[i];
    }
  } else if (type_ptr == SemiAlignTypeFactory::getInternalPtr()) {
    for (size_t i = 0; i < raw_proteo_lens_.size(); i++) {
      cand_num = cand_num + raw_proteo_lens_[i] * raw_proteo_lens_[i];
    }
  }
  return cand_num;
}

double CountTestNum::compPtmRestrictCandNum (SemiAlignTypePtr type_ptr, 
                                             int shift_num, double ori_mass) {
  double shift = max_ptm_mass_ * shift_num;
  int low = std::floor((ori_mass - shift) * convert_ratio_);
  int high = std::ceil((ori_mass + shift) * convert_ratio_);
  double cand_num = compSeqNum(type_ptr, low, high);
  return cand_num;
}

double CountTestNum::compSeqNum(SemiAlignTypePtr type_ptr, int low, int high) {
  double candNum = 0;
  if (type_ptr == SemiAlignTypeFactory::getCompletePtr()) {
    candNum = compMassNum(comp_mass_cnts_, low, high);
  } else if (type_ptr == SemiAlignTypeFactory::getPrefixPtr()) {
    candNum = compMassNum(pref_mass_cnts_, low, high);
  } else if (type_ptr == SemiAlignTypeFactory::getSuffixPtr()) {
    candNum = compMassNum(suff_mass_cnts_, low, high);
  } else if (type_ptr == SemiAlignTypeFactory::getInternalPtr()) {
	LOG_DEBUG("internal_mass_cnts_ " << internal_mass_cnts_[low]);
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

int computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio) {
  double mass_sum = 0;
  double freq_sum = 0;
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    freq_sum = freq_sum + residue_ptrs[i]->getFreq();
    mass_sum = mass_sum + residue_ptrs[i]->getFreq() * residue_ptrs[i]->getMass();
  }
  return (int)std::round(mass_sum/freq_sum * convert_ratio);
}

}
