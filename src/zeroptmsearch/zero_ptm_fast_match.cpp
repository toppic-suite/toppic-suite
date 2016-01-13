#include "base/logger.hpp"
#include "base/base_data.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"

namespace prot {

ZeroPtmFastMatch::ZeroPtmFastMatch (ProteoformPtr proteo_ptr, 
                                    double score, int begin, int end) {
  proteo_ptr_ = proteo_ptr;
  score_ = score;
  begin_ = begin;
  end_ = end;
}

/*
 * in the computing of diagonal score in fast filtering, we allow to use n
 * terminal large error tolerance
 */
double compDiagScr(ExtendMsPtr ms_ptr,
                   const std::vector<double> &masses, double center) {
  size_t i = 0;
  size_t j = 0;
  double s = 0;
  while (i < ms_ptr->size() && j < masses.size()) {
    ExtendPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    double distance = peak_ptr->getMonoMass() - masses[j];
    if (std::abs(center - distance) <= peak_ptr->getOrigTolerance()) {
      s += peak_ptr->getScore();
      i++;
      j++;
    }
    /*
     * we use 1 here since the difference between consecutive mass_b is
     * at least 50. and some time masses[i+1] may have small error tolerance
     * than masses[i]
     */
    if (distance > center + 1) {
      j++;
    } else {
      i++;
    }
  }
  return s;
}

double compScore(ExtendMsPtr ms_ptr, ProteoformPtr proteo_ptr, double n_shift,
                 double c_shift) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  ActivationPtr activation = header_ptr->getActivationPtr();
  IonTypePtr n_ion_type_ptr = activation->getNIonTypePtr();
  std::vector<double> masses 
      = proteo_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type_ptr);
  double score = compDiagScr(ms_ptr, masses, n_shift);

  IonTypePtr c_ion_type_ptr = activation->getCIonTypePtr();
  masses = proteo_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type_ptr);
  score += compDiagScr(ms_ptr, masses, c_shift);
  LOG_TRACE("score " << score);
  return score;
}

double compScore(const ExtendMsPtrVec &ms_ptr_vec, ProteoformPtr proteo_ptr, 
                 double n_shift, double c_shift) {
  double score = 0;
  for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
    score = score + compScore(ms_ptr_vec[i], proteo_ptr, n_shift, c_shift);
  }
  return score;
}


ZpFastMatchPtr computeCompMatch(const ExtendMsPtrVec &ms_ptr_vec, 
                                ProteoformPtr proteo_ptr, double ppo) {
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error = header_ptr->getErrorTolerance(ppo);
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();
  double prot_mass = proteo_ptr->getResSeqPtr()->getResMassSum();
  double error = std::abs(res_sum_mass - prot_mass);
  LOG_TRACE("complete protein mass " << prot_mass 
            << " precursor mass " << res_sum_mass 
            << " proteoform name " << proteo_ptr->getSeqName()
            << " error " << error << " error tolerance " << max_error);
  double score = 0;
  if (error <= max_error) {
    score = compScore(ms_ptr_vec, proteo_ptr, 0, 0);
  }
  int end = proteo_ptr->getResSeqPtr()->getLen() - 1;
  return ZpFastMatchPtr(
      new ZeroPtmFastMatch(proteo_ptr, score, 0, end));
}

ZpFastMatchPtr computePrefixMatch(const ExtendMsPtrVec &ms_ptr_vec, 
                                  ProteoformPtr proteo_ptr, double ppo) {
  /* check if there is a matched prefix */
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error = header_ptr->getErrorTolerance(ppo);
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();

  bool is_prefix = false;
  int seq_end = 0;
  for (size_t i = 0; i < prms.size() - 1; i++) {
    /*
    LOG_TRACE("residue sum mass " << res_sum_mass << " prsm  " << prms[i] 
              << " residue " << proteo_ptr->getResSeqPtr()->getResiduePtr(i)->toString()
              << " error " << std::abs(res_sum_mass - prms[i])
              << " max error " << max_error << " "
              << proteo_ptr->getDbResSeqPtr()->getName());
              */
    if (std::abs(res_sum_mass - prms[i]) <= max_error) {
      is_prefix = true;
      seq_end = i - 1;
      break;
    } else {
      if (prms[i] > res_sum_mass) {
        break;
      }
    }
  }
  double score = 0;
  if (is_prefix) {
    double c_term_shift = prms[seq_end+1] - prms[prms.size()-1];
    score = compScore(ms_ptr_vec, proteo_ptr, 0, c_term_shift);
  }
  return ZpFastMatchPtr(new ZeroPtmFastMatch(proteo_ptr, score, 0, seq_end));
}

ZpFastMatchPtr computeSuffixMatch(const ExtendMsPtrVec &ms_ptr_vec, 
                                  ProteoformPtr proteo_ptr, double ppo) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error = header_ptr->getErrorTolerance(ppo);
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();
  double diff = prms[prms.size()-1] - res_sum_mass;

  bool is_suffix = false;
  int seq_bgn = 0;
  for (size_t i = 1; i < prms.size(); i++) {
    if (std::abs(prms[i] - diff) <= max_error) {
      is_suffix = true;
      seq_bgn = i;
      break;
    } else {
      if (prms[i] > diff) {
        break;
      }
    }
  }
  double score = 0;
  if (is_suffix) {
    double n_term_shift = - prms[seq_bgn];
    score = compScore(ms_ptr_vec, proteo_ptr, n_term_shift, 0);
  }
  int seq_end = proteo_ptr->getResSeqPtr()->getLen() - 1;
  return ZpFastMatchPtr(
      new ZeroPtmFastMatch(proteo_ptr, score, seq_bgn, seq_end));
}

ZpFastMatchPtr computeInternalMatch(const ExtendMsPtrVec &ms_ptr_vec, 
                                    ProteoformPtr proteo_ptr, double ppo) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error = header_ptr->getErrorTolerance(ppo);
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();

  ActivationPtr activation = header_ptr->getActivationPtr();
  IonTypePtr n_ion_type_ptr = activation->getNIonTypePtr();
  std::vector<double> n_masses 
      = proteo_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type_ptr);
  IonTypePtr c_ion_type_ptr = activation->getCIonTypePtr();
  std::vector<double> c_masses 
      = proteo_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type_ptr);

  double best_score = 0;
  int best_bgn = -1;
  int best_end = -1;

  size_t mass_bgn = 1;
  size_t mass_end = 1;
  while (mass_end < prms.size()-1 && mass_bgn < prms.size()-1) {
    double diff = prms[mass_end] - prms[mass_bgn] - res_sum_mass;
    if (std::abs(diff) <= max_error) {
      double n_term_shift = -prms[mass_bgn];
      double c_term_shift = prms[mass_end] - prms[prms.size() -1];
      double cur_score = 0;
      for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
        cur_score += compDiagScr(ms_ptr_vec[i], n_masses, n_term_shift);
        cur_score += compDiagScr(ms_ptr_vec[i], c_masses, c_term_shift);
      }

      if (cur_score > best_score) {
        best_score = cur_score;
        best_bgn = mass_bgn;
        best_end = mass_end;
      }
    }
    if (diff < 0) {
      mass_end++;
    } else {
      mass_bgn++;
    }
  }
  return ZpFastMatchPtr(new ZeroPtmFastMatch(proteo_ptr, best_score, 
                                             best_bgn, best_end-1));
}

ZpFastMatchPtrVec ZeroPtmFastMatch::filter(AlignTypePtr align_type_ptr,
                                           const ExtendMsPtrVec &ms_ptr_vec,
                                           const ProteoformPtrVec &proteo_ptrs,
                                           int report_num, double ppo) {
  
  ZpFastMatchPtrVec match_vec;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    if (align_type_ptr == AlignType::COMPLETE) { 
        match_vec.push_back(computeCompMatch(ms_ptr_vec, proteo_ptrs[i], ppo));
    }
    else if (align_type_ptr == AlignType::PREFIX) { 
        match_vec.push_back(computePrefixMatch(ms_ptr_vec, proteo_ptrs[i], ppo));
    }
    else if (align_type_ptr == AlignType::SUFFIX) { 
        match_vec.push_back(computeSuffixMatch(ms_ptr_vec, proteo_ptrs[i], ppo));
    }
    else if (align_type_ptr == AlignType::INTERNAL) { 
        match_vec.push_back(computeInternalMatch(ms_ptr_vec, proteo_ptrs[i], ppo));
    }
  }

  // sort 
  std::sort(match_vec.begin(), match_vec.end(), ZeroPtmFastMatch::cmpScoreDec);
  //LOG_DEBUG("sort  finished BEST SCORE " << match_vec[0]->getScore());

  size_t num = report_num;
  if (num > proteo_ptrs.size()) {
    num = proteo_ptrs.size();
  }
  ZpFastMatchPtrVec report_vec;
  for (size_t i = 0; i < num; i++) {
    if (match_vec[i]->getScore() > 0) {
      report_vec.push_back(match_vec[i]);
    } else {
      break;
    }
  }
  return report_vec;
}

}
