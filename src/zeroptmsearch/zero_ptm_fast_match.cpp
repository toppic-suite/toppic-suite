#include <base/logger.hpp>
#include "base/base_data.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"

namespace prot {

ZpFastMatchPtrVec zeroPtmFastFilter(int semi_align_type,
                                    ExtendMsPtr ms_ptr,
                                    ProteoformPtrVec &form_ptr_vec,
                                    int report_num) {
  ZpFastMatchPtrVec match_vec;
  for (unsigned int i = 0; i < form_ptr_vec.size(); i++) {
    switch (semi_align_type) {
      case SEMI_ALIGN_TYPE_COMPLETE: 
        match_vec.push_back(computeCompMatch(ms_ptr, form_ptr_vec[i]));
        break;
      case SEMI_ALIGN_TYPE_PREFIX:
        match_vec.push_back(computePrefixMatch(ms_ptr, form_ptr_vec[i]));
        break;
      case SEMI_ALIGN_TYPE_SUFFIX:
        match_vec.push_back(computeSuffixMatch(ms_ptr, form_ptr_vec[i]));
        break;
      case SEMI_ALIGN_TYPE_INTERNAL:
        match_vec.push_back(computeInternalMatch(ms_ptr, form_ptr_vec[i]));
        break;
    }
  }

  /* sort */
  std::sort(match_vec.begin(), match_vec.end(), compareZeroPtmFastMatchDown);
  LOG_DEBUG("sort  finished BEST SCORE " << match_vec[0]->getScore());

  unsigned int num = report_num;
  if (num > form_ptr_vec.size()) {
    num = form_ptr_vec.size();
  }
  ZpFastMatchPtrVec report_vec;
  for (unsigned int i = 0; i < num; i++) {
    if (match_vec[i]->getScore() > 0) {
      report_vec.push_back(match_vec[i]);
    } else {
      break;
    }
  }
  return report_vec;
}

ZpFastMatchPtr computeCompMatch(ExtendMsPtr ms_ptr, ProteoformPtr form_ptr) {
  MsHeaderPtr header_ptr = ms_ptr->getHeaderPtr();
  double max_error = header_ptr->getErrorTolerance();
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();
  double prot_mass = form_ptr->getResSeqPtr()->getResMassSum();
  double error = abs(res_sum_mass - prot_mass);
  LOG_TRACE("Protein mass " << prot_mass 
            << " precursor mass " << res_sum_mass 
            << " proteoform name " << form_ptr->getName()
            << " error " << error << " error tolerance " << max_error);
  double score = 0;
  if (error <= max_error) {
    ActivationPtr activation = header_ptr->getActivationPtr();
    IonTypePtr n_ion_type = activation->getNIonType();
    std::vector<double> masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type);
    score = compDiagScr(ms_ptr, masses, 0);

    IonTypePtr c_ion_type = activation->getCIonType();
    masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type);
    score += compDiagScr(ms_ptr, masses, 0);
    LOG_TRACE("score " << score);
  }
  return ZpFastMatchPtr(
      new ZeroPtmFastMatch(form_ptr, score, 0, form_ptr->getResSeqPtr()->getLen() - 1));
}

ZpFastMatchPtr computePrefixMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr) {
  /* check if there is a matched prefix */
  std::vector<double> prms = form_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr->getHeaderPtr();
  double max_error = header_ptr->getErrorTolerance();
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();

  bool is_prefix = false;
  int seq_end = 0;
  for (unsigned int i = 0; i < prms.size() - 1; i++) {
    if (abs(res_sum_mass - prms[i]) <= max_error) {
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
    ActivationPtr activation = header_ptr->getActivationPtr();
    IonTypePtr n_ion_type = activation->getNIonType();
    std::vector<double> masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type);
    score = compDiagScr(ms_ptr, masses, 0);

    IonTypePtr c_ion_type = activation->getCIonType();
    masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type);
    score += compDiagScr(ms_ptr, masses, c_term_shift);
  }
  return ZpFastMatchPtr(new ZeroPtmFastMatch(form_ptr, score, 0, seq_end));
}

ZpFastMatchPtr computeSuffixMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr) {
  std::vector<double> prms = form_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr->getHeaderPtr();
  double max_error = header_ptr->getErrorTolerance();
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();
  double diff = prms[prms.size()-1] - res_sum_mass;

  bool is_suffix = false;
  int seq_bgn = 0;
  for (unsigned int i = 1; i < prms.size(); i++) {
    if (abs(prms[i] - diff) <= max_error) {
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
    ActivationPtr activation = header_ptr->getActivationPtr();
    IonTypePtr n_ion_type = activation->getNIonType();
    std::vector<double> masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type);
    score = compDiagScr(ms_ptr, masses, n_term_shift);

    IonTypePtr c_ion_type = activation->getCIonType();
    masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type);
    score += compDiagScr(ms_ptr, masses, 0);
  }
  return ZpFastMatchPtr(
      new ZeroPtmFastMatch(form_ptr, score, seq_bgn, form_ptr->getResSeqPtr()->getLen() - 1));
}

ZpFastMatchPtr computeInternalMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr) {
  std::vector<double> prms = form_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr->getHeaderPtr();
  double max_error = header_ptr->getErrorTolerance();
  double res_sum_mass = header_ptr->getPrecMonoMassMinusWater();

  ActivationPtr activation = header_ptr->getActivationPtr();
  IonTypePtr n_ion_type = activation->getNIonType();
  std::vector<double> n_masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type);
  IonTypePtr c_ion_type = activation->getCIonType();
  std::vector<double> c_masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type);

  double best_score = 0;
  int best_bgn = -1;
  int best_end = -1;

  unsigned int mass_bgn = 1;
  unsigned int mass_end = 1;
  while (mass_end < prms.size()-1 && mass_bgn < prms.size()-1) {
    double diff = prms[mass_end] - prms[mass_bgn] - res_sum_mass;
    if (abs(diff) <= max_error) {
      double n_term_shift = -prms[mass_bgn];
      double cur_score = compDiagScr(ms_ptr, n_masses, n_term_shift);
      double c_term_shift = prms[mass_end] - prms[prms.size() -1];
      cur_score += compDiagScr(ms_ptr, c_masses, c_term_shift);

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
  return ZpFastMatchPtr(new ZeroPtmFastMatch(form_ptr, best_score, best_bgn, best_end-1));
}


double compDiagScr(ExtendMsPtr ms_ptr,
                  std::vector<double> &masses, double center) {
  unsigned int i = 0;
  unsigned int j = 0;
  double s = 0;
  while (i < ms_ptr->size() && j < masses.size()) {
    ExtendPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    double distance = peak_ptr->getMonoMass() - masses[j];
    //LOG_DEBUG("peak " << peak_ptr->getMonoMass() << " score " << peak_ptr->getScore());
    if (abs(center - distance) <= peak_ptr->getOrigTolerance()) {
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

}
