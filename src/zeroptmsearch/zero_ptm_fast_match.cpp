
#include <log4cxx/logger.h>

#include "base/base_data.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("ComputeMatch"));

std::vector<ZeroPtmFastMatch> zeroPtmFastFilter(int semi_align_type,
                                                ExtendMsPtr ms_ptr,
                                                ProteoformPtrVec &form_ptr_vec,
                                                int report_num) {
  std::vector<ZeroPtmFastMatch> match_vec;
  for (unsigned int i = 0; i < form_ptr_vec.size(); i++) {
    switch (semi_align_type) {
      case SEMI_ALIGN_TYPE_COMPLETE: 
        match_vec.push_back(computeCompMatch(ms_ptr, form_ptr_vec[i]));
        break;
      case SEMI_ALIGN_TYPE_PREFIX:
        match_vec.push_back(computePrefixMatch(ms_ptr, form_ptr_vec[i]));
        break;
      case SEMI_ALIGN_TYPE_SUFFIX:
        //matches[i] = CompMatch
        //    .computeSuffixMatch(msThree, sequences[i]);
        break;
      case SEMI_ALIGN_TYPE_INTERNAL:
        //matches[i] = CompMatch.computeInternalMatch(msThree,
        //                                            sequences[i]);
        break;
    }
  }

  /* sort */
  std::sort(match_vec.begin(), match_vec.end(), compareZeroPtmFastMatchDown);

  unsigned int num = report_num;
  if (num > form_ptr_vec.size()) {
    num = form_ptr_vec.size();
  }
  std::vector<ZeroPtmFastMatch> report_vec;
  for (unsigned int i = 0; i < num; i++) {
    if (match_vec[i].getScore() > 0) {
      report_vec.push_back(match_vec[i]);
    } else {
      break;
    }
  }
  return report_vec;
}

ZeroPtmFastMatch computeCompMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr) {
  //to do: double maxError = ms_ptr->getHeaderPtr()->getHeader().getErrorTolerance();
  MsHeaderPtr header_ptr = ms_ptr->getHeaderPtr();
  double max_error = 0.01;
  double prec_mass = header_ptr->getPrecMonoMassMinusWater();
  double prot_mass = form_ptr->getResSeqPtr()->getResMassSum();
  double error = abs(prec_mass - prot_mass);
  LOG4CXX_TRACE(logger, "Protein mass " << prot_mass 
                << " precursor mass " << prec_mass 
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
  }
  return ZeroPtmFastMatch(form_ptr, score, 0, form_ptr->getResSeqPtr()->getLen() - 1);
}

ZeroPtmFastMatch computePrefixMatch(
    ExtendMsPtr ms_ptr, ProteoformPtr form_ptr) {
  /* check if there is a matched prefix */
  bool is_prefix = false;
  std::vector<double> prms = form_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr->getHeaderPtr();
  // to do double maxError = msThree.getHeader().getErrorTolerance();
  double max_error = 0.01;
  double prec_mass = header_ptr->getPrecMonoMassMinusWater();
  int seq_end = 0;
  double c_term_shift = 0;
  for (unsigned int i = 0; i < prms.size() - 1; i++) {
    if (abs(prec_mass - prms[i]) <= max_error) {
      is_prefix = true;
      seq_end = i - 1;
      c_term_shift = prms[i] - prms[prms.size() -1];
      break;
    } else {
      if (prms[i] > prec_mass) {
        break;
      }
    }
  }
  double score = 0;
  if (is_prefix) {
    ActivationPtr activation = header_ptr->getActivationPtr();
    IonTypePtr n_ion_type = activation->getNIonType();
    std::vector<double> masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type);
    score = compDiagScr(ms_ptr, masses, 0);

    IonTypePtr c_ion_type = activation->getCIonType();
    masses = form_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type);
    score += compDiagScr(ms_ptr, masses, c_term_shift);
  }
  return ZeroPtmFastMatch(form_ptr, score, 0, seq_end);
}


double compDiagScr(ExtendMsPtr ms_ptr,
                  std::vector<double> &masses, double center) {
  unsigned int i = 0;
  unsigned int j = 0;
  double s = 0;
  while (i < ms_ptr->size() && j < masses.size()) {
    ExtendPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    double distance = peak_ptr->getMonoMass() - masses[j];
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
