// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "search/zeroptmsearch/zero_ptm_fast_search.hpp"

#include <algorithm>
#include <cmath>

namespace toppic {

namespace zero_ptm_fast_search {

/*
 * In the computing of diagonal score in fast filtering, we use a
 * large error tolerance for N-terminal fragment masses
 */
double compDiagScr(const ExtendMsPtr& ms_ptr, const std::vector<double>& masses,
                   double center) {
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

double compScore(const ExtendMsPtr& ms_ptr, const ProteoformPtr& proteo_ptr,
                 double n_shift, double c_shift) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  ActivationPtr activation = header_ptr->getActivationPtr();
  IonTypePtr n_ion_type_ptr = activation->getNIonTypePtr();
  std::vector<double> masses =
      proteo_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type_ptr);
  // the N-terminal fragments of a labeled run (toppic -N) carry the label
  double score = compDiagScr(
      ms_ptr, masses, n_shift + header_ptr->getNTermLabelMass());

  IonTypePtr c_ion_type_ptr = activation->getCIonTypePtr();
  masses = proteo_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type_ptr);
  score += compDiagScr(ms_ptr, masses, c_shift);
  return score;
}

double compScore(const ExtendMsPtrVec& ms_ptr_vec,
                 const ProteoformPtr& proteo_ptr, double n_shift,
                 double c_shift) {
  double score = 0;
  for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
    score = score + compScore(ms_ptr_vec[i], proteo_ptr, n_shift, c_shift);
  }
  return score;
}

ZpFastMatchPtr computeCompMatch(const ExtendMsPtrVec& ms_ptr_vec,
                                const ProteoformPtr& proteo_ptr,
                                const PeakTolerancePtr& prec_tole_ptr) {
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error =
      prec_tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass());
  double res_sum_mass = header_ptr->getFirstPrecMonoMassMinusWater();
  double prot_mass = proteo_ptr->getResSeqPtr()->getResMassSum();
  double error = std::abs(res_sum_mass - prot_mass);
  double score = 0;
  if (error <= max_error) {
    score = compScore(ms_ptr_vec, proteo_ptr, 0, 0);
  }
  int end = proteo_ptr->getResSeqPtr()->getLen() - 1;
  return std::make_shared<ZeroPtmFastMatch>(proteo_ptr, score, 0, end);
}

ZpFastMatchPtr computePrefixMatch(const ExtendMsPtrVec& ms_ptr_vec,
                                  const ProteoformPtr& proteo_ptr,
                                  const PeakTolerancePtr& prec_tole_ptr) {
  /* check if there is a matched prefix */
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error =
      prec_tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass());
  double res_sum_mass = header_ptr->getFirstPrecMonoMassMinusWater();

  bool is_prefix = false;
  int seq_end = 0;
  for (size_t i = 0; i < prms.size() - 1; i++) {
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
    double c_term_shift = prms[seq_end + 1] - prms[prms.size() - 1];
    score = compScore(ms_ptr_vec, proteo_ptr, 0, c_term_shift);
  }
  return std::make_shared<ZeroPtmFastMatch>(proteo_ptr, score, 0, seq_end);
}

ZpFastMatchPtr computeSuffixMatch(const ExtendMsPtrVec& ms_ptr_vec,
                                  const ProteoformPtr& proteo_ptr,
                                  const PeakTolerancePtr& prec_tole_ptr) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error =
      prec_tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass());
  double res_sum_mass = header_ptr->getFirstPrecMonoMassMinusWater();
  double diff = prms[prms.size() - 1] - res_sum_mass;

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
    double n_term_shift = -prms[seq_bgn];
    score = compScore(ms_ptr_vec, proteo_ptr, n_term_shift, 0);
  }
  int seq_end = proteo_ptr->getResSeqPtr()->getLen() - 1;
  return std::make_shared<ZeroPtmFastMatch>(proteo_ptr, score, seq_bgn,
                                            seq_end);
}

ZpFastMatchPtr computeInternalMatch(const ExtendMsPtrVec& ms_ptr_vec,
                                    const ProteoformPtr& proteo_ptr,
                                    const PeakTolerancePtr& prec_tole_ptr) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error =
      prec_tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass());
  double res_sum_mass = header_ptr->getFirstPrecMonoMassMinusWater();

  ActivationPtr activation = header_ptr->getActivationPtr();
  IonTypePtr n_ion_type_ptr = activation->getNIonTypePtr();
  std::vector<double> n_masses =
      proteo_ptr->getBpSpecPtr()->getBreakPointMasses(n_ion_type_ptr);
  IonTypePtr c_ion_type_ptr = activation->getCIonTypePtr();
  std::vector<double> c_masses =
      proteo_ptr->getBpSpecPtr()->getBreakPointMasses(c_ion_type_ptr);

  double best_score = 0;
  int best_bgn = -1;
  int best_end = -1;

  size_t mass_bgn = 1;
  size_t mass_end = 1;
  while (mass_end < prms.size() - 1 && mass_bgn < prms.size() - 1) {
    double diff = prms[mass_end] - prms[mass_bgn] - res_sum_mass;
    if (std::abs(diff) <= max_error) {
      double n_term_shift = -prms[mass_bgn];
      double c_term_shift = prms[mass_end] - prms[prms.size() - 1];
      double cur_score = 0;
      for (size_t i = 0; i < ms_ptr_vec.size(); i++) {
        double n_term_label_mass =
            ms_ptr_vec[i]->getMsHeaderPtr()->getNTermLabelMass();
        cur_score += compDiagScr(ms_ptr_vec[i], n_masses,
                                 n_term_shift + n_term_label_mass);
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
  return std::make_shared<ZeroPtmFastMatch>(proteo_ptr, best_score, best_bgn,
                                            best_end - 1);
}

// The match of a proteoform on the range [bgn, end] of its residues: the
// range the type allows (as the searches above find it: PREFIX not to the
// last residue, SUFFIX not from the first, INTERNAL neither), its mass within
// the precursor error tolerance, and the score of its fragment masses; a
// score of 0 otherwise.
ZpFastMatchPtr computeRangeMatch(const ProteoformTypePtr& align_type_ptr,
                                 const ExtendMsPtrVec& ms_ptr_vec,
                                 const ProteoformPtr& proteo_ptr, int bgn,
                                 int end,
                                 const PeakTolerancePtr& prec_tole_ptr) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  int len = static_cast<int>(prms.size()) - 1;
  MsHeaderPtr header_ptr = ms_ptr_vec[0]->getMsHeaderPtr();
  double max_error =
      prec_tole_ptr->compPrecErrorTole(header_ptr->getFirstPrecMonoMass());
  double res_sum_mass = header_ptr->getFirstPrecMonoMassMinusWater();
  bool allowed = bgn >= 0 && bgn <= end && end < len;
  if (align_type_ptr == ProteoformType::COMPLETE) {
    allowed = allowed && bgn == 0 && end == len - 1;
  } else if (align_type_ptr == ProteoformType::PREFIX) {
    allowed = allowed && bgn == 0 && end < len - 1;
  } else if (align_type_ptr == ProteoformType::SUFFIX) {
    allowed = allowed && bgn > 0 && end == len - 1;
  } else {
    allowed = allowed && bgn > 0 && end < len - 1;
  }
  double score = 0;
  if (allowed &&
      std::abs(prms[end + 1] - prms[bgn] - res_sum_mass) <= max_error) {
    double n_term_shift = -prms[bgn];
    double c_term_shift = prms[end + 1] - prms[len];
    score = compScore(ms_ptr_vec, proteo_ptr, n_term_shift, c_term_shift);
  }
  return std::make_shared<ZeroPtmFastMatch>(proteo_ptr, score, bgn, end);
}

// The best report_num of the matches with a positive score.
ZpFastMatchPtrVec selectTopMatches(ZpFastMatchPtrVec& match_vec,
                                   int report_num) {
  std::sort(match_vec.begin(), match_vec.end(),
            ZeroPtmFastMatch::cmpScoreDecProtNameInc);

  size_t num = std::min(static_cast<size_t>(report_num), match_vec.size());
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

ZpFastMatchPtr computeMatch(const ProteoformTypePtr& align_type_ptr,
                            const ExtendMsPtrVec& ms_ptr_vec,
                            const ProteoformPtr& proteo_ptr,
                            const PeakTolerancePtr& prec_tole_ptr) {
  if (align_type_ptr == ProteoformType::COMPLETE) {
    return computeCompMatch(ms_ptr_vec, proteo_ptr, prec_tole_ptr);
  } else if (align_type_ptr == ProteoformType::PREFIX) {
    return computePrefixMatch(ms_ptr_vec, proteo_ptr, prec_tole_ptr);
  } else if (align_type_ptr == ProteoformType::SUFFIX) {
    return computeSuffixMatch(ms_ptr_vec, proteo_ptr, prec_tole_ptr);
  }
  return computeInternalMatch(ms_ptr_vec, proteo_ptr, prec_tole_ptr);
}

ZpFastMatchPtrVec filterRanges(const ProteoformTypePtr& align_type_ptr,
                               const ExtendMsPtrVec& ms_ptr_vec,
                               const ProteoformPtrVec& proteo_ptrs,
                               const std::vector<std::pair<int, int>>& ranges,
                               int report_num,
                               const PeakTolerancePtr& prec_tole_ptr) {
  ZpFastMatchPtrVec match_vec;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    if (ranges[i].first < 0 || ranges[i].second < 0) {
      match_vec.push_back(computeMatch(align_type_ptr, ms_ptr_vec,
                                       proteo_ptrs[i], prec_tole_ptr));
    } else {
      match_vec.push_back(computeRangeMatch(align_type_ptr, ms_ptr_vec,
                                            proteo_ptrs[i], ranges[i].first,
                                            ranges[i].second, prec_tole_ptr));
    }
  }
  return selectTopMatches(match_vec, report_num);
}

ZpFastMatchPtrVec filter(const ProteoformTypePtr& align_type_ptr,
                         const ExtendMsPtrVec& ms_ptr_vec,
                         const ProteoformPtrVec& proteo_ptrs, int report_num,
                         const PeakTolerancePtr& prec_tole_ptr) {
  ZpFastMatchPtrVec match_vec;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    if (align_type_ptr == ProteoformType::COMPLETE) {
      match_vec.push_back(
          computeCompMatch(ms_ptr_vec, proteo_ptrs[i], prec_tole_ptr));
    } else if (align_type_ptr == ProteoformType::PREFIX) {
      match_vec.push_back(
          computePrefixMatch(ms_ptr_vec, proteo_ptrs[i], prec_tole_ptr));
    } else if (align_type_ptr == ProteoformType::SUFFIX) {
      match_vec.push_back(
          computeSuffixMatch(ms_ptr_vec, proteo_ptrs[i], prec_tole_ptr));
    } else if (align_type_ptr == ProteoformType::INTERNAL) {
      match_vec.push_back(
          computeInternalMatch(ms_ptr_vec, proteo_ptrs[i], prec_tole_ptr));
    }
  }

  // sort
  std::sort(match_vec.begin(), match_vec.end(),
            ZeroPtmFastMatch::cmpScoreDecProtNameInc);

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

}  // namespace zero_ptm_fast_search

}  // namespace toppic
