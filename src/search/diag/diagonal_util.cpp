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

#include "search/diag/diagonal_util.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "common/base/neutral_loss_base.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/extend_ms_util.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/theo_peak_util.hpp"

namespace toppic {

namespace diagonal_util {

// header_ptr: the header of the extended spectrum the peaks are matched
// with; its activation gives the ion types and its N-terminal label (toppic
// -N) is carried by the N-terminal ions
inline TheoPeakPtrVec getDiagonalTheoPeak(const ProteoformPtr& proteo_ptr,
                                          const MsHeaderPtr& header_ptr,
                                          const DiagHeaderPtrVec& header_ptrs,
                                          int i, double min_mass) {
  ActivationPtr activation_ptr = header_ptr->getActivationPtr();
  DiagHeaderPtr first_header_ptr = header_ptrs[0];
  DiagHeaderPtr last_header_ptr = header_ptrs[header_ptrs.size() - 1];
  int first_res_pos = first_header_ptr->getTruncFirstResPos();
  int last_res_pos = last_header_ptr->getTruncLastResPos();
  ResSeqPtr subseq_ptr =
      proteo_ptr->getResSeqPtr()->getSubResidueSeq(first_res_pos, last_res_pos);
  BpSpecPtr pep_ptr = std::make_shared<BpSpec>(subseq_ptr);
  double pep_n_term_shift = header_ptrs[i]->getProtNTermShift() -
                            first_header_ptr->getProtNTermShift() +
                            first_header_ptr->getPepNTermShift() +
                            header_ptr->getNTermLabelMass();
  double pep_c_term_shift = header_ptrs[i]->getProtCTermShift() -
                            last_header_ptr->getProtCTermShift() +
                            last_header_ptr->getPepCTermShift();
  double max_mass =
      subseq_ptr->getSeqMass() + pep_n_term_shift + pep_c_term_shift - min_mass;
  return theo_peak_util::geneTheoPeak(
      pep_ptr, activation_ptr, NeutralLossBase::getNeutralLossPtr_NONE(),
      pep_n_term_shift, pep_c_term_shift,
      header_ptrs[i]->getMatchFirstBpPos() - first_res_pos,
      header_ptrs[i]->getMatchLastBpPos() - first_res_pos, min_mass, max_mass);
}

namespace {

// a precursor-dependent theoretical ion and an extended peak within the
// peak's fragment tolerance plus the precursor tolerance: the difference
// theoretical - observed and the peak's fragment tolerance
struct DepIonPair {
  int ion_id;
  double diff;
  double frag_tole;
};

// the precursor-dependent ions of the headers (the C-terminal ions of every
// diagonal but the last and the N-terminal ions of the last: the ions whose
// theoretical masses move with the precursor mass) paired with the peaks of
// each spectrum. A single diagonal (TopMG's graph alignment without shifts)
// is anchored at its N-terminus, so its C-terminal ions are the dependent
// ones. ion_num gets the number of dependent ions.
std::vector<DepIonPair> getDepIonPairs(const ProteoformPtr& proteo_ptr,
                                       const ExtendMsPtrVec& ms_three_ptr_vec,
                                       const DiagHeaderPtrVec& header_ptrs,
                                       double prec_tole, double min_mass,
                                       int& ion_num) {
  std::vector<DepIonPair> pairs;
  ion_num = 0;
  int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
  size_t last = header_ptrs.size() - 1;
  // the diagonal whose N-terminal ions are the dependent ones; none when
  // there is one diagonal
  size_t n_term_idx = last > 0 ? last : header_ptrs.size();
  for (size_t j = 0; j < ms_three_ptr_vec.size(); j++) {
    const ExtendMsPtr& ms_ptr = ms_three_ptr_vec[j];
    std::vector<double> ms_masses = extend_ms_util::getExtendMassVec(ms_ptr);
    double max_frag_tole = 0;
    for (size_t k = 0; k < ms_ptr->size(); k++) {
      max_frag_tole =
          std::max(max_frag_tole, ms_ptr->getPeakPtr(k)->getOrigTolerance());
    }
    for (size_t i = 0; i < header_ptrs.size(); i++) {
      TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
          proteo_ptr, ms_ptr->getMsHeaderPtr(), header_ptrs, i, min_mass);
      int bgn = header_ptrs[i]->getMatchFirstBpPos() - first_res_pos;
      int end = header_ptrs[i]->getMatchLastBpPos() - first_res_pos;
      for (size_t t = 0; t < theo_peak_ptrs.size(); t++) {
        const IonPtr& ion_ptr = theo_peak_ptrs[t]->getIonPtr();
        if (ion_ptr->getIonTypePtr()->isNTerm() != (i == n_term_idx) ||
            ion_ptr->getPos() < bgn || ion_ptr->getPos() > end) {
          continue;
        }
        double theo_mass = theo_peak_ptrs[t]->getPosition();
        auto it = std::lower_bound(ms_masses.begin(), ms_masses.end(),
                                   theo_mass - max_frag_tole - prec_tole);
        for (size_t k = it - ms_masses.begin(); k < ms_masses.size(); k++) {
          if (ms_masses[k] > theo_mass + max_frag_tole + prec_tole) {
            break;
          }
          double frag_tole = ms_ptr->getPeakPtr(k)->getOrigTolerance();
          double diff = theo_mass - ms_masses[k];
          if (std::abs(diff) <= frag_tole + prec_tole) {
            pairs.push_back(DepIonPair{ion_num, diff, frag_tole});
          }
        }
        ion_num++;
      }
    }
  }
  return pairs;
}

// the distance of [lo, hi] to 0
double distToZero(double lo, double hi) {
  if (lo <= 0 && 0 <= hi) {
    return 0;
  }
  return std::min(std::abs(lo), std::abs(hi));
}

}  // namespace

double refinePrecursorByMatchedIons(
    const ProteoformPtr& proteo_ptr, const ExtendMsPtrVec& ms_three_ptr_vec,
    DiagHeaderPtrVec& header_ptrs, const PeakTolerancePtr& prec_tole_ptr,
    double min_mass, double delta_min, double delta_max) {
  double prec_mass =
      ms_three_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
  if (header_ptrs.size() == 0) {
    return prec_mass;
  }
  double tole = prec_tole_ptr->compPrecErrorTole(prec_mass);
  // the delta range [-tole, tole], cut to [delta_min, delta_max]
  double range_lo = std::max(-tole, delta_min);
  double range_hi = std::min(tole, delta_max);
  int ion_num = 0;
  std::vector<DepIonPair> pairs = getDepIonPairs(
      proteo_ptr, ms_three_ptr_vec, header_ptrs, tole, min_mass, ion_num);

  // After the refinement by delta, a pair matches within the fragment
  // tolerance if |diff + delta| <= frag_tole. Sweep the intervals of delta
  // in [-tole, tole] and find the segment covered by the most ions (each
  // ion counted once); among equal segments the one closest to 0.
  struct Event {
    double pos;
    int type;  // 0 start, 1 end
    int ion_id;
  };
  std::vector<Event> events;
  for (size_t p = 0; p < pairs.size(); p++) {
    double lo = std::max(-pairs[p].diff - pairs[p].frag_tole, range_lo);
    double hi = std::min(-pairs[p].diff + pairs[p].frag_tole, range_hi);
    if (lo <= hi) {
      events.push_back(Event{lo, 0, pairs[p].ion_id});
      events.push_back(Event{hi, 1, pairs[p].ion_id});
    }
  }
  std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
    return a.pos < b.pos || (a.pos == b.pos && a.type < b.type);
  });
  std::vector<int> cover(ion_num, 0);
  int count = 0;
  int best_count = 0;
  double best_lo = 0;
  double best_hi = 0;
  for (size_t e = 0; e + 1 < events.size(); e++) {
    if (events[e].type == 0) {
      if (cover[events[e].ion_id]++ == 0) {
        count++;
      }
    } else if (--cover[events[e].ion_id] == 0) {
      count--;
    }
    // the coverage of [events[e].pos, events[e + 1].pos]
    double lo = events[e].pos;
    double hi = events[e + 1].pos;
    if (count > best_count ||
        (count == best_count && count > 0 &&
         distToZero(lo, hi) < distToZero(best_lo, best_hi))) {
      best_count = count;
      best_lo = lo;
      best_hi = hi;
    }
  }
  double delta = 0;
  if (best_count > 0) {
    // the ions matched in the segment, each with its pair closest to the
    // segment's midpoint; delta is the median of their errors, clipped to
    // the delta range that keeps all of them matched. (1e-9 Da: a segment
    // of one point is the end of an interval.)
    double mid = (best_lo + best_hi) / 2;
    std::vector<int> best_pair(ion_num, -1);
    for (size_t p = 0; p < pairs.size(); p++) {
      if (std::abs(pairs[p].diff + mid) > pairs[p].frag_tole + 1e-9) {
        continue;
      }
      int& b = best_pair[pairs[p].ion_id];
      if (b < 0 || std::abs(pairs[p].diff + mid) <
                       std::abs(pairs[b].diff + mid)) {
        b = static_cast<int>(p);
      }
    }
    double lo = range_lo;
    double hi = range_hi;
    std::vector<double> errors;
    for (int i = 0; i < ion_num; i++) {
      if (best_pair[i] >= 0) {
        const DepIonPair& pair = pairs[best_pair[i]];
        lo = std::max(lo, -pair.diff - pair.frag_tole);
        hi = std::min(hi, -pair.diff + pair.frag_tole);
        errors.push_back(-pair.diff);
      }
    }
    std::sort(errors.begin(), errors.end());
    size_t n = errors.size();
    double median = n % 2 == 1 ? errors[n / 2]
                               : (errors[n / 2 - 1] + errors[n / 2]) / 2;
    delta = std::min(std::max(median, lo), hi);
  }

  if (header_ptrs.size() == 1) {
    header_ptrs[0]->changeOnlyCTermShift(delta);
    return prec_mass + delta;
  }
  for (size_t i = 0; i < header_ptrs.size() - 1; i++) {
    header_ptrs[i]->changeOnlyCTermShift(delta);
  }
  header_ptrs[header_ptrs.size() - 1]->changeOnlyNTermShift(delta);
  return prec_mass + delta;
}

double refineVarPtmPrecursorAndHeaderShift(const ProteoformPtr& proteo_ptr,
                                           const DiagHeaderPtrVec& header_ptrs,
                                           double ptm_mass) {
  double refine_prec_mass = proteo_ptr->getMass() + ptm_mass;
  // as VarPtmSlowMatch::geneVarPtmNTermShiftHeaders sets them, with the
  // refined precursor mass in place of the observed one. With the observed
  // one, a precursor 0.02 Da above the proteoform mass put every C-terminal
  // ion below 2 kDa outside the fragment tolerance, and the shift range
  // between two diagonals grew over the ions (toppic_test experiment 120)
  double seq_mass = proteo_ptr->getResSeqPtr()->getSeqMass();
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    double c_shift =
        refine_prec_mass - seq_mass - header_ptrs[i]->getProtNTermShift();
    header_ptrs[i]->initHeader(c_shift, proteo_ptr);
  }
  return refine_prec_mass;
}

DiagHeaderPtrVec refineHeadersBgnEnd(const ProteoformPtr& proteo_ptr,
                                     const ExtendMsPtrVec& ms_three_ptr_vec,
                                     const DiagHeaderPtrVec& header_ptrs,
                                     double min_mass) {
  DiagHeaderPtrVec result_list;
  int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
  int last_res_pos = header_ptrs[header_ptrs.size() - 1]->getTruncLastResPos();
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    int bgn = header_ptrs[i]->getMatchFirstBpPos() - first_res_pos;
    int end = header_ptrs[i]->getMatchLastBpPos() - first_res_pos;
    PeakIonPairPtrVec pair_ptrs;
    for (size_t j = 0; j < ms_three_ptr_vec.size(); j++) {
      TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
          proteo_ptr, ms_three_ptr_vec[j]->getMsHeaderPtr(),
          header_ptrs, i, min_mass);
      PeakIonPairPtrVec cur_pair_ptrs = peak_ion_pair_util::findPairs(
          ms_three_ptr_vec[j], theo_peak_ptrs, bgn, end, 0);
      pair_ptrs.insert(pair_ptrs.end(), cur_pair_ptrs.begin(),
                       cur_pair_ptrs.end());
    }
    if (pair_ptrs.size() < 1) {
      int pair_size = pair_ptrs.size();
      LOG_DEBUG("Empty Segment is found " + std::to_string(pair_size));
      if (i == 0) {
        int new_bgn = first_res_pos;
        int new_end = first_res_pos;
        header_ptrs[i]->setMatchFirstBpPos(new_bgn);
        header_ptrs[i]->setMatchLastBpPos(new_end);
        result_list.push_back(header_ptrs[i]);
      } else if (i == header_ptrs.size() - 1) {
        int new_bgn = last_res_pos + 1;
        int new_end = last_res_pos + 1;
        header_ptrs[i]->setMatchFirstBpPos(new_bgn);
        header_ptrs[i]->setMatchLastBpPos(new_end);
        result_list.push_back(header_ptrs[i]);
      }
    } else {
      int new_bgn = first_res_pos + getNewBgn(pair_ptrs);
      int new_end = first_res_pos + getNewEnd(pair_ptrs);
      if (i == 0) {
        new_bgn = first_res_pos;
      }
      if (i == header_ptrs.size() - 1) {
        new_end = last_res_pos + 1;
      }
      header_ptrs[i]->setMatchFirstBpPos(new_bgn);
      header_ptrs[i]->setMatchLastBpPos(new_end);
      result_list.push_back(header_ptrs[i]);
    }
  }
  return result_list;
}

DiagHeaderPtrVec refineVarPtmHeadersBgnEnd(
    const ProteoformPtr& proteo_ptr, const ExtendMsPtrVec& ms_three_ptr_vec,
    const DiagHeaderPtrVec& header_ptrs, double min_mass) {
  DiagHeaderPtrVec result_list;
  int first_res_pos = header_ptrs[0]->getTruncFirstResPos();
  int last_res_pos = header_ptrs[header_ptrs.size() - 1]->getTruncLastResPos();
  for (size_t i = 0; i < header_ptrs.size(); i++) {
    int bgn = header_ptrs[i]->getMatchFirstBpPos() - first_res_pos;
    int end = header_ptrs[i]->getMatchLastBpPos() - first_res_pos;
    PeakIonPairPtrVec pair_ptrs;
    for (size_t j = 0; j < ms_three_ptr_vec.size(); j++) {
      TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
          proteo_ptr, ms_three_ptr_vec[j]->getMsHeaderPtr(),
          header_ptrs, i, min_mass);
      PeakIonPairPtrVec cur_pair_ptrs = peak_ion_pair_util::findPairs(
          ms_three_ptr_vec[j], theo_peak_ptrs, bgn, end, 0);
      pair_ptrs.insert(pair_ptrs.end(), cur_pair_ptrs.begin(),
                       cur_pair_ptrs.end());
    }
    if (pair_ptrs.size() < 1) {
      int pair_size = pair_ptrs.size();
      LOG_DEBUG("Empty Segment is found " + std::to_string(pair_size));
      if (i == 0) {
        int new_bgn = first_res_pos;
        int new_end = first_res_pos;
        header_ptrs[i]->setMatchFirstBpPos(new_bgn);
        header_ptrs[i]->setMatchLastBpPos(new_end);
        result_list.push_back(header_ptrs[i]);
      } else if (i == header_ptrs.size() - 1) {
        int new_bgn = last_res_pos + 1;
        int new_end = last_res_pos + 1;
        header_ptrs[i]->setMatchFirstBpPos(new_bgn);
        header_ptrs[i]->setMatchLastBpPos(new_end);
        result_list.push_back(header_ptrs[i]);
      } else {
        // add empty header_ptr;
        result_list.push_back(nullptr);
      }
    } else {
      int new_bgn = first_res_pos + getNewBgn(pair_ptrs);
      int new_end = first_res_pos + getNewEnd(pair_ptrs);
      if (i == 0) {
        new_bgn = first_res_pos;
      }
      if (i == header_ptrs.size() - 1) {
        new_end = last_res_pos + 1;
      }
      header_ptrs[i]->setMatchFirstBpPos(new_bgn);
      header_ptrs[i]->setMatchLastBpPos(new_end);
      result_list.push_back(header_ptrs[i]);
    }
  }
  return result_list;
}

// used in graph alignment, each vector contains a set of diagonals that can be
// explained by variable PTMs.
DiagHeaderPtrVec2D refineHeadersBgnEnd(const ProteoformPtr& proteo_ptr,
                                       const ExtendMsPtrVec& ms_three_ptr_vec,
                                       const DiagHeaderPtrVec2D& header_ptrs_2d,
                                       const DiagHeaderPtrVec& header_ptrs_1d,
                                       double min_mass) {
  DiagHeaderPtrVec2D result_list;
  DiagHeaderPtr first_header = header_ptrs_1d[0];
  DiagHeaderPtr last_header = header_ptrs_1d[header_ptrs_1d.size() - 1];
  int first_res_pos = first_header->getTruncFirstResPos();
  int last_res_pos = last_header->getTruncLastResPos();
  int index = 0;
  for (size_t i = 0; i < header_ptrs_2d.size(); i++) {
    DiagHeaderPtrVec cur_vec;
    for (size_t j = 0; j < header_ptrs_2d[i].size(); j++) {
      int bgn = header_ptrs_2d[i][j]->getMatchFirstBpPos() - first_res_pos;
      int end = header_ptrs_2d[i][j]->getMatchLastBpPos() - first_res_pos;
      PeakIonPairPtrVec pair_ptrs;
      for (size_t k = 0; k < ms_three_ptr_vec.size(); k++) {
        TheoPeakPtrVec theo_peak_ptrs = getDiagonalTheoPeak(
            proteo_ptr,
            ms_three_ptr_vec[k]->getMsHeaderPtr(),
            header_ptrs_1d, index, min_mass);
        PeakIonPairPtrVec cur_pair_ptrs = peak_ion_pair_util::findPairs(
            ms_three_ptr_vec[k], theo_peak_ptrs, bgn, end, 0);
        pair_ptrs.insert(pair_ptrs.end(), cur_pair_ptrs.begin(),
                         cur_pair_ptrs.end());
      }
      index++;
      if (pair_ptrs.size() < 1) {
        int pair_size = pair_ptrs.size();
        LOG_DEBUG("Empty Segment is found " + std::to_string(pair_size));
        if (header_ptrs_2d[i][j] == first_header) {
          int new_bgn = first_res_pos;
          int new_end = first_res_pos;
          header_ptrs_2d[i][j]->setMatchFirstBpPos(new_bgn);
          header_ptrs_2d[i][j]->setMatchLastBpPos(new_end);
          cur_vec.push_back(header_ptrs_2d[i][j]);
        } else if (header_ptrs_2d[i][j] == last_header) {
          int new_bgn = last_res_pos + 1;
          int new_end = last_res_pos + 1;
          header_ptrs_2d[i][j]->setMatchFirstBpPos(new_bgn);
          header_ptrs_2d[i][j]->setMatchLastBpPos(new_end);
          cur_vec.push_back(header_ptrs_2d[i][j]);
        }
      } else {
        int new_bgn = first_res_pos + getNewBgn(pair_ptrs);
        int new_end = first_res_pos + getNewEnd(pair_ptrs);
        if (header_ptrs_2d[i][j] == first_header) {
          new_bgn = first_res_pos;
        }
        if (header_ptrs_2d[i][j] == last_header) {
          new_end = last_res_pos + 1;
        }
        header_ptrs_2d[i][j]->setMatchFirstBpPos(new_bgn);
        header_ptrs_2d[i][j]->setMatchLastBpPos(new_end);
        cur_vec.push_back(header_ptrs_2d[i][j]);
      }
    }
    if (cur_vec.size() > 0) {
      result_list.push_back(cur_vec);
    }
  }
  return result_list;
}

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs) {
  int new_bgn = std::numeric_limits<int>::max();
  for (size_t i = 0; i < pair_ptrs.size(); i++) {
    if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos() < new_bgn) {
      new_bgn = pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
    }
  }
  return new_bgn;
}

int getNewEnd(const PeakIonPairPtrVec& pair_ptrs) {
  int new_end = 0;
  for (size_t i = 0; i < pair_ptrs.size(); i++) {
    if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos() > new_end) {
      new_end = pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
    }
  }
  return new_end;
}

}  // namespace diagonal_util

}  // namespace toppic
