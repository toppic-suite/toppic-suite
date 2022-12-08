//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <algorithm>

#include "prsm/prsm_algo.hpp"
#include "search/diag/diag_pair_util.hpp"

namespace toppic {

namespace diag_pair_util {

inline DiagPairPtrVec compDiagPair(const PrmPeakPtrVec &prm_peaks,
                                   int group_spec_num, const std::vector<double> &seq_masses,
                                   DiagHeaderPtr header_ptr) {

  std::vector<std::vector<double>> n_term_scores(group_spec_num, 
                                                 std::vector<double>(seq_masses.size(), 0));
  std::vector<std::vector<double>> c_term_scores(group_spec_num, 
                                                 std::vector<double>(seq_masses.size(), 0));
  std::vector<double> top_scores(seq_masses.size(), 0);
  std::vector<int> top_positions(seq_masses.size(), -1);

  // i starts from 0 and ends at size - 1 to include the first prm 0 and the
  // last prm precursor_mass - water_mass
  size_t i = 0;
  size_t j = 0;
  double n_term_shift = header_ptr->getProtNTermShift();
  std::vector<double> real_masses;
  for (size_t k = 0; k < prm_peaks.size(); k++) {
    real_masses.push_back(prm_peaks[k]->getPosition());
  }

  while (i < prm_peaks.size() && j < seq_masses.size()) {
    PrmPeakPtr peak = prm_peaks[i];
    BasePeakTypePtr type_ptr = peak->getBaseTypePtr();
    int spec_id = peak->getSpectrumId();
    double error = 0;
    if (header_ptr->isNStrict() && !header_ptr->isCStrict()) {
      error = peak->getNStrictCRelaxTolerance();
    } else {
      error = peak->getNRelaxCStrictTolerance();
    }
    double deviation = peak->getPosition()-seq_masses[j] - n_term_shift;
    if (std::abs(deviation) <= error) {
      double peak_score = peak->getScore();
      if (type_ptr == BasePeakType::ORIGINAL) {
        if (n_term_scores[spec_id][j] < peak_score) {
          n_term_scores[spec_id][j] = peak_score;
        }
      } else {
        if (c_term_scores[spec_id][j] < peak_score) {
          c_term_scores[spec_id][j] = peak_score;
        }
      }
      // update top position
      if (top_scores[j] < peak_score) {
        top_scores[j] = peak_score;
        top_positions[j] = i;
      }
    }
    if (prsm_algo::increaseIJ(i, j, deviation, peak->getNRelaxCStrictTolerance(),
                              real_masses, seq_masses)) {
      i++;
    } else {
      j++;
    }
  }
  std::vector<double> sum_scores(seq_masses.size(), 0);
  for (size_t p = 0; p < seq_masses.size(); p++) {
    for (int m = 0; m < group_spec_num; m++) {
      sum_scores[p] += n_term_scores[m][p];
      sum_scores[p] += c_term_scores[m][p];
    }
  }

  // add pairs
  DiagPairPtrVec  pair_list;
  for (j = 0; j < seq_masses.size(); j++) {
    int pos = top_positions[j];
    if (pos >= 0) {
      double diff = prm_peaks[pos]->getPosition() - seq_masses[j];
      double score = sum_scores[j];
      DiagPairPtr diag_pair_ptr
        = std::make_shared<DiagPair>(pos, j, score, pair_list.size(), diff);
      pair_list.push_back(diag_pair_ptr);
    }
  }
  // compare y (seq_mass) first, then x (peak position)
  std::sort(pair_list.begin(), pair_list.end(), Pair::cmpPosInc);
  return pair_list;
}

DiagonalPtr getDiagonalPtrWithoutEmptyList(DiagHeaderPtr header_ptr,
                                           const PrmPeakPtrVec &prm_peaks,
                                           int group_spec_num,
                                           ProteoformPtr proteo_ptr) {
  BpSpecPtr bp_spec_ptr = proteo_ptr->getBpSpecPtr();

  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  DiagPairPtrVec diag_pair_list = compDiagPair(prm_peaks, group_spec_num, 
                                               prm_masses, header_ptr);
  if (diag_pair_list.size() > 0) {
    DiagonalPtr diagonal_ptr
      = std::make_shared<Diagonal>(header_ptr, diag_pair_list);
    for (size_t i = 0; i < diag_pair_list.size(); i++) {
      diag_pair_list[i]->setDiagonalPtr(diagonal_ptr);
    }
    return diagonal_ptr;
  }
  return nullptr;
}

DiagonalPtr getDiagonalPtrWithEmptyList(DiagHeaderPtr header_ptr,
                                        const PrmPeakPtrVec &prm_peaks,
                                        int group_spec_num,
                                        ProteoformPtr proteo_ptr) {
  BpSpecPtr bp_spec_ptr = proteo_ptr->getBpSpecPtr();

  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  DiagPairPtrVec diag_pair_list = compDiagPair(prm_peaks, group_spec_num, 
                                               prm_masses, header_ptr);
  DiagonalPtr diagonal_ptr
    = std::make_shared<Diagonal>(header_ptr, diag_pair_list);
  for (size_t i = 0; i < diag_pair_list.size(); i++) {
    diag_pair_list[i]->setDiagonalPtr(diagonal_ptr);
  }
  return diagonal_ptr;
}


DiagonalPtrVec geneDiagonalsWithoutEmptyList(const DiagHeaderPtrVec& header_ptr_vec,
                                             const PrmPeakPtrVec &prm_peaks,
                                             int group_spec_num, ProteoformPtr proteo_ptr) {
  DiagonalPtrVec diagonal_list;
  for (size_t i = 0; i < header_ptr_vec.size(); i++) {
    DiagonalPtr diagonal_ptr = getDiagonalPtrWithoutEmptyList(header_ptr_vec[i], prm_peaks,
                                                              group_spec_num, proteo_ptr);
    if (diagonal_ptr != nullptr) {
      diagonal_list.push_back(diagonal_ptr);
    }
  }
  // important set id for headers
  for (size_t i = 0; i < diagonal_list.size(); i++) {
    diagonal_list[i]->getHeader()->setId(i);
  }
  return diagonal_list;
}

DiagonalPtrVec geneDiagonalsWithEmptyList(const DiagHeaderPtrVec& header_ptr_vec,
                                          const PrmPeakPtrVec &prm_peaks,
                                          int group_spec_num, ProteoformPtr proteo_ptr) {
  DiagonalPtrVec diagonal_list;
  for (size_t i = 0; i < header_ptr_vec.size(); i++) {
    DiagonalPtr diagonal_ptr = getDiagonalPtrWithEmptyList(header_ptr_vec[i], prm_peaks,
                                                           group_spec_num, proteo_ptr);
    diagonal_list.push_back(diagonal_ptr);
  }
  // important set id for headers
  for (size_t i = 0; i < diagonal_list.size(); i++) {
    diagonal_list[i]->getHeader()->setId(i);
  }
  return diagonal_list;
}


}

} /* namespace toppic */
