//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include <utility>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "base/amino_acid_base.hpp"
#include "base/activation_base.hpp"
#include "base/mod_util.hpp"
#include "base/residue_base.hpp"
#include "base/residue_util.hpp"
#include "base/proteoform_factory.hpp"
#include "base/mass_constant.hpp"

#include "prsm/peak_ion_pair_util.hpp"

#include "spec/theo_peak_util.hpp"
#include "spec/extend_ms.hpp"

#include "local/local_mng.hpp"
#include "local/local_util.hpp"

namespace toppic {

namespace local_util {

void scrFilter(std::vector<double> & scr, int & bgn, int & end, double & conf, double threshold) {
  conf = 0.0;
  bgn = std::distance(scr.begin(), std::max_element(scr.begin(), scr.end()));
  if (scr[bgn] < threshold) {
    bgn = -1;
    return;
  }

  end = bgn;

  for (int i = bgn; i >= 0; i--) {
    if (scr[i] > threshold) {
      bgn = i;
    } else {
      scr[i] = 0;
    }
  }

  for (int i = end; i < static_cast<int>(scr.size()); i++) {
    if (scr[i] > threshold) {
      end = i;
    } else {
      scr[i] = 0;
    }
  }

  for (int i = bgn; i <= end; i++) { conf += scr[i];}

  std::vector<double> scr2;
  scr2.insert(scr2.end(), scr.begin() + bgn, scr.begin() + end + 1);
  scr = scr2;
}

PtmPtrVec getPtmPtrVecByMass(double mass, double err, const PtmPtrVec & ptm_vec) {
  PtmPtrVec res;
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    if (std::abs(ptm_vec[i]->getMonoMass() - mass) < err
        || std::abs(std::abs(ptm_vec[i]->getMonoMass() - mass) - mass_constant::getIsotopeMass()) < err) {
      res.push_back(ptm_vec[i]);
    }
  }
  return res;
}

PtmPairVec getPtmPairVecByMass(double mass1, double mass2, double err, const PtmPairVec & ptm_pair_vec) {
  PtmPairVec res;
  double mass = mass1 + mass2;
  for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
    double pair_mass = ptm_pair_vec[i].first->getMonoMass() + ptm_pair_vec[i].second->getMonoMass();

    if (std::abs(pair_mass - mass) < err
        || std::abs(std::abs(pair_mass - mass) - mass_constant::getIsotopeMass() ) < err)
      res.push_back(ptm_pair_vec[i]);
  }

  return res;
}

void compSupPeakNum(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                    MassShiftPtr mass_shift, double min_mass, int & left, int & right) {
  left = right = 0;
  PeakIonPairPtrVec pair_ptrs =
      peak_ion_pair_util::genePeakIonPairs(proteoform, extend_ms_ptr_vec, min_mass);

  for (size_t i = 0; i < pair_ptrs.size(); i++) {
    std::string ion_name =
        pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->getName();
    // A, B, C are for N-terminal ions
    if (ion_name == "A" || ion_name == "B" || ion_name == "C") {
      if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
          <= mass_shift->getLeftBpPos()) {
        left++;
      } else if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
                 >= mass_shift->getRightBpPos()) {
        right++;
      }
    } else {
      if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
          >= proteoform->getLen() - mass_shift->getLeftBpPos()) {
        left++;
      } else if (pair_ptrs[i]->getTheoPeakPtr()->getIonPtr()->getDisplayPos()
                 <= proteoform->getLen() - mass_shift->getRightBpPos()) {
        right++;
      }
    }
  }
}

void ptmMassAdjust(double & mass1, double & mass2, PtmPtr p1, PtmPtr p2) {
  if (p1 == nullptr || p2 == nullptr) return;
  double err = mass1 + mass2 - p1->getMonoMass() - p2->getMonoMass();
  if (std::abs(err) < mass_constant::getIsotopeMass()) {
    mass1 = p1->getMonoMass() + err / 2;
    mass2 = p2->getMonoMass() + err / 2;
  } else if (err > mass_constant::getIsotopeMass()) {
    err = err - mass_constant::getIsotopeMass();
    mass1 = p1->getMonoMass() + mass_constant::getIsotopeMass() + err / 2;
    mass2 = p2->getMonoMass() + err / 2;
  } else {
    err = err + mass_constant::getIsotopeMass();
    mass1 = p1->getMonoMass() - mass_constant::getIsotopeMass() + err / 2;
    mass2 = p2->getMonoMass() + err / 2;
  }
}

void fillTableB(std::vector<std::vector<double> > & b_table, double mass1, double mass2) { 
  for (size_t i = 1; i < b_table.size(); i++) {
    b_table[i].resize(b_table[0].size());
    std::fill(b_table[i].begin(), b_table[i].end(), 0);
  }

  for (size_t i = 0; i < b_table[0].size(); i++) {
    b_table[1][i] = b_table[0][i] + mass1;
    b_table[2][i] = b_table[1][i] + mass2;
  }
}

void compNumMatch(const std::vector<double> & b, std::vector<int> & s,
                  const ExtendMsPtr & extend_ms_ptr, double prec_mass) {
  std::vector<std::pair<double, double> > spec_peak
      = extend_ms::getExtendMassToleranceList(extend_ms_ptr);
  double n_shift = extend_ms_ptr->getMsHeaderPtr()->getActivationPtr()->getNShift();

  size_t i = 0, j = 0;
  while (i < b.size() && j < spec_peak.size()) {
    if (std::abs(spec_peak[j].first - b[i] - n_shift) <= spec_peak[j].second) {
      s[i]++; i++; j++;
    } else if (b[i] + n_shift > spec_peak[j].first) {
      j++;
    } else {
      i++;
    }
  }

  i = b.size() - 1, j = 0;
  while (static_cast<int>(i) >= 0 && j < spec_peak.size()) {
    if (std::abs(spec_peak[j].first - prec_mass + b[i] + n_shift) <= spec_peak[j].second) {
      s[i]++; i--; j++;
    } else if (prec_mass - b[i] - n_shift > spec_peak[j].first) {
      j++;
    } else {
      i--;
    }
  }
}

void fillTableS(std::vector<std::vector<double> > & b_table,
                std::vector<std::vector<int> > & s_table,
                ExtendMsPtr extend_ms_ptr, double prec_mass) {
  for (size_t i = 0; i < s_table.size(); i++) {
    s_table[i].resize(b_table[i].size() + 1);
    std::fill(s_table[i].begin(), s_table[i].end(), 0);
  }

  for (int j = 0; j < 3; j++) {
    compNumMatch(b_table[j], s_table[j], extend_ms_ptr, prec_mass);
  }
}

std::vector<double> geneNTheoMass(ProteoformPtr proteoform, ExtendMsPtr extend_ms_ptr,
                                  double min_mass) {
  TheoPeakPtrVec theo_peaks =
      theo_peak_util::geneProteoformTheoPeak(proteoform,
                                             extend_ms_ptr->getMsHeaderPtr()->getActivationPtr(),
                                             min_mass);

  std::vector<double> res;
  for (size_t i = 0; i < theo_peaks.size(); i++) {
    if (theo_peaks[i]->getIonPtr()->getIonTypePtr()->isNTerm()) {
      res.push_back(theo_peaks[i]->getModMass());
    }
  }
  std::sort(res.begin(), res.end());
  return res;
}

MassShiftPtrVec massShiftFilter(const MassShiftPtrVec & mass_shift_vec,
                                MassShiftTypePtr type) {
  MassShiftPtrVec res;
  for (size_t k = 0; k < mass_shift_vec.size(); k++) {
    if (mass_shift_vec[k]->getTypePtr() != type) {
      res.push_back(mass_shift_vec[k]);
    }
  }
  return res;
}

MassShiftPtrVec copyMassShiftVec(const MassShiftPtrVec & mass_shift_vec) {
  MassShiftPtrVec new_mass_shift_vec;
  for (size_t k = 0; k < mass_shift_vec.size(); k++) {
    MassShiftPtr mass_shift
        = std::make_shared<MassShift>(mass_shift_vec[k]->getLeftBpPos(),
                                      mass_shift_vec[k]->getRightBpPos(),
                                      mass_shift_vec[k]->getTypePtr());

    mass_shift->setChangePtr(mass_shift_vec[k]->getChangePtr(0));

    new_mass_shift_vec.push_back(mass_shift);
  }
  return new_mass_shift_vec;
}

double compMassShift(const MassShiftPtrVec & mass_shift_vec) {
  double mass = 0.0;
  for (size_t k = 0; k < mass_shift_vec.size(); k++) {
    mass += mass_shift_vec[k]->getMassShift();
  }
  return mass;
}

MassShiftPtr geneMassShift(MassShiftPtr shift, double mass, MassShiftTypePtr type) {
  ChangePtr change = std::make_shared<Change>(shift->getLeftBpPos(),
                                              shift->getRightBpPos(),
                                              type, mass,
                                              std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(),
                                                                    ResidueBase::getEmptyResiduePtr()));
  MassShiftPtr mass_shift = std::make_shared<MassShift>(shift->getLeftBpPos(),
                                                        shift->getRightBpPos(),
                                                        type);
  mass_shift->setChangePtr(change);
  return mass_shift;
}

void normalize(std::vector<double> & scr) {
  // to avoid overflow if the scores are too large
  double max = *std::max_element(scr.begin(), scr.end());
  for (size_t i = 0; i < scr.size(); i++) {
    scr[i] /= max;
  }
  double sum = std::accumulate(scr.begin(), scr.end(), 0.0);
  for (size_t i = 0; i < scr.size(); i++) {
    scr[i] /= sum;
  }
}

int compMatchFragNum(ProteoformPtr proteoform_ptr, const ExtendMsPtrVec &ms_ptr_vec, double min_mass) {
  return static_cast<int>(peak_ion_pair_util::genePeakIonPairs(proteoform_ptr, ms_ptr_vec, min_mass).size());
}

}  // namespace local_util

}  // namespace toppic
