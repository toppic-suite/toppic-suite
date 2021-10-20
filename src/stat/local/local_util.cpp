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
#include <numeric>
#include <iomanip>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/factory/extend_ms_util.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_algo.hpp"
#include "prsm/theo_peak_util.hpp"
#include "stat/local/local_util.hpp"

namespace toppic {

namespace local_util {

MassShiftPtrVec massShiftFilter(const MassShiftPtrVec & mass_shift_vec,
                                AlterTypePtr type) {
  MassShiftPtrVec res;
  for (size_t k = 0; k < mass_shift_vec.size(); k++) {
    if (mass_shift_vec[k]->getTypePtr() != type) {
      res.push_back(mass_shift_vec[k]);
    }
  }
  return res;
}

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

PtmPairVec getPtmPairVecByMass(double mass, double err, const PtmPairVec & ptm_pair_vec) {
  PtmPairVec res;
  for (size_t i = 0; i < ptm_pair_vec.size(); i++) {
    double pair_mass = ptm_pair_vec[i].first->getMonoMass() + ptm_pair_vec[i].second->getMonoMass();

    if (std::abs(pair_mass - mass) < err
        || std::abs(std::abs(pair_mass - mass) - mass_constant::getIsotopeMass() ) < err)
      res.push_back(ptm_pair_vec[i]);
  }

  return res;
}

// compute match table
void compMTable(const std::vector<double> &ori_theo_masses, double shift, 
    const std::vector<double> &exp_masses, 
    PeakTolerancePtr tole_ptr,  
    std::vector<int> &m_table) {

  std::vector<double> theo_masses;
  for (size_t k = 0; k < ori_theo_masses.size(); k++) {
    theo_masses.push_back(ori_theo_masses[k] + shift);
  }
  size_t i = 0;
  size_t j = 0;
  while (i < exp_masses.size() && j < theo_masses.size()) {
    double deviation = exp_masses[i] - theo_masses[j];
    double err = tole_ptr->compStrictErrorTole(exp_masses[i]); 
    if (std::abs(deviation) <= err) {
      if (j > 0 && j < theo_masses.size() - 1) {
        if (m_table[j] == 0) {
          m_table[j] = 1;
        }
      }
    }
    if (prsm_algo::increaseIJ(i, j, deviation, err, exp_masses, theo_masses)) {
      i++;
    } else {
      j++;
    }
  }
}

void addTwoVectors(std::vector<int> &row_1, std::vector<int> &row_2) {
  for (size_t i = 0; i < row_1.size(); i++) {
    row_1[i] = row_1[i] + row_2[i];
  }
}

void addReverseTwoVectors(std::vector<int> &row_1, std::vector<int> &row_2) {
  int len = row_1.size();
  for (size_t i = 0; i < row_1.size(); i++) {
    row_1[i] = row_1[i] + row_2[len - 1 - i];
  }
}


std::vector<int> compPrefScore(std::vector<int> & row) {
  std::vector<int> result(row.size(), 0);
  int total_score = 0;
  for (size_t i = 0; i < row.size(); i++) {
    total_score = total_score + row[i];
    result[i] = total_score;
  }
  return result;
}

std::vector<int> compSuffScore(std::vector<int> & row) {
  std::vector<int> result(row.size(), 0);
  int total_score = 0;
  for (int i = row.size() - 1; i >=0; i--) {
    total_score = total_score + row[i];
    result[i] = total_score; 
  }
  return result;
}

void compOnePtmSTable(std::vector<int> &s_table, ProteoformPtr form_ptr, 
                      const ExtendMsPtrVec & extend_ms_ptr_vec,
                      PtmPtr ptm_ptr, LocalMngPtr mng_ptr) {

  int len = form_ptr->getLen();
  std::vector<int> row_1(len + 1, 0);
  std::vector<int> row_2(len + 1, 0);

  //LOG_DEBUG("s table 1");
  BpSpecPtr bp_spec_ptr = form_ptr->getBpSpecPtr();
  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  // srm is in the increasing order 
  std::vector<double> srm_masses = bp_spec_ptr->getSrmMasses();
  //LOG_DEBUG("prm size " << prm_masses.size() << " len " << len);
  //LOG_DEBUG("srm size " << srm_masses.size() << " len " << len);
  PeakTolerancePtr tole_ptr = mng_ptr->peak_tole_ptr_;

  double shift_mass = ptm_ptr->getMonoMass();
  for (size_t i = 0; i < extend_ms_ptr_vec.size(); i++) {
    std::vector<double> ms_masses = extend_ms_util::getExtendMassVec(extend_ms_ptr_vec[i]);
    // updated Match table using prm masses 
    double n_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getN_BYShift();
    //LOG_DEBUG("n shift " << n_shift << " shift " << shift_mass);
    std::vector<int> row_1_n(len + 1, 0);
    local_util::compMTable(prm_masses, n_shift, ms_masses, tole_ptr, row_1_n);
    local_util::addTwoVectors(row_1, row_1_n);

    std::vector<int> row_2_n(len + 1, 0);
    local_util::compMTable(prm_masses, n_shift + shift_mass, ms_masses, tole_ptr, row_2_n);  
    local_util::addTwoVectors(row_2, row_2_n);

    // updated Match table using srm masses 
    double c_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getC_BYShift() 
      + mass_constant::getWaterMass();
    //ActivationPtr act_ptr = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr();
    //LOG_DEBUG("activation type " << act_ptr->getName());
    //LOG_DEBUG(std::setprecision(8) << "c shift " << c_shift);
    std::vector<int> row_1_c(len + 1, 0);
    LOG_DEBUG(std::setprecision(8) << "c shift " << c_shift);
    local_util::compMTable(srm_masses, c_shift + shift_mass, ms_masses, tole_ptr, row_1_c);  
    local_util::addReverseTwoVectors(row_1, row_1_c);

    std::vector<int> row_2_c(len + 1, 0);
    local_util::compMTable(srm_masses, c_shift, ms_masses, tole_ptr, row_2_c);  
    local_util::addReverseTwoVectors(row_2, row_2_c);
  }

  std::vector<int> n_prec_score = local_util::compPrefScore(row_1);
  //LOG_DEBUG("prec_score size " << n_prec_score.size() << " len " << len);
  std::vector<int> c_suff_score = local_util::compSuffScore(row_2);
  //LOG_DEBUG("suff score size " << c_suff_score.size() << " len " << len);

  // get result table
  for (int i = 0; i < len; i++) {
    s_table.push_back(n_prec_score[i] + c_suff_score[i+1]);
  }
}

void compTwoPtmMTable(std::vector<std::vector<int>> &m_table, ProteoformPtr form_ptr, 
                      const ExtendMsPtrVec & extend_ms_ptr_vec,
                      PtmPtr ptm_ptr_1, PtmPtr ptm_ptr_2, LocalMngPtr mng_ptr) {

  int len = form_ptr->getLen();
  for (int i = 0; i < 3; i++) {
    std::vector<int> row(len + 1, 0);
    m_table.push_back(row);
  }

  BpSpecPtr bp_spec_ptr = form_ptr->getBpSpecPtr();
  std::vector<double> prm_masses = bp_spec_ptr->getPrmMasses();
  //srm mass is in the increasing order
  std::vector<double> srm_masses = bp_spec_ptr->getSrmMasses();
  PeakTolerancePtr tole_ptr = mng_ptr->peak_tole_ptr_;

  double mass_1 = ptm_ptr_1->getMonoMass();
  double mass_2 = ptm_ptr_2->getMonoMass();

  for (size_t i = 0; i < extend_ms_ptr_vec.size(); i++) {
    std::vector<double> ms_masses = extend_ms_util::getExtendMassVec(extend_ms_ptr_vec[i]);
    // updated S table using prm masses 
    double n_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getN_BYShift();
    std::vector<int> row_1_n(len + 1, 0);
    compMTable(prm_masses, n_shift, ms_masses, tole_ptr, row_1_n);
    addTwoVectors(m_table[0], row_1_n);

    std::vector<int> row_2_n(len + 1, 0);
    compMTable(prm_masses, n_shift + mass_1, ms_masses, tole_ptr, row_2_n);
    addTwoVectors(m_table[1], row_2_n);

    std::vector<int> row_3_n(len + 1, 0);
    compMTable(prm_masses, n_shift + mass_1 + mass_2, ms_masses, tole_ptr, row_3_n);  
    addTwoVectors(m_table[2], row_3_n);

    // updated S table using srm masses 
    double c_shift = extend_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getC_BYShift()
          + mass_constant::getWaterMass();
    std::vector<int> row_1_c(len + 1, 0);
    compMTable(srm_masses, c_shift + mass_1 + mass_2, ms_masses, tole_ptr, row_1_c);  
    addReverseTwoVectors(m_table[0], row_1_c);

    std::vector<int> row_2_c(len + 1, 0);
    compMTable(srm_masses, c_shift + mass_2, ms_masses, tole_ptr, row_2_c);  
    addReverseTwoVectors(m_table[1], row_2_c);

    std::vector<int> row_3_c(len + 1, 0);
    compMTable(srm_masses, c_shift, ms_masses, tole_ptr, row_3_c);  
    addReverseTwoVectors(m_table[2], row_3_c);
  }

  /*
  for (int j = 0; j < len+1; j++) {
    LOG_DEBUG(std::setprecision(8) << " m_table[0]  " << j << " " << prm_masses[j] << " "  << m_table[0][j]);
  }

  for (int j = 0; j < len+1; j++) {
    LOG_DEBUG(std::setprecision(8) << " m_table[1]  " << j << " " << prm_masses[j] << " "  << m_table[1][j]);
  }

  for (int j = 0; j < len+1; j++) {
    LOG_DEBUG(std::setprecision(8) << " m_table[2]  " << j << " " << prm_masses[j] << " "  << m_table[2][j]);
  }
  */
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

}

}
