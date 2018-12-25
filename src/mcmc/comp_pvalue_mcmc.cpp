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

#include <vector>
#include <string>
#include <algorithm>

#include "common/base/residue_util.hpp"
#include "prsm/prsm_algo.hpp"

#include "mcmc/comp_pvalue_mcmc.hpp"

namespace toppic {

void getTheoMassVec(const ResiduePtrVec &residues,
                    IonTypePtr n_ion_type_ptr,
                    IonTypePtr c_ion_type_ptr,
                    double min_mass,
                    std::vector<double> & n_theo_masses,
                    std::vector<double> & c_theo_masses) {
  n_theo_masses.clear();

  c_theo_masses.clear();

  double res_mass = residue_util::compResiduePtrVecMass(residues);
  double max_mass = res_mass + mass_constant::getWaterMass() - min_mass;

  double prm = 0;

  for (size_t i = 0; i < residues.size(); i++) {
    prm += residues[i]->getMass();
    double srm = res_mass - prm;

    double n_mass = prm + n_ion_type_ptr->getShift();
    if (min_mass <= n_mass && n_mass <= max_mass) {
      n_theo_masses.push_back(n_mass);
    }

    double c_mass = srm + c_ion_type_ptr->getShift();
    if (min_mass <= c_mass && c_mass <= max_mass) {
      c_theo_masses.push_back(c_mass);
    }
  }

  std::sort(n_theo_masses.begin(), n_theo_masses.end());

  std::sort(c_theo_masses.begin(), c_theo_masses.end());
}

ResiduePtrVec CompPValueMCMC::randomTrans(ResiduePtrVec residues) {
  double mass = 0.0;
  bool is_small = (residue_util::compResiduePtrVecMass(residues) < this->pep_mass_);
  std::uniform_int_distribution<int> res_dist(0, residues.size() / 2 - 1);
  int pos1 = res_dist(*generator_);
  mass += residues[pos1]->getMass();

  int pos2 = pos1 + 1;

  mass += residues[pos2]->getMass();

  res_dist = std::uniform_int_distribution<int>(residues.size() / 2 + 1, residues.size() - 1);

  int pos3 = res_dist(*generator_);

  mass += residues[pos3]->getMass();

  ResiduePtrVec ori_res = {residues[pos1], residues[pos2], residues[pos3]};

  std::random_shuffle(ori_res.begin(), ori_res.end());

  std::vector<std::string> res_vec = mass_table_[std::round(mass * mng_ptr_->convert_ratio_)];

  ResiduePtrVec new_res_vec;

  if (res_vec.size() > 1) {
    res_dist = std::uniform_int_distribution<int>(0, res_vec.size() - 1);
    std::string res_seq = res_vec[res_dist(*generator_)];
    std::random_shuffle(res_seq.begin(), res_seq.end());
    new_res_vec
        = residue_util::convertStrToResiduePtrVec(res_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    new_res_vec = ori_res;
  }

  double new_mass = residue_util::compResiduePtrVecMass(new_res_vec);

  if (is_small) {
    if (new_mass < mass) {
      new_res_vec = ori_res;
    }
  } else {
    if (new_mass > mass) {
      new_res_vec = ori_res;
    }
  }

  residues[pos1] = new_res_vec[0];
  residues[pos2] = new_res_vec[1];
  residues[pos3] = new_res_vec[2];

  return residues;
}

double CompPValueMCMC::compOneProbMCMC(PrsmPtr prsm_ptr, ActivationPtr act,
                                       const std::vector<int> & ms_mass_int) {
  this->act_ = act;

  this->ms_mass_int_ = ms_mass_int;

  ProteoformPtr prot_form = prsm_ptr->getProteoformPtr();

  pep_mass_ = residue_util::compResiduePtrVecMass(prot_form->getResSeqPtr()->getResidues());

  ptm_vec_ = prot_form->getPtmVec(MassShiftType::VARIABLE);

  ptm_mass_vec_.resize(ptm_vec_.size());

  for (size_t i = 0; i < ptm_vec_.size(); i++) {
    ptm_mass_vec_[i] = ptm_vec_[i]->getMonoMass();
  }

  MassShiftPtrVec unknown_shift_vec = prot_form->getMassShiftPtrVec(MassShiftType::UNEXPECTED);

  if (unknown_shift_vec.size() > 0) {
    for (size_t k = 0; k < unknown_shift_vec.size(); k++) {
      ptm_vec_.push_back(nullptr);
      ptm_mass_vec_.push_back(unknown_shift_vec[k]->getMassShift());
    }
  }

  ResiduePtrVec residues = prot_form->getResSeqPtr()->getResidues();

  int scr = getMaxScore(residues);

  score_vec_.resize(mng_ptr_->N_ + 1);
  std::fill(score_vec_.begin(), score_vec_.end(), 0);

  std::vector<double> p(mng_ptr_->n_, 1.0);
  std::vector<int> n(mng_ptr_->n_, 0);
  double one_prob = 0.0;

  for (int k = 0; k < mng_ptr_->k_; k++) {
    this->generator_ = new std::default_random_engine(42);
    std::fill(n.begin(), n.end(), 0);
    std::fill(score_vec_.begin(), score_vec_.end(), 0);
    score_vec_[0] = scr;
    residues = prot_form->getResSeqPtr()->getResidues();

    if (k != 0) {
      for (size_t i = 0; i < mu_.size(); i++) {
        if (p[i] != 0.0) {
          mu_[i] = static_cast<long long>(1 / p[i]);
        } else {
          mu_[i] = mu_[i - 1];
        }
        // LOG_DEBUG("mu[" << i << "] " << mu_[i]);
      }
    }
    this->z_ = 0;
    // LOG_DEBUG("mu min " << *std::min_element(mu_.begin(), mu_.end()));

    long mu_min = std::round(*std::min_element(mu_.begin(), mu_.end()));

    simulateDPR(residues, mu_min, scr, k);

    for (size_t i = 0; i < score_vec_.size(); i++) {
      n[score_vec_[i]]++;
    }

    std::fill(p.begin(), p.end(), 0);

    for (size_t i = 0; i < n.size(); i++) {
      p[i] = n[i] * 1.0 / mu_[i];
    }

    double sum = std::accumulate(p.begin(), p.end(), 0.0);
    for (size_t i = 0; i < p.size(); i++) {
      p[i] = p[i] / sum;
      // LOG_DEBUG("p[" << i << "] " << p[i]);
    }

    if (prsm_ptr->getMatchFragNum() < 10) {
      int corrected_scr = std::min(scr, static_cast<int>(prsm_ptr->getMatchFragNum()));
      one_prob = std::accumulate(p.begin() + corrected_scr, p.end(), 0.0);
    } else {
      one_prob = std::accumulate(p.begin() + scr, p.end(), 0.0);
    }

    // LOG_DEBUG("k " << k << " one prob: " << one_prob);

    if (one_prob > 0.005 || one_prob < std::pow(10, -6)) {
      return one_prob;
    }
  }

  return std::max(one_prob, std::pow(10, -20));
}

int CompPValueMCMC::compScoreNoPtm() {
  std::vector<int> n_theo_masses_int(n_theo_masses_.size());

  for (size_t k = 0; k < n_theo_masses_.size(); k++) {
    n_theo_masses_int[k] = static_cast<int>(n_theo_masses_[k] * mng_ptr_->convert_ratio_) >> 5;
  }

  std::vector<int> n_match;

  std::set_intersection(ms_mass_int_.begin(), ms_mass_int_.end(),
                        n_theo_masses_int.begin(), n_theo_masses_int.end(),
                        std::back_inserter(n_match));

  std::vector<int> c_theo_masses_int(c_theo_masses_.size());

  for (size_t k = 0; k < c_theo_masses_.size(); k++) {
    c_theo_masses_int[k] = static_cast<int>(c_theo_masses_[k] * mng_ptr_->convert_ratio_) >> 5;
  }

  std::vector<int> c_match;

  std::set_intersection(ms_mass_int_.begin(), ms_mass_int_.end(),
                        c_theo_masses_int.begin(), c_theo_masses_int.end(),
                        std::back_inserter(c_match));

  return n_match.size() + c_match.size();
}

int CompPValueMCMC::getMaxScore(const ResiduePtrVec &residues) {
  getTheoMassVec(residues, act_->getNIonTypePtr(), act_->getCIonTypePtr(),
                 min_mass_, n_theo_masses_, c_theo_masses_);

  // no ptm. almost never happen
  if (this->ptm_vec_.size() == 0) return compScoreNoPtm();

  int max_scr = 0;

  for (int i = 0; i < 3; i++) {
    max_scr = std::max(getMaxScoreN(residues), max_scr);
  }

  return max_scr;
}

// update n_theo_masses, c_theo_masses
// this should be called first
void CompPValueMCMC::initTheoMassWithPtm(const std::vector<size_t> & change_pos) {
  std::vector<double> change_masses(n_theo_masses_.size(), 0.0);

  for (size_t i = 0; i < ptm_vec_.size(); i++) {
    // n theo mass won't change if it modifies the last residue
    if (change_pos[i] <= change_masses.size()) {
      double m = ptm_mass_vec_[i];
      std::for_each(change_masses.begin() + change_pos[i],
                    change_masses.end(), [m](double& d) { d += m;});
    }
  }

  for (size_t i = 0; i < n_theo_masses_.size(); i++) {
    n_theo_masses_[i] += change_masses[i];
  }

  change_masses.resize(c_theo_masses_.size());

  std::fill(change_masses.begin(), change_masses.end(), 0.0);

  for (size_t i = 0; i < ptm_vec_.size(); i++) {
    // c theo mass won't change if it modifies the first residue
    if (change_pos[i] > 0) {
      double m = ptm_mass_vec_[i];
      std::for_each(change_masses.begin() + change_masses.size() - change_pos[i],
                    change_masses.end(), [m](double& d) { d += m;});
    }
  }

  for (size_t i = 0; i < c_theo_masses_.size(); i++) {
    c_theo_masses_[i] += change_masses[i];
  }
}

std::vector<int> CompPValueMCMC::compTheoMassPpos(const std::vector<double> &theo_masses) {
  std::vector<int> theo_mass_int(theo_masses.size());
  for (size_t k = 0; k < theo_masses.size(); k++) {
    theo_mass_int[k] = static_cast<int>(theo_masses[k] * mng_ptr_->convert_ratio_) >> 5;
  }

  std::vector<int> results(theo_masses.size(), 0);

  // extendMsThree do not have 0 and precursor mass
  size_t i = 0;
  size_t j = 0;
  while (i < ms_mass_int_.size() && j < theo_mass_int.size()) {
    if (ms_mass_int_[i] == theo_mass_int[j]) {
      results[j] = 1;
      i++;
      j++;
    } else if (ms_mass_int_[i] > theo_mass_int[j]) {
      j++;
    } else {
      i++;
    }
  }

  return results;
}

void CompPValueMCMC::geneScrVec(std::vector<int> & n_scr_no_ptm,
                                std::vector<int> & n_scr_with_ptm,
                                std::vector<int> & c_scr_no_ptm,
                                std::vector<int> & c_scr_with_ptm,
                                double mass) {
  // n-term
  n_scr_no_ptm = compTheoMassPpos(n_theo_masses_);

  std::vector<double> n_theo_masses_ptm(n_theo_masses_);

  std::for_each(n_theo_masses_ptm.begin(), n_theo_masses_ptm.end(), [mass](double& d) { d += mass;});

  n_scr_with_ptm = compTheoMassPpos(n_theo_masses_ptm);

  // c-term
  c_scr_no_ptm = compTheoMassPpos(c_theo_masses_);

  std::vector<double> c_theo_masses_ptm(c_theo_masses_);

  std::for_each(c_theo_masses_ptm.begin(), c_theo_masses_ptm.end(), [mass](double& d) { d += mass;});

  c_scr_with_ptm = compTheoMassPpos(c_theo_masses_ptm);
}

int getMaxPosScrVec(const std::vector<size_t> & possible_change_pos,
                    const std::vector<int> & n_scr_no_ptm,
                    const std::vector<int> & n_scr_with_ptm,
                    const std::vector<int> & c_scr_no_ptm,
                    const std::vector<int> & c_scr_with_ptm, size_t & p) {
  int max_scr = 0;

  size_t prev_n_pos = 0;

  size_t prev_c_pos = c_scr_no_ptm.size() - prev_n_pos;

  int n_scr = std::accumulate(n_scr_with_ptm.begin(), n_scr_with_ptm.end(), 0);

  int c_scr = std::accumulate(c_scr_no_ptm.begin(), c_scr_no_ptm.end(), 0);

  for (size_t k = 0; k < possible_change_pos.size(); k++) {
    size_t n_pos = possible_change_pos[k];

    size_t c_pos = c_scr_with_ptm.size() - n_pos;

    for (size_t i = prev_n_pos; i < n_pos; i++) {
      n_scr -= n_scr_with_ptm[i];
      n_scr += n_scr_no_ptm[i];
    }

    for (size_t i = c_pos; i < prev_c_pos; i++) {
      c_scr += c_scr_with_ptm[i];
      c_scr -= c_scr_no_ptm[i];
    }

    int new_scr = n_scr + c_scr;

    if (new_scr > max_scr) {
      max_scr = new_scr;
      p = possible_change_pos[k];
    }

    prev_n_pos = n_pos;

    prev_c_pos = c_pos;
  }

  return max_scr;
}

void CompPValueMCMC::rmMassTheoMass(size_t pos, double mass) {
  // n theo mass won't change if it modifies the last residue
  if (pos <= n_theo_masses_.size()) {
    std::for_each(n_theo_masses_.begin() + pos,
                  n_theo_masses_.end(), [mass](double& d) { d -= mass;});
  }

  // c theo mass won't change if it modifies the first residue
  if (pos > 0) {
    std::for_each(c_theo_masses_.begin() + c_theo_masses_.size() - pos,
                  c_theo_masses_.end(), [mass](double& d) { d -= mass;});
  }
}

void CompPValueMCMC::addMassTheoMass(size_t pos, double mass) {
  // n theo mass won't change if it modifies the last residue
  if (pos <= n_theo_masses_.size()) {
    std::for_each(n_theo_masses_.begin() + pos,
                  n_theo_masses_.end(), [mass](double& d) { d += mass;});
  }

  // c theo mass won't change if it modifies the first residue
  if (pos > 0) {
    std::for_each(c_theo_masses_.begin() + c_theo_masses_.size() - pos,
                  c_theo_masses_.end(), [mass](double& d) { d += mass;});
  }
}

int CompPValueMCMC::getMaxScoreN(const ResiduePtrVec &residues) {
  std::vector<std::vector<size_t> > possible_change_pos(ptm_vec_.size());
  std::vector<size_t> change_pos(ptm_vec_.size());

  for (size_t i = 0; i < ptm_vec_.size(); i++) {
    if (ptm_vec_[i] != nullptr) {
      ResiduePtrVec possible_res = ptm_residue_map_[ptm_vec_[i]];
      for (size_t k = 0; k < residues.size(); k++) {
        if (std::find(possible_res.begin(), possible_res.end(), residues[k]) != possible_res.end()) {
          possible_change_pos[i].push_back(k);
        }
      }
    } else {
      for (size_t k = 0; k < residues.size(); k++) {
        possible_change_pos[i].push_back(k);
      }
    }
  }

  // random init positions
  for (size_t i = 0; i < possible_change_pos.size(); i++) {
    // no possible modified position on the peptides
    if (possible_change_pos[i].size() == 0) {
      return 0;
    } else {
      std::uniform_int_distribution<size_t> dis(0, possible_change_pos[i].size() - 1);
      change_pos[i] = possible_change_pos[i][dis(*generator_)];
    }
  }

  int max_scr = 0;

  std::vector<int> n_scr_no_ptm;

  std::vector<int> n_scr_with_ptm;

  std::vector<int> c_scr_no_ptm;

  std::vector<int> c_scr_with_ptm;

  initTheoMassWithPtm(change_pos);

  for (size_t p = 0; p < this->ptm_vec_.size(); p++) {
    rmMassTheoMass(change_pos[p], ptm_mass_vec_[p]);
    geneScrVec(n_scr_no_ptm, n_scr_with_ptm, c_scr_no_ptm, c_scr_with_ptm, ptm_mass_vec_[p]);
    size_t new_pos = change_pos[p];
    int new_max_scr = getMaxPosScrVec(possible_change_pos[p],
                                      n_scr_no_ptm,
                                      n_scr_with_ptm,
                                      c_scr_no_ptm,
                                      c_scr_with_ptm,
                                      new_pos);
    if (new_max_scr > max_scr) {
      change_pos[p] = new_pos;
      max_scr = new_max_scr;
    }

    addMassTheoMass(change_pos[p], ptm_mass_vec_[p]);
  }

  return max_scr;
}

void CompPValueMCMC::simulateDPR(ResiduePtrVec &residues, long omega, int scr_init, int k) {
  size_t CT_LIMIT = 20000;

  residues_stack_.reserve(CT_LIMIT);
  omega_stack_.reserve(CT_LIMIT);
  score_stack_.reserve(CT_LIMIT);

  residues_stack_.push_back(residues);
  omega_stack_.push_back(omega);
  score_stack_.push_back(scr_init);

  while (this->z_ < mng_ptr_->N_) {
    while (!residues_stack_.empty()) {
      ResiduePtrVec residues1 = residues_stack_.back();
      residues_stack_.pop_back();
      long omega1 = omega_stack_.back();
      omega_stack_.pop_back();
      int score1 = score_stack_.back();
      score_stack_.pop_back();

      ResiduePtrVec residues2 = randomTrans(residues1);

      int score2 = getMaxScore(residues2);

      if (mu_[score2] < omega1) {
        continue;
      }

      if (mu_[score2] > mu_[score1] && residues_stack_.size() < CT_LIMIT) {
        long Y = std::round(mu_[score2] / mu_[score1]);
        Y = std::min(Y, 100L);
        std::uniform_int_distribution<long> omega_dist(static_cast<long>(mu_[score1]), static_cast<long>(mu_[score2]));
        for (long i = 1; i < Y; i++) {
          long omega2 = omega_dist(*generator_);
          residues_stack_.push_back(residues2);
          omega_stack_.push_back(omega2);
          score_stack_.push_back(getMaxScore(residues2));
        }
      } else {
        this->z_++;
        // if (this->z_ % 500 == 0) {
        // LOG_DEBUG(this->z_);
        // }

        if (this->z_ > mng_ptr_->N_) {
          residues_stack_.clear();
          omega_stack_.clear();
          score_stack_.clear();
          return;
        }

        this->score_vec_[this->z_] = score2;

        residues_stack_.push_back(residues2);

        omega_stack_.push_back(omega1);

        score_stack_.push_back(score2);
      }
    }

    ResiduePtrVec residues3 = randomTrans(residues);
    residues_stack_.push_back(residues3);
    omega_stack_.push_back(omega);
    score_stack_.push_back(getMaxScore(residues3));
  }
}

}  // namespace toppic
