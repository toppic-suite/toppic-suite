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

#include "base/residue_util.hpp"

#include "mcmc/comp_pvalue_mcmc.hpp"

namespace prot {

std::vector<double> getNTheoMassVec(const ResiduePtrVec &residues,
                                    IonTypePtr n_ion_type_ptr,
                                    double min_mass) {
  std::vector<double> theo_masses;
  ResSeqPtr res_seq = std::make_shared<ResidueSeq>(residues);
  BpSpecPtr bp_spec = std::make_shared<BpSpec>(res_seq);
  BreakPointPtrVec bps = bp_spec->getBreakPointPtrVec();
  double max_mass = res_seq->getSeqMass() - min_mass;
  for (size_t i = 0; i < bps.size(); i++) {
    double n_mass = bps[i]->getNTermMass(n_ion_type_ptr);
    if (min_mass <= n_mass && n_mass <= max_mass) {
      theo_masses.push_back(n_mass);
    }
  }
  return theo_masses;
}

std::vector<double> getCTheoMassVec(const ResiduePtrVec &residues,
                                    IonTypePtr c_ion_type_ptr,
                                    double min_mass) {
  std::vector<double> theo_masses;
  ResSeqPtr res_seq = std::make_shared<ResidueSeq>(residues);
  BpSpecPtr bp_spec = std::make_shared<BpSpec>(res_seq);
  BreakPointPtrVec bps = bp_spec->getBreakPointPtrVec();
  double max_mass = res_seq->getSeqMass() - min_mass;
  for (size_t i = 0; i < bps.size(); i++) {
    double n_mass = bps[i]->getCTermMass(c_ion_type_ptr);
    if (min_mass <= n_mass && n_mass <= max_mass) {
      theo_masses.push_back(n_mass);
    }
  }
  return theo_masses;
}

ResiduePtrVec CompPValueMCMC::randomTrans(ResiduePtrVec residues) {
  double mass = 0.0;
  bool is_small = (residue_util::compResiduePtrVecMass(residues) < this->pep_mass_);
  std::uniform_int_distribution<int> res_dist(0, residues.size() / 2 - 1);
  int pos1 = res_dist(mt_);
  mass += residues[pos1]->getMass();

  int pos2 = pos1 + 1;

  mass += residues[pos2]->getMass();

  res_dist = std::uniform_int_distribution<int>(residues.size() / 2 + 1, residues.size() - 1);

  int pos3 = res_dist(mt_);

  mass += residues[pos3]->getMass();

  ResiduePtrVec ori_res = {residues[pos1], residues[pos2], residues[pos3]};

  std::random_shuffle(ori_res.begin(), ori_res.end());

  std::vector<std::string> res_vec = mass_table_[std::round(mass * mng_ptr_->convert_ratio_)];

  ResiduePtrVec new_res_vec;

  if (res_vec.size() > 1) {
    res_dist = std::uniform_int_distribution<int>(0, res_vec.size() - 1);
    std::string res_seq = res_vec[res_dist(mt_)];
    std::random_shuffle(res_seq.begin(), res_seq.end());
    new_res_vec = residue_util::convertStrToResiduePtrVec(res_seq);
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

int CompPValueMCMC::compNumMatched(const std::vector<int> &ms_masses,
                                   const std::vector<double> &theo_masses) {
  std::vector<int> theo_mass_int(theo_masses.size());

  for (size_t k = 0; k < theo_masses.size(); k++) {
    theo_mass_int[k] = static_cast<int>(theo_masses[k] * mng_ptr_->convert_ratio_) >> 5;
  }

  std::sort(theo_mass_int.begin(), theo_mass_int.end());
  theo_mass_int.erase(std::unique(theo_mass_int.begin(), theo_mass_int.end()), theo_mass_int.end());

  std::vector<int> match;

  std::set_intersection(ms_masses.begin(), ms_masses.end(),
                        theo_mass_int.begin(), theo_mass_int.end(),
                        std::back_inserter(match));

  return static_cast<int>(match.size());
}

double CompPValueMCMC::compPValueMCMC(PrsmPtr prsm_ptr, ActivationPtr act,
                                      const std::vector<int> & ms_masses) {
  ProteoformPtr prot_form = prsm_ptr->getProteoformPtr();

  pep_mass_ = residue_util::compResiduePtrVecMass(prot_form->getResSeqPtr()->getResidues());
  PtmPtrVec ptm_vec = prot_form->getPtmVec();

  ResiduePtrVec residues = prot_form->getResSeqPtr()->getResidues();

  int scr = getMaxScore(residues, ms_masses, act, ptm_vec);

  LOG_DEBUG("matching score: " << scr);

  score_vec_.resize(mng_ptr_->N_ + 1);
  std::fill(score_vec_.begin(), score_vec_.end(), 0);

  std::vector<long long> mu(mng_ptr_->n_, 1);
  std::vector<double> p(mng_ptr_->n_, 1.0);
  std::vector<int> n(mng_ptr_->n_, 0);
  double p_value = 0.0;

  for (int k = 0; k < mng_ptr_->k_; k++) {
    std::fill(n.begin(), n.end(), 0);
    std::fill(score_vec_.begin(), score_vec_.end(), 0);
    score_vec_[0] = scr;
    residues = prot_form->getResSeqPtr()->getResidues();

    if (k != 0) {
      for (size_t i = 0; i < mu.size(); i++) {
        if (p[i] != 0.0) {
          mu[i] = static_cast<long long>(1 / p[i]);
        } else {
          mu[i] = mu[i - 1];
        }
      }
    }
    this->z_ = 0;
    LOG_DEBUG("mu min " << *std::min_element(mu.begin(), mu.end()));

    long mu_min = std::round(*std::min_element(mu.begin(), mu.end()));

    simulateDPR(residues, ms_masses, act, ptm_vec, mu_min, mu);

    for (size_t i = 0; i < score_vec_.size(); i++) {
      n[score_vec_[i]]++;
    }

    for (size_t i = 0; i < n.size(); i++) {
      p[i] = n[i] * 1.0 / mu[i];
    }

    double sum = std::accumulate(p.begin(), p.end(), 0.0);
    for (size_t i = 0; i < p.size(); i++) {
      p[i] = p[i] / sum;
    }

    p_value = std::accumulate(p.begin() + scr, p.end(), 0.0);

    LOG_DEBUG("k " << k << " p-value: " << p_value);

    if (p_value > 0.5) {
      return 1.0;
    }
  }

  return p_value;
}

int CompPValueMCMC::compScore(const std::vector<int> & ms_masses,
                              std::vector<double> n_theo_masses,
                              std::vector<double> c_theo_masses,
                              const std::vector<size_t> & change_pos,
                              const PtmPtrVec & ptm_vec) {
  std::vector<double> change_masses(n_theo_masses.size(), 0.0);
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    double m = ptm_vec[i]->getMonoMass();
    std::for_each(change_masses.begin() + change_pos[i],
                  change_masses.end(), [m](double& d) { d += m;});
  }

  for (size_t i = 0; i < n_theo_masses.size(); i++) {
    n_theo_masses[i] += change_masses[i];
  }

  int scr = compNumMatched(ms_masses, n_theo_masses);

  change_masses.resize(c_theo_masses.size());
  std::fill(change_masses.begin(), change_masses.end(), 0.0);
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    double m = ptm_vec[i]->getMonoMass();
    std::for_each(change_masses.begin() + change_masses.size() - change_pos[i],
                  change_masses.end(), [m](double& d) { d += m;});
  }

  for (size_t i = 0; i < c_theo_masses.size(); i++) {
    c_theo_masses[i] += change_masses[i];
  }

  scr += compNumMatched(ms_masses, c_theo_masses);

  return scr;
}

int CompPValueMCMC::getMaxScore(const ResiduePtrVec &residues,
                                const std::vector<int> & ms_masses,
                                ActivationPtr act, const PtmPtrVec & ptm_vec) {
  std::vector<double> n_theo_masses = getNTheoMassVec(residues, act->getNIonTypePtr(), min_mass_);
  std::sort(n_theo_masses.begin(), n_theo_masses.end());

  std::vector<double> c_theo_masses = getCTheoMassVec(residues, act->getCIonTypePtr(), min_mass_);
  std::sort(c_theo_masses.begin(), c_theo_masses.end());

  std::vector<std::vector<size_t> > possible_change_pos(ptm_vec.size());
  std::vector<size_t> change_pos(ptm_vec.size());

  if (ptm_vec.size() == 0) {
    return compScore(ms_masses, n_theo_masses, c_theo_masses, change_pos, ptm_vec);
  }

  for (size_t i = 0; i < ptm_vec.size(); i++) {
    ResiduePtrVec possible_res = ptm_residue_map_[ptm_vec[i]];
    for (size_t k = 0; k < residues.size(); k++) {
      if (std::find(possible_res.begin(), possible_res.end(), residues[k]) != possible_res.end()) {
        possible_change_pos[i].push_back(k);
      }
    }
  }

  // random init positions
  for (size_t i = 0; i < possible_change_pos.size(); i++) {
    if (possible_change_pos[i].size() == 0) {
      return 0;
    } else {
      std::uniform_int_distribution<size_t> dis(0, possible_change_pos[i].size() - 1);
      change_pos[i] = dis(mt_);
    }
  }

  int max_scr = compScore(ms_masses, n_theo_masses, c_theo_masses, change_pos, ptm_vec);

  // search the max score using greedy method
  for (size_t i = 0; i < possible_change_pos.size(); i++) {
    size_t max_idx = change_pos[i];
    for (size_t k = 0; k < possible_change_pos[i].size(); k++) {
      change_pos[i] = possible_change_pos[i][k];
      int scr = compScore(ms_masses, n_theo_masses, c_theo_masses, change_pos, ptm_vec);
      if (scr > max_scr) {
        max_idx = possible_change_pos[i][k];
        max_scr = scr;
      }
    }
    change_pos[i] = max_idx;
  }

  return max_scr;
}

void CompPValueMCMC::simulateDPR(ResiduePtrVec &residues, const std::vector<int> & ms_masses,
                                 ActivationPtr act, const PtmPtrVec & ptm_vec,
                                 long omega, const std::vector<long long> & mu) {
  int score = getMaxScore(residues, ms_masses, act, ptm_vec);
  while (this->z_ < mng_ptr_->N_) {
    ResiduePtrVec residues2 = randomTrans(residues);

    int score2 = getMaxScore(residues2, ms_masses, act, ptm_vec);

    if (mu[score2] < omega) return;

    if (mu[score2] > mu[score]) {
      long Y = std::round(mu[score2] / mu[score]);
      Y = std::min(Y, 100L);
      std::uniform_int_distribution<long> omega_dist(static_cast<long>(mu[score]), static_cast<long>(mu[score2]));
      for (long i = 1; i < Y; i++) {
        long omega2 = omega_dist(mt_);
        ResiduePtrVec residues3 = residues2;
        simulateDPR(residues3, ms_masses, act, ptm_vec, omega2, mu);
      }
    }

    if (this->z_ > mng_ptr_->N_) {
      return;
    }

    this->z_++;

    this->score_vec_[this->z_] = score2;

    residues = residues2;

    score = score2;
  }
}

}  // namespace prot
