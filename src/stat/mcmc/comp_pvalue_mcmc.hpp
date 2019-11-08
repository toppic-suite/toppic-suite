//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_STAT_MCMC_COMP_PVALUE_MCMC_HPP_
#define TOPPIC_STAT_MCMC_COMP_PVALUE_MCMC_HPP_

#include <random>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include "common/base/activation.hpp"
#include "prsm/prsm.hpp"
#include "stat/mcmc/mcmc_mng.hpp"

namespace toppic {

class CompPValueMCMC{
 public:
  CompPValueMCMC(MCMCMngPtr mng_ptr,
                 std::map<PtmPtr, std::vector<ResiduePtr> > ptm_residue_map,
                 std::map<int, std::vector<std::string> > mass_table):
      mng_ptr_(mng_ptr),
      generator_(new std::default_random_engine(42)),
      min_mass_(mng_ptr->prsm_para_ptr_->getSpParaPtr()->getMinMass()),
      ptm_residue_map_(ptm_residue_map),
      mass_table_(mass_table),
      ppo_(mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo()) {
        mu_.resize(mng_ptr_->n_);
        std::fill(mu_.begin(), mu_.end(), 1);
      }

  double compOneProbMCMC(PrsmPtr prsm_ptr, ActivationPtr act,
                         const std::vector<int> & ms_mass_int);

 private:
  void simulateDPR(ResiduePtrVec &residues, long omega, int scr_init, int k);

  ResiduePtrVec randomTrans(ResiduePtrVec residues);

  int getMaxScore(const ResiduePtrVec &residues);

  // if no ptm in the prsm
  int compScoreNoPtm();

  int getMaxScoreN(const ResiduePtrVec &residues);

  // this should be called first
  void initTheoMassWithPtm(const std::vector<size_t> & change_pos);

  void geneScrVec(std::vector<int> & n_scr_no_ptm,
                  std::vector<int> & n_scr_with_ptm,
                  std::vector<int> & c_scr_no_ptm,
                  std::vector<int> & c_scr_with_ptm,
                  double mass);

  void rmMassTheoMass(size_t pos, double mass);

  void addMassTheoMass(size_t pos, double mass);

  std::vector<int> compTheoMassPpos(const std::vector<double> &theo_masses);

  MCMCMngPtr mng_ptr_;

  std::default_random_engine * generator_;

  double min_mass_;

  std::map<PtmPtr, std::vector<ResiduePtr> > ptm_residue_map_;

  std::map<int, std::vector<std::string> > mass_table_;

  std::vector<int> score_vec_;

  int z_;

  double pep_mass_;

  double ppo_;

  std::vector<long long> mu_;

  ActivationPtr act_;

  PtmPtrVec ptm_vec_;

  std::vector<double> ptm_mass_vec_;

  std::vector<int> ms_mass_int_;

  std::vector<double> n_theo_masses_;

  std::vector<double> c_theo_masses_;

  std::vector<ResiduePtrVec> residues_stack_;

  std::vector<long> omega_stack_;

  std::vector<int> score_stack_;
};

typedef std::shared_ptr<CompPValueMCMC> CompPValueMCMCPtr;

}  // namespace toppic

#endif
