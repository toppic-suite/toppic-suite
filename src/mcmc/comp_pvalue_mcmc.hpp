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

#ifndef PROT_COMP_PVALUE_MCMC_HPP_
#define PROT_COMP_PVALUE_MCMC_HPP_

#include <random>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include "base/activation.hpp"

#include "prsm/prsm.hpp"

#include "mcmc_mng.hpp"

namespace prot {

class CompPValueMCMC{
 public:
  CompPValueMCMC(MCMCMngPtr mng_ptr,
                 std::map<PtmPtr, std::vector<ResiduePtr> > ptm_residue_map,
                 std::map<int, std::vector<std::string> > mass_table):
      mng_ptr_(mng_ptr),
      mt_(42),
      min_mass_(mng_ptr->prsm_para_ptr_->getSpParaPtr()->getMinMass()),
      ptm_residue_map_(ptm_residue_map),
      mass_table_(mass_table) {
        mu_.resize(mng_ptr_->n_);
        std::fill(mu_.begin(), mu_.end(), 1);
      }

  double compPValueMCMC(PrsmPtr prsm_ptr, ActivationPtr act,
                        const std::vector<int> & ms_masses);

 private:
  void simulateDPR(ResiduePtrVec &residues, long omega);

  ResiduePtrVec randomTrans(ResiduePtrVec residues);

  int getMaxScore(const ResiduePtrVec &residues);

  int compScore(std::vector<double> n_theo_masses, std::vector<double> c_theo_masses,
                const std::vector<size_t> & change_pos);

  int compNumMatched(const std::vector<int> &ms_masses, const std::vector<double> &theo_masses);

  MCMCMngPtr mng_ptr_;

  std::mt19937 mt_;

  double min_mass_;

  std::map<PtmPtr, std::vector<ResiduePtr> > ptm_residue_map_;

  std::map<int, std::vector<std::string> > mass_table_;

  std::vector<int> score_vec_;

  int z_;

  double pep_mass_;

  std::vector<long long> mu_;

  ActivationPtr act_;

  PtmPtrVec ptm_vec_;

  std::vector<int> ms_masses_;
};

typedef std::shared_ptr<CompPValueMCMC> CompPValueMCMCPtr;

}  // namespace prot

#endif
