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

#include "prsm/prsm_algo.hpp"
#include "search/diag/diag_header.hpp"

namespace toppic {

DiagHeader::DiagHeader(double n_term_shift,
                       bool n_strict, bool c_strict,
                       bool prot_n_match, bool prot_c_match,
                       bool pep_n_match, bool pep_c_match):

  prot_N_term_shift_(n_term_shift),
  n_strict_(n_strict),
  c_strict_(c_strict),
  prot_N_term_match_(prot_n_match),
  prot_C_term_match_(prot_c_match),
  pep_N_term_match_(pep_n_match),
  pep_C_term_match_(pep_c_match) {}

void DiagHeader::changeOnlyNTermShift(double s) {
  prot_N_term_shift_ += s;
  pep_N_term_shift_  += s;
}

void DiagHeader::changeOnlyCTermShift(double s) {
  prot_C_term_shift_ += s;
  pep_C_term_shift_  += s;
}

void DiagHeader::initHeader(double c_shift, ProteoformPtr proteo_ptr) {
  // set protein c term shift
  prot_C_term_shift_ = c_shift;
  std::vector<double> prm_masses = proteo_ptr->getBpSpecPtr()->getPrmMasses();

  trunc_first_res_pos_ = prsm_algo::getFirstResPos(prot_N_term_shift_,
                                                   prm_masses);

  trunc_last_res_pos_ = prsm_algo::getLastResPos(prot_C_term_shift_, prm_masses);
  pep_N_term_shift_ = prot_N_term_shift_ + prm_masses[trunc_first_res_pos_];
  pep_C_term_shift_ = prot_C_term_shift_ + prm_masses[prm_masses.size() - 1]
    - prm_masses[trunc_last_res_pos_ + 1];
}

}  // namespace toppic
