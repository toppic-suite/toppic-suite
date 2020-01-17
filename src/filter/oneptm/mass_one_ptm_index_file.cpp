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

#include "seq/proteoform_util.hpp"
#include "ms/spec/prm_ms.hpp"
#include "filter/massmatch/filter_protein.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/oneptm/mass_one_ptm_filter.hpp"
#include "filter/oneptm/mass_one_ptm_index_file.hpp"

namespace toppic {

MassOnePtmIndex::MassOnePtmIndex(const ProteoformPtrVec &proteo_ptrs,
                                   OnePtmFilterMngPtr mng_ptr, std::vector<std::string> file_vec) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  std::vector<std::vector<double>> shift_2d
      = proteoform_util::getNTermShift2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
  std::vector<std::vector<double>> n_term_acet_2d
      = proteoform_util::getNTermAcet2D(proteo_ptrs, mng_ptr->prsm_para_ptr_->getProtModPtrVec());
  term_index_ptr_ = MassMatchFactory::getPrmTermMassMatchPtr(proteo_ptrs, shift_2d,
                                                             mng_ptr->max_proteoform_mass_,
                                                             mng_ptr->filter_scale_);
  diag_index_ptr_ = MassMatchFactory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                             mng_ptr->max_proteoform_mass_,
                                                             mng_ptr->filter_scale_);
  std::vector<std::vector<double>> rev_shift_2d;
  std::vector<double> shift_1d(1, 0);
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    rev_shift_2d.push_back(shift_1d);
  }
  rev_term_index_ptr_ = MassMatchFactory::getSrmTermMassMatchPtr(proteo_ptrs, rev_shift_2d, n_term_acet_2d,
                                                                 mng_ptr->max_proteoform_mass_,
                                                                 mng_ptr->filter_scale_);
  rev_diag_index_ptr_ = MassMatchFactory::getSrmDiagMassMatchPtr(proteo_ptrs, n_term_acet_2d,
                                                                 mng_ptr->max_proteoform_mass_,
                                                                 mng_ptr->filter_scale_);
  term_index_ptr_->setfileName(file_vec[0]);
  diag_index_ptr_->setfileName(file_vec[1]);
  rev_term_index_ptr_->setfileName(file_vec[2]);
  rev_diag_index_ptr_->setfileName(file_vec[3]);

  term_index_ptr_->serializeMassMatch();
  diag_index_ptr_->serializeMassMatch();
  rev_term_index_ptr_->serializeMassMatch();
  rev_diag_index_ptr_->serializeMassMatch();  
}
} /* namespace toppic */
