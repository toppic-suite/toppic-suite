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
#include <iostream>

#include "ms/spec/prm_ms.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "filter/massmatch/filter_protein.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/zeroindex/topindex_file_name.hpp"
#include "filter/diag/mass_diag_filter.hpp"
#include "filter/diagindex/diag_index_file.hpp"

namespace toppic {

DiagIndexFile::DiagIndexFile(const ProteoformPtrVec &proteo_ptrs,
                               DiagFilterMngPtr mng_ptr, std::string block_str) {
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;

  index_ptr_ = MassMatchFactory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                        mng_ptr->max_proteoform_mass_,
                                                        mng_ptr->filter_scale_);
                                                        
  TopIndexFileNamePtr file_name_ptr = std::make_shared<TopIndexFileName>();
  std::string parameters = file_name_ptr->geneFileName(prsm_para_ptr);

  std::string dir_name = mng_ptr_->prsm_para_ptr_->getOriDbName() + "_idx";
  std::string file_name = file_name_ptr->multi_ptm_file_vec_[0] + parameters + block_str;

  index_ptr_->serializeMassMatch(dir_name, file_name);

}

} /* namespace toppic */
