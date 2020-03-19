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

#include <iostream>

#include "filter/massmatch/mass_match_factory.hpp"
#include "filter/diagindex/diag_index_file.hpp"

namespace toppic {

namespace diag_index_file {

void geneDiagIndexFile(const ProteoformPtrVec &proteo_ptrs,
                       DiagFilterMngPtr mng_ptr, std::string block_str) {

  MassMatchPtr index_ptr = MassMatchFactory::getPrmDiagMassMatchPtr(proteo_ptrs,
                                                                    mng_ptr->max_proteoform_mass_,
                                                                    mng_ptr->filter_scale_);
                                                        
  std::string parameters = mng_ptr->getIndexFilePara();
  std::string dir_name = mng_ptr->prsm_para_ptr_->getOriDbName() + "_idx";
  std::string file_name = mng_ptr->multi_ptm_file_vec_[0] + parameters + block_str;
  index_ptr->serializeMassMatch(file_name, dir_name);
}

}

} /* namespace toppic */
