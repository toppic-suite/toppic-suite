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

#ifndef TOPPIC_FILTER_ZERO_PTM_ZERO_PTM_INDEX_FILE_HPP_
#define TOPPIC_FILTER_ZERO_PTM_ZERO_PTM_INDEX_FILE_HPP_

#include "filter/massmatch/mass_match.hpp"
#include "filter/zeroptm/zero_ptm_filter_mng.hpp"

namespace toppic {

class ZeroPtmIndexFile {
 public:
  ZeroPtmIndexFile(const ProteoformPtrVec &proteo_ptrs, ZeroPtmFilterMngPtr
                    mng_ptr, std::vector<std::string> file_vec);

 private:
  ZeroPtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;

  MassMatchPtr diag_index_ptr_;
  MassMatchPtr rev_diag_index_ptr_;
  MassMatchPtr term_index_ptr_;
  MassMatchPtr rev_term_index_ptr_;
  
  std::string block_str; //which db_block it is. Used when merging the output files by order.
};

typedef std::shared_ptr<ZeroPtmIndexFile> ZeroPtmIndexFilePtr;
} /* namespace toppic */

#endif 
