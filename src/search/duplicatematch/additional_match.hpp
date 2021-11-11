//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEARCH_DUPLICATEMATCH_ADDITIONAL_MATCH_HPP_
#define TOPPIC_SEARCH_DUPLICATEMATCH_ADDITIONAL_MATCH_HPP_

#include "prsm/prsm.hpp"

namespace toppic {

class AdditionalMatch {
 public:
  AdditionalMatch(std::string db_file_name, std::vector<PrsmPtr> &prsm_ptr_vec, SpParaPtr sp_para_ptr):
      db_file_name_(db_file_name),
      prsm_ptr_vec_(prsm_ptr_vec),
      sp_para_ptr_(sp_para_ptr) {}

  std::string trimSequence(std::string raw_prot);
  void process(PrsmPtr prsm_ptr_);

 private:
  std::string db_file_name_;
  std::vector<PrsmPtr> &prsm_ptr_vec_;
  SpParaPtr sp_para_ptr_;
  int prsm_cnt_ = 0;
};

typedef std::shared_ptr<AdditionalMatch> AdditionalMatchPtr;

} //namespace toppic

#endif