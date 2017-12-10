// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef PROT_TAG_LONG_FILTER_MNG_HPP
#define PROT_TAG_LONG_FILTER_MNG_HPP

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class TagFilterMng {
 public:
  TagFilterMng(PrsmParaPtr prsm_para_ptr,
               const std::string & output_file_ext):
      prsm_para_ptr_(prsm_para_ptr),
      output_file_ext_(output_file_ext) {}

  size_t top_num_ = 20;

  size_t L = 150;

  size_t tag_min_len_ = 4;

  PrsmParaPtr prsm_para_ptr_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<TagFilterMng> TagFilterMngPtr;

}  // namespace prot
#endif
