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

#ifndef TOPPIC_PRSM_PRSM_TOP_SELECTOR_HPP_
#define TOPPIC_PRSM_PRSM_TOP_SELECTOR_HPP_

#include <memory>
#include <string>

namespace toppic {

class PrsmTopSelector {
 public:
  PrsmTopSelector(const std::string &db_file_name,
                  const std::string &spec_file_name,
                  const std::string &in_file_ext, 
                  const std::string &out_file_ext, int n_top);

  void process();
 private:
  std::string spec_file_name_;

  std::string db_file_name_;

  std::string input_file_ext_;

  std::string output_file_ext_;

  int n_top_;
};

typedef std::shared_ptr<PrsmTopSelector> PrsmTopSelectorPtr;

} /* namespace toppic */

#endif /* PRSM_SELECTOR_HPP_ */
