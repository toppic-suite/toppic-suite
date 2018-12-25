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


#ifndef PROT_PRSM_SIMPLE_PRSM_READER_HPP_
#define PROT_PRSM_SIMPLE_PRSM_READER_HPP_

#include <iostream>
#include <fstream>

#include "common/xml/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace toppic {

class SimplePrsmReader {
 public:
  SimplePrsmReader(const std::string &file_name);

  std::vector<std::string> readOnePrsmLines();

  SimplePrsmStrPtr readOnePrsmStr();

  SimplePrsmPtr readOnePrsm();

  void close();

  static SimplePrsmPtrVec readSimplePrsms(const std::string &file_name);

 private:
  std::ifstream input_;
};

typedef std::shared_ptr<SimplePrsmReader> SimplePrsmReaderPtr;
typedef std::vector<SimplePrsmReaderPtr> SimplePrsmReaderPtrVec;

} /* namespace toppic */

#endif /* PROT_SIMPLE_PRSM_READER_HPP_ */
