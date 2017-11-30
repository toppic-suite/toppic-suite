//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_PRSM_PRSM_READER_HPP_
#define PROT_PRSM_PRSM_READER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmReader {
 public:
  explicit PrsmReader(const std::string &file_name);

  std::vector<std::string> readOnePrsmLines();

  PrsmStrPtr readOnePrsmStr();

  PrsmPtr readOnePrsm(FastaIndexReaderPtr reader_ptr, const ModPtrVec fix_mod_list);

  void close();

  static PrsmStrPtrVec readAllPrsmStrs(const std::string &input_file_name);

  static PrsmStrPtrVec readAllPrsmStrsMatchSeq(const std::string &input_file_name,
                                               FastaIndexReaderPtr fasta_reader_ptr,
                                               const ModPtrVec fix_mod_list);

  static PrsmPtrVec readAllPrsms(const std::string &prsm_file_name,
                                 const std::string &db_file_name,
                                 const ModPtrVec  &fix_mod_list);

 private:
  std::ifstream input_;
};

typedef std::shared_ptr<PrsmReader> PrsmReaderPtr;
typedef std::vector<PrsmReaderPtr> PrsmReaderPtrVec;

} /* namespace prot */

#endif /* PROT_PRSM_READER_HPP_ */
