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

#include "common/util/logger.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_reader_util.hpp"

namespace toppic {

namespace prsm_reader_util {

PrsmStrPtrVec readAllPrsmStrs(const std::string &input_file_name) {
  PrsmReader reader(input_file_name);
  PrsmStrPtrVec prsm_str_ptrs;
  PrsmStrPtr prsm_str_ptr = reader.readOnePrsmStr();
  while (prsm_str_ptr != nullptr) {
    prsm_str_ptrs.push_back(prsm_str_ptr);
    prsm_str_ptr = reader.readOnePrsmStr();
  }
  reader.close();
  return prsm_str_ptrs;
}

PrsmStrPtrVec readAllPrsmStrsMatchSeq(const std::string &input_file_name) {
  PrsmReaderPtr str_reader = std::make_shared<PrsmReader>(input_file_name);
  PrsmStrPtrVec prsm_str_ptrs;
  PrsmStrPtr prsm_str_ptr = str_reader->readOnePrsmStr();
  while (prsm_str_ptr != nullptr) {
    prsm_str_ptrs.push_back(prsm_str_ptr);
    prsm_str_ptr = str_reader->readOnePrsmStr();
  }
  str_reader->close();
  return prsm_str_ptrs;
}

PrsmPtrVec readAllPrsms(const std::string &prsm_file_name,
                        const std::string &db_file_name,
                        const ModPtrVec  &fix_mod_list) {
  FastaIndexReaderPtr fasta_reader_ptr 
      = std::make_shared<FastaIndexReader>(db_file_name);
  PrsmReader reader(prsm_file_name);
  PrsmPtrVec prsm_ptrs;
  PrsmPtr prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  while (prsm_ptr != nullptr) {
    prsm_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  }
  reader.close();
  return prsm_ptrs;
}

PrsmPtrVec readAllPrsms(const std::string &prsm_file_name,
                        FastaIndexReaderPtr fasta_reader_ptr, 
                        const ModPtrVec  &fix_mod_list) {
  PrsmReader reader(prsm_file_name);
  PrsmPtrVec prsm_ptrs;
  PrsmPtr prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  while (prsm_ptr != nullptr) {
    prsm_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  }
  reader.close();
  return prsm_ptrs;
}


}

}  // namespace toppic
