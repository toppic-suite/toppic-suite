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

#ifndef TOPPIC_PRSM_PRSM_XML_WRITER_SET_HPP_
#define TOPPIC_PRSM_PRSM_XML_WRITER_SET_HPP_

#include <vector>

#include "prsm/prsm_xml_writer.hpp"

namespace toppic {

class PrsmXmlWriterSet {
 public:
  PrsmXmlWriterSet(const std::string & output_file_name, int num_mass_shift);

  PrsmXmlWriterPtr getCompleteWriterPtr(int n_shift) {return complete_writer_ptrs_[n_shift-2];}

  PrsmXmlWriterPtr getPrefixWriterPtr(int n_shift) {return prefix_writer_ptrs_[n_shift-2];}

  PrsmXmlWriterPtr getSuffixWriterPtr(int n_shift) {return suffix_writer_ptrs_[n_shift-2];}

  PrsmXmlWriterPtr getInternalWriterPtr(int n_shift) {return internal_writer_ptrs_[n_shift-2];}

  void close();

 private:
  std::vector<PrsmXmlWriterPtr> complete_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> prefix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> suffix_writer_ptrs_;
  std::vector<PrsmXmlWriterPtr> internal_writer_ptrs_;
  // maximum number of mass shifts;
  int n_mass_shift_;

};

typedef std::shared_ptr<PrsmXmlWriterSet> PrsmXmlWriterSetPtr;
typedef std::vector<PrsmXmlWriterSetPtr> PrsmXmlWriterSetPtrVec;

} /* namespace toppic */

#endif 
