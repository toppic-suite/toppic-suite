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

#ifndef TOPPIC_PRSM_SIMPLE_PRSM_XML_WRITER_SET_HPP_
#define TOPPIC_PRSM_SIMPLE_PRSM_XML_WRITER_SET_HPP_

#include "prsm/simple_prsm_xml_writer.hpp"

namespace toppic {

class SimplePrsmXmlWriterSet {
 public:
  SimplePrsmXmlWriterSet(const std::string & output_file_name);

  SimplePrsmXmlWriterPtr getCompleteWriterPtr() {return complete_writer_ptr_;}

  SimplePrsmXmlWriterPtr getPrefixWriterPtr() {return prefix_writer_ptr_;}

  SimplePrsmXmlWriterPtr getSuffixWriterPtr() {return suffix_writer_ptr_;}

  SimplePrsmXmlWriterPtr getInternalWriterPtr() {return internal_writer_ptr_;}

  void close();

 private:
  SimplePrsmXmlWriterPtr complete_writer_ptr_;
  SimplePrsmXmlWriterPtr prefix_writer_ptr_;
  SimplePrsmXmlWriterPtr suffix_writer_ptr_;
  SimplePrsmXmlWriterPtr internal_writer_ptr_;
};

typedef std::shared_ptr<SimplePrsmXmlWriterSet> SimplePrsmXmlWriterSetPtr;

} /* namespace toppic */

#endif 
