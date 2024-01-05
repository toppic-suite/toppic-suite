//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_PRSM_SIMPLE_PRSM_XML_WRITER_HPP_
#define TOPPIC_PRSM_SIMPLE_PRSM_XML_WRITER_HPP_

#include <fstream>

#include "common/xml/xml_dom_impl.hpp"

#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace toppic {

class SimplePrsmXmlWriter {
 public:
  SimplePrsmXmlWriter();

  explicit SimplePrsmXmlWriter(const std::string &file_name);

  ~SimplePrsmXmlWriter();

  void close();

  void write(SimplePrsmStrPtr prsm_str_ptr);

  void write(const SimplePrsmPtrVec &simple_prsm_ptrs);

  void write(SimplePrsmPtr simple_prsm_ptrs);

 private:
  std::ofstream file_;

  XmlDOMImpl* impl_; 

  xercesc::DOMLSSerializer* serializer_;

  std::string file_name_;
};

typedef std::shared_ptr<SimplePrsmXmlWriter>   SimplePrsmXmlWriterPtr;
typedef std::vector<SimplePrsmXmlWriterPtr>    SimplePrsmXmlWriterPtrVec;

}  // namespace toppic

#endif /* SIMPLE_PRSM_WRITER_HPP_ */
