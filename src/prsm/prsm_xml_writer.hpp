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


#ifndef PROT_PRSM_PRSM_XML_WRITER_HPP_
#define PROT_PRSM_PRSM_XML_WRITER_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmXmlWriter {
 public:
  explicit PrsmXmlWriter(const std::string &file_name);

  void close();

  void write(PrsmStrPtr prsm_str_ptr);

  void writeVector(const PrsmStrPtrVec &prsm_str_ptr_vec);

  void write(PrsmPtr prsm_ptr);

  void writeVector(const PrsmPtrVec &prsm_ptrs);

  void writeVector2D(const PrsmPtrVec2D &prsm_ptrs);

  void writeVector3D(const PrsmPtrVec3D &prsm_ptrs);

 private:
  // XmlDOMDocument* doc_;
  std::ofstream file_;

  std::string file_name_;
};

typedef std::shared_ptr<PrsmXmlWriter> PrsmXmlWriterPtr;
typedef std::vector<PrsmXmlWriterPtr> PrsmXmlWriterPtrVec;

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
