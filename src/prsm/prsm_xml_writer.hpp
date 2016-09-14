// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_PRSM_PRSM_XML_WRITER_HPP_
#define PROT_PRSM_PRSM_XML_WRITER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

class PrsmXmlWriter {
public:
 PrsmXmlWriter(const std::string &file_name);

 void close();

 void write(PrsmStrPtr prsm_str_ptr);
 void writeVector(const PrsmStrPtrVec &prsm_str_ptr_vec);

 void write(PrsmPtr prsm_ptr);
 void writeVector(const PrsmPtrVec &prsm_ptrs);
 void writeVector2D(const PrsmPtrVec2D &prsm_ptrs);
 void writeVector3D(const PrsmPtrVec3D &prsm_ptrs);
private:
  //XmlDOMDocument* doc_;
  std::ofstream file_;
};

typedef std::shared_ptr<PrsmXmlWriter> PrsmXmlWriterPtr;
typedef std::vector<PrsmXmlWriterPtr> PrsmXmlWriterPtrVec;

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
