// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_BASE_PROT_MOD_HPP_
#define PROT_BASE_PROT_MOD_HPP_

#include "base/ptm_base.hpp"
#include "base/mod.hpp"
#include "base/trunc.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class ProtMod {
 public:
  ProtMod(const std::string &name, const std::string &type,
          TruncPtr trunc_ptr, ModPtr mod_ptr);

  ProtMod(xercesc::DOMElement* element); 

  const std::string& getName() { return name_;};

  const std::string& getType() { return type_;};

  TruncPtr getTruncPtr() { return trunc_ptr_;}

  ModPtr getModPtr() { return mod_ptr_;}

  int getModPos() {return mod_pos_;}

  double getProtShift() { return prot_shift_;}

  double getPepShift() { return pep_shift_;}

  bool isAcetylation();

  void appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "prot_mod";}

  static std::string getNameFromXml(xercesc::DOMElement * element);

 private:
  std::string name_;
  std::string type_;
  TruncPtr trunc_ptr_;
  ModPtr mod_ptr_;
  int mod_pos_;
  double prot_shift_;
  double pep_shift_;
};

typedef std::shared_ptr<ProtMod> ProtModPtr;
typedef std::vector<ProtModPtr> ProtModPtrVec;

}
#endif
