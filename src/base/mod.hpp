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


#ifndef PROT_BASE_MOD_HPP_
#define PROT_BASE_MOD_HPP_

#include "base/residue.hpp"
#include "base/residue_base.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Mod;
typedef std::shared_ptr<Mod> ModPtr;

class Mod {
 public:
  Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr);

  Mod(xercesc::DOMElement* element); 

  ResiduePtr getOriResiduePtr() { return ori_residue_ptr_;}

  ResiduePtr getModResiduePtr() { return mod_residue_ptr_;}

  bool isSame(ModPtr mod_ptr) {
    return ori_residue_ptr_ == mod_ptr->getOriResiduePtr() 
        && mod_residue_ptr_ == mod_ptr->getModResiduePtr();
  }

  double getShift() {return mod_residue_ptr_->getMass() - ori_residue_ptr_->getMass();}

  void appendToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "mod";}

 private:
  ResiduePtr ori_residue_ptr_;
  ResiduePtr mod_residue_ptr_;
};

typedef std::vector<ModPtr> ModPtrVec;

}
#endif
