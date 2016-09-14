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


#ifndef PROT_BASE_PTM_HPP_
#define PROT_BASE_PTM_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/acid.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;

class Ptm {
 public:
  Ptm(const std::string &name, const std::string &abbr_name,
      double mono_mass, int unimod_id, 
      const std::string &n_term_residue_str,
      const std::string &c_term_residue_str, 
      const std::string &anywhere_residue_str);

  Ptm(xercesc::DOMElement* element); 

  const std::string& getName() {return name_;}

  const std::string& getAbbrName() {return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  int getUnimodId() {return unimod_id_;}

  void appendAbbrNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getAbbrNameFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "ptm";}

  // comparison function
  static bool cmpMassInc(const PtmPtr &a, const PtmPtr &b) {
    return a->getMonoMass() < b->getMonoMass();
  }

  bool isSame(PtmPtr ptm_ptr) { return abbr_name_ == ptm_ptr->getAbbrName();}

 private:
  /* Full name */
  std::string name_;
  // abbrevation name 
  std::string abbr_name_;
  /* monoisotopic mass */
  double mono_mass_;
  // unimod id
  int unimod_id_;
};

typedef std::vector<PtmPtr> PtmPtrVec;

}

#endif
