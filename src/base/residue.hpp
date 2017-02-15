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


#ifndef PROT_BASE_RESIDUE_HPP_
#define PROT_BASE_RESIDUE_HPP_

#include <string>
#include <memory>
#include <map>

#include "base/acid.hpp"
#include "base/ptm.hpp"
#include "base/logger.hpp"

namespace prot {

class Residue;
typedef std::shared_ptr<Residue> ResiduePtr;

class Residue {
 public:
  Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr); 

  Residue(const std::string &acid_name, const std::string &abbr_name);

  Residue(xercesc::DOMElement* element); 

  /** Get amino acid. */
  AcidPtr getAcidPtr() {return acid_ptr_; }

  /** Get residue mass. */
  double getMass() { return mass_; }

  /** Get post-translational modification. */
  PtmPtr getPtmPtr() { return ptm_ptr_; }

  /**
   * Checks if the residue contains the same amino acid and ptm.
   */
  bool isSame(ResiduePtr residue_ptr) {
    return acid_ptr_ == residue_ptr->getAcidPtr() 
        && ptm_ptr_ == residue_ptr->getPtmPtr();
  }

  /** Get string representation */
  std::string toString(const std::string &delim_bgn, 
                       const std::string &delim_end);

  std::string toString() {return toString("[", "]");}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                 const std::string &element_name);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "residue";}

private:
  /** amino acid */
  AcidPtr acid_ptr_;
  /** post-translational modification */
  PtmPtr ptm_ptr_;
  /** residue mass */
  double mass_;
};

typedef std::vector<ResiduePtr> ResiduePtrVec;
typedef std::vector<ResiduePtrVec> ResiduePtrVec2D;

typedef std::vector<std::pair<std::string,std::string>> StringPairVec;

}

#endif
