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


#ifndef PROT_BASE_ACID_HPP_
#define PROT_BASE_ACID_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class Acid {
 public:
  Acid(const std::string &name, const std::string &one_letter,
       const std::string &three_letter, const std::string &composition,
       double mono_mass, double avg_mass);

  explicit Acid(xercesc::DOMElement* element);

  // Get amino acid composition
  std::string getAcidComposition() {return composition_;}

  // Get average mass
  double getAvgMass() {return average_mass_;}

  // Get monoisotopic mass
  double getMonoMass() {return mono_mass_;}

  // Get amino acid name
  std::string getName() {return name_;}

  // Get amino acid one letter representation
  std::string getOneLetter() {return one_letter_;}

  // Get amino acid three letter representation
  std::string getThreeLetter() {return three_letter_;}

  void appendNameToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getNameFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "amino_acid";}

 private:
  // Name of amino acid
  std::string name_;
  // One letter representation
  std::string one_letter_;
  // Three letter representation
  std::string three_letter_;
  // amino acid chemical composition
  std::string composition_;
  // residue monoisotopic mass
  double mono_mass_;
  // residue average mass
  double average_mass_;
};

typedef std::shared_ptr<Acid> AcidPtr;
typedef std::vector<AcidPtr> AcidPtrVec;

}  // namespace prot

#endif

