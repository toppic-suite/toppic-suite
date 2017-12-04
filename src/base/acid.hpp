//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
       double mono_mass, double avg_mass):
      name_(name),
      one_letter_(one_letter),
      three_letter_(three_letter),
      composition_(composition),
      mono_mass_(mono_mass),
      average_mass_(avg_mass) {}

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

