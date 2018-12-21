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

#ifndef TOPPIC_BASE_AMINO_ACID_HPP_
#define TOPPIC_BASE_AMINO_ACID_HPP_

#include <string>
#include <memory>
#include <vector>

#include "xml/xml_dom_element.hpp"

namespace toppic {

class XmlDOMDocument;

class AminoAcid {
 public:
  AminoAcid(const std::string &name, const std::string &one_letter,
            const std::string &three_letter, const std::string &composition,
            double mono_mass, double avg_mass);

  explicit AminoAcid(XmlDOMElement* element);

  // Get amino acid composition
  std::string getComposition() {return composition_;}
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

  static std::string getXmlElementName() {return "amino_acid";}

  void appendNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getNameFromXml(XmlDOMElement * element);

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

typedef std::shared_ptr<AminoAcid> AminoAcidPtr;
typedef std::vector<AminoAcidPtr> AminoAcidPtrVec;

}  // namespace toppic

#endif

