//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_BASE_PTM_HPP_
#define TOPPIC_COMMON_BASE_PTM_HPP_

#include <string>
#include <vector>
#include <memory>

#include "common/xml/xml_dom_element.hpp"

namespace toppic {

class XmlDOMDocument;

class Ptm;
using PtmPtr = std::shared_ptr<Ptm>;

class Ptm {
 public:
  Ptm(const std::string &name, const std::string &abbr_name,
      double mono_mass, int unimod_id);

  explicit Ptm(XmlDOMElement* element);

  const std::string& getName() const {return name_;}

  const std::string& getAbbrName() const {return abbr_name_;}

  // Get monoisotopic mass.
  double getMonoMass() const {return mono_mass_;}

  int getUnimodId() const {return unimod_id_;}

  void appendAbbrNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) const;

  void appendAbbrNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent, const std::string &element_name) const;

  // Add mass for visualization
  void appendAbbrNameMassToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) const;

  static std::string getAbbrNameFromXml(XmlDOMElement * element);

  static std::string getXmlElementName() {return "ptm";}

  // comparison function
  static bool cmpMassInc(const PtmPtr &a, const PtmPtr &b) {
    return a->getMonoMass() < b->getMonoMass();
  }

  bool isSame(PtmPtr ptm_ptr) const {return abbr_name_ == ptm_ptr->getAbbrName();}

 private:
  /* Full name */
  std::string name_;
  // abbreviation name
  std::string abbr_name_;
  /* monoisotopic mass */
  double mono_mass_;
  // unimod id
  int unimod_id_;
};

using PtmPtrVec = std::vector<PtmPtr>;
using PtmPair = std::pair<PtmPtr, PtmPtr>;
using PtmPairVec = std::vector<PtmPair>;

}  // namespace toppic

#endif
