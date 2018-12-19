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


#ifndef TOPPIC_BASE_PTM_HPP_
#define TOPPIC_BASE_PTM_HPP_

#include <string>
#include <vector>
#include <memory>

#include <xercesc/dom/DOMElement.hpp>

namespace toppic {

class XmlDOMDocument;

class Ptm;
typedef std::shared_ptr<Ptm> PtmPtr;

class Ptm {
 public:
  Ptm(const std::string &name, const std::string &abbr_name,
      double mono_mass, int unimod_id);

  explicit Ptm(xercesc::DOMElement* element);

  const std::string& getName() {return name_;}

  const std::string& getAbbrName() {return abbr_name_;}

  /* Get  monoisotopic mass. */
  double getMonoMass() {return mono_mass_;}

  int getUnimodId() {return unimod_id_;}

  void appendAbbrNameToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getAbbrNameFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "ptm";}

  // comparison function
  static bool cmpMassInc(const PtmPtr &a, const PtmPtr &b) {
    return a->getMonoMass() < b->getMonoMass();
  }

  bool isSame(PtmPtr ptm_ptr) {return abbr_name_ == ptm_ptr->getAbbrName();}

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
typedef std::pair<PtmPtr, PtmPtr> PtmPair;
typedef std::vector<PtmPair> PtmPairVec;
}  // namespace toppic

#endif
