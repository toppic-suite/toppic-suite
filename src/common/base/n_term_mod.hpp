//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_BASE_N_TERM_MOD_HPP_
#define TOPPIC_COMMON_BASE_N_TERM_MOD_HPP_

#include "common/xml/xml_dom_element.hpp"
#include "common/base/residue.hpp"

namespace toppic {

class XmlDOMDocument;

class NTermMod;
typedef std::shared_ptr<NTermMod> NTermModPtr;

class NTermMod {
 public:
  NTermMod(ResiduePtr residue_ptr, PtmPtr ori_ptm_ptr, PtmPtr mod_ptm_ptr);

  explicit NTermMod(XmlDOMElement* element);

  ResiduePtr getResiduePtr() { return residue_ptr_;}

  PtmPtr getOriPtmPtr() { return ori_ptm_ptr_;}

  PtmPtr getModPtmPtr() { return mod_ptm_ptr_;}

  bool isSame(NTermModPtr mod_ptr);

  double getMass();

  double getShift();

  void appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getXmlElementName() {return "n_term_mod";}

 private:
  ResiduePtr residue_ptr_;
  PtmPtr ori_ptm_ptr_;
  PtmPtr mod_ptm_ptr_;
};

typedef std::vector<NTermModPtr> NTermModPtrVec;

typedef std::vector<NTermModPtrVec> NTermModPtrVec2D;

}  // namespace toppic

#endif
