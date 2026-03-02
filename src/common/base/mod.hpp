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

#ifndef TOPPIC_COMMON_BASE_MOD_HPP_
#define TOPPIC_COMMON_BASE_MOD_HPP_

#include "common/xml/xml_dom_element.hpp"
#include "common/base/residue.hpp"
#include "common/base/mod_type.hpp"

namespace toppic {

class XmlDOMDocument;

class Mod;
typedef std::shared_ptr<Mod> ModPtr;

class Mod {
 public:
  Mod(ResiduePtr ori_residue_ptr, 
      ResiduePtr mod_residue_ptr, 
      ModTypePtr mod_type_ptr);

  explicit Mod(XmlDOMElement* element);

  ResiduePtr getOriResiduePtr() const { return ori_residue_ptr_;}

  ResiduePtr getModResiduePtr() const { return mod_residue_ptr_;}

  ModTypePtr getModTypePtr() const { return mod_type_ptr_;}

  bool isSame(ModPtr mod_ptr) const;

  double getReplaceShift() const;

  double getShift() const;

  void appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) const;

  static std::string getXmlElementName() {return "mod";}

 private:
  ResiduePtr ori_residue_ptr_;
  ResiduePtr mod_residue_ptr_;
  ModTypePtr mod_type_ptr_;
};

typedef std::vector<ModPtr> ModPtrVec;

typedef std::vector<ModPtrVec> ModPtrVec2D;

}  // namespace toppic

#endif
