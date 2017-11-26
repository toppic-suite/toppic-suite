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


#ifndef PROT_BASE_CHANGE_TYPE_HPP_
#define PROT_BASE_CHANGE_TYPE_HPP_

#include <memory>
#include <vector>
#include <string>

#include "base/xml_dom_document.hpp"

namespace prot {

class ChangeType;
typedef std::shared_ptr<ChangeType> ChangeTypePtr;

class ChangeType {
 public:
  static const ChangeTypePtr INPUT;

  static const ChangeTypePtr FIXED;

  static const ChangeTypePtr PROTEIN_VARIABLE;

  static const ChangeTypePtr VARIABLE;

  static const ChangeTypePtr UNEXPECTED;

  ChangeType(int id, std::string name): id_(id), name_(name) {}

  int getId() {return id_;}

  std::string getName() {return name_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static ChangeTypePtr getChangeTypePtrFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "change_type";}

 private:
  int id_;
  std::string name_;
};

typedef std::vector<ChangeTypePtr> ChangeTypePtrVec;

}  // namespace prot

#endif

