//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEQ_ALTER_TYPE_HPP_
#define TOPPIC_SEQ_ALTER_TYPE_HPP_

#include <memory>
#include <vector>
#include <string>

#include "common/xml/xml_dom_element.hpp"

namespace toppic {

class XmlDOMDocument;

class AlterType;
typedef std::shared_ptr<AlterType> AlterTypePtr;

class AlterType {
 public:
  static const AlterTypePtr INPUT;

  static const AlterTypePtr FIXED;

  static const AlterTypePtr PROTEIN_VARIABLE;

  static const AlterTypePtr VARIABLE;

  static const AlterTypePtr UNEXPECTED;

  AlterType(int id, std::string name);

  int getId() {return id_;}

  std::string getName() {return name_;}

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static AlterTypePtr getTypePtrFromXml(XmlDOMElement* element);

  static std::string getXmlElementName() {return "alter_type";}

 private:
  int id_;
  std::string name_;
};

typedef std::vector<AlterTypePtr> AlterTypePtrVec;

}  // namespace toppic

#endif

