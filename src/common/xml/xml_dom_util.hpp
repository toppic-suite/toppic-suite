//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_COMMON_XML_XML_DOM_UTIL_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_UTIL_HPP_

#include <string>

#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>

namespace toppic {

namespace xml_dom_util {

xercesc::DOMNodeList* getChildElements(xercesc::DOMElement *parent, const char* tag);

xercesc::DOMElement* getChildElement(xercesc::DOMElement* parent, const char* tag, int index);

std::string getChildValue(xercesc::DOMElement* parent, const char* child_tag, int index);

double getScientificChildValue(xercesc::DOMElement* parent, const char* child_tag, int index);

double getDoubleChildValue(xercesc::DOMElement* parent, const char* child_tag, int index);

int getIntChildValue(xercesc::DOMElement* parent, const char* child_tag, int index);

bool getBoolChildValue(xercesc::DOMElement* parent, const char* child_tag, int index);

int getChildCount(xercesc::DOMElement * parent, const char* child_tag);

std::string getAttributeValue(xercesc::DOMElement * parent, const char* attribute_tag);

std::string writeToString(xercesc::DOMLSSerializer* serializer, xercesc::DOMNode *node);

void writeToStreamByRemovingDoubleLF(std::ofstream & file, std::string &str);

}  // namespace xml_dom_util

}  // namespace toppic
#endif
