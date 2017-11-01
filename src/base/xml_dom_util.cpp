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


#include <string>
#include <sstream>
#include <exception>
#include <algorithm>

#include "base/logger.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

xercesc::DOMElement* XmlDomUtil::getChildElement(xercesc::DOMElement *parent, 
                                                 const char* tag, int index) { 
  xercesc::DOMNodeList* list = parent->getElementsByTagName(X(tag));
  xercesc::DOMElement* element = 
      dynamic_cast<xercesc::DOMElement*>(list->item(index));
  if (element == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << tag << " return null";
    LOG_WARN(stream.str());
    throw stream.str();
  }
  return element;
}

std::string XmlDomUtil::getChildValue(xercesc::DOMElement* parent,  
                                      const char* child_tag, int i) {
  xercesc::DOMNodeList* node_list;
  node_list= parent->getElementsByTagName(X(child_tag));
  if (node_list == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << child_tag << " return null";
    LOG_WARN( stream.str());
    throw stream.str();
  }
  xercesc::DOMElement* child = 
  dynamic_cast<xercesc::DOMElement*>(node_list->item(i));
  if (child == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << child_tag << " return null";
    LOG_WARN( stream.str());
    throw stream.str();
  }

  std::string value;
  if (child) {
    value = Y(child->getTextContent());
  }
  else {
    value = "";
  }
  return value;
}

double XmlDomUtil::getDoubleChildValue(xercesc::DOMElement* parent,  
                                       const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  //LOG_DEBUG("tag " << child_tag << "double value " << value);
  return std::stod(value);
}

int XmlDomUtil::getIntChildValue(xercesc::DOMElement* parent,  
                                 const char* child_tag, int i) {
  try {
    std::string value = getChildValue(parent, child_tag, i);
    return std::stoi(value);
  }
  catch (std::string s) {
    return 0;
  }
}

bool XmlDomUtil::getBoolChildValue(xercesc::DOMElement* parent,
                                   const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  if(value == "true") {
    return true;
  }
  return false;
}

int XmlDomUtil::getChildCount(xercesc::DOMElement* parent, 
                              const char* child_tag) {
  xercesc::DOMNodeList* childList = parent->getElementsByTagName(X(child_tag));
  return (int)childList->getLength();
}

std::string XmlDomUtil::getAttributeValue(xercesc::DOMElement* element, 
                                          const char* attribute_tag) {
  std::string value = Y(element->getAttribute(X(attribute_tag)));
  return value;
}

std::string XmlDomUtil::writeToString(xercesc::DOMLSSerializer* serializer, 
                                      xercesc::DOMNode *node) {
  XMLCh* ch = serializer->writeToString(node, 0);
  std::string result = Y(ch);
  xercesc::XMLString::release(&ch);
  return result;
}

void XmlDomUtil::writeToStreamByRemovingDoubleLF(std::ofstream &file, 
                                                 std::string &str) {
  int pos = 0;
  std::size_t found = str.find("\n\n", pos);
  while (found != std::string::npos) {
    std::string sub = str.substr(pos, found - pos);
    file << sub << std::endl;
    pos = found + 2;
    found = str.find("\n\n", pos);
  }
  if (pos < (int)str.length()) {
    std::string sub = str.substr(pos);
    file << sub << std::endl; 
  }
}

}
