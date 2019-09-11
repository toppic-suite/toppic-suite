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

#include "common/xml/xml_dom_util.hpp"
#include "prsmview/anno_file_list.hpp"

namespace toppic {

xercesc::DOMElement* AnnoFileList::geneFileList(XmlDOMDocument* xml_doc) {
  xercesc::DOMElement* element = xml_doc->createElement("file_list");
  for (size_t i = 0; i < file_list_.size(); i++) {
    xercesc::DOMElement* file = xml_doc->createElement("file");
    xml_doc->addElement(file, "xml", file_list_[i][0].c_str());
    xml_doc->addElement(file, "json", file_list_[i][1].c_str());
    element->appendChild(file);
  }
  return element;
}

std::vector<std::vector<std::string>> AnnoFileList::readFromXml(const std::string &file_name) {
  std::vector<std::vector<std::string>> file_list;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name.c_str());
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      int file_num = xml_dom_util::getChildCount(root, "file");
      for (int i = 0; i < file_num; i++) {
        xercesc::DOMElement* file_element = xml_dom_util::getChildElement(root, "file", i);
        std::vector<std::string> file_info;
        file_info.push_back(xml_dom_util::getChildValue(file_element, "xml", 0));
        file_info.push_back(xml_dom_util::getChildValue(file_element, "json", 0));
        file_list.push_back(file_info);
      }
    }
    delete doc;
  }
  return file_list;
}

}  // namespace toppic
