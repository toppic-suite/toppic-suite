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


#include <set>
#include <utility>
#include <limits>
#include <string>
#include <algorithm>
#include <vector>

#include "boost/algorithm/string.hpp"

#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"
#include "base/fasta_reader.hpp"
#include "base/proteoform_factory.hpp"
#include "spec/peak.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_util.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_prsm.hpp"
#include "prsmview/anno_view.hpp"

namespace prot {

xercesc::DOMElement* AnnoView::geneFileList(XmlDOMDocument* xml_doc) {
  xercesc::DOMElement* element = xml_doc->createElement("file_list");
  for (size_t i = 0; i < file_list_.size(); i++) {
    xercesc::DOMElement* file = xml_doc->createElement("file");
    xml_doc->addElement(file, "xml", file_list_[i][0].c_str());
    xml_doc->addElement(file, "xsl", file_list_[i][1].c_str());
    xml_doc->addElement(file, "html", file_list_[i][2].c_str());
    element->appendChild(file);
  }
  return element;
}

std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name) {
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
        file_info.push_back(xml_dom_util::getChildValue(file_element, "xsl", 0));
        file_info.push_back(xml_dom_util::getChildValue(file_element, "html", 0));
        file_list.push_back(file_info);
      }
    }
    delete doc;
  }
  return file_list;
}

xercesc::DOMElement* proteoformToXml(XmlDOMDocument* xml_doc, const PrsmPtrVec &prsm_ptrs,
                                     PrsmViewMngPtr mng_ptr, bool detail, bool add_ms) {
  xercesc::DOMElement* proteoform_element = xml_doc->createElement("compatible_proteoform");
  std::string str = string_util::convertToString(prsm_ptrs[0]->getProteoformPtr()->getProtId());
  xml_doc->addElement(proteoform_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(proteoform_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(proteoform_element, "sequence_description", str.c_str());
  str = string_util::convertToString(prsm_ptrs[0]->getProteoformPtr()->getProteoClusterId());
  xml_doc->addElement(proteoform_element, "proteoform_id", str.c_str());
  int count = prsm_ptrs.size();
  str = string_util::convertToString(count);
  xml_doc->addElement(proteoform_element, "prsm_number", str.c_str());
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    proteoform_element->appendChild(geneAnnoPrsm(xml_doc, prsm_ptrs[i], mng_ptr, detail, add_ms));
  }
  return proteoform_element;
}

xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  int prot_id,
                                  const std::vector<int> &cluster_ids,
                                  PrsmViewMngPtr mng_ptr,
                                  bool detail, bool add_ms) {
  xercesc::DOMElement* prot_element = xml_doc->createElement("protein");
  std::string str = string_util::convertToString(prot_id);
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  int count = cluster_ids.size();
  str = string_util::convertToString(count);
  xml_doc->addElement(prot_element, "compatible_proteoform_number", str.c_str());
  for (size_t i = 0; i < cluster_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, cluster_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(), Prsm::cmpEValueInc);
    prot_element->appendChild(proteoformToXml(xml_doc, select_prsm_ptrs, mng_ptr, detail, add_ms));
  }
  return prot_element;
}

void writeProteinToXml(XmlWriterPtr xml_writer,
                       const PrsmPtrVec &prsm_ptrs,
                       int prot_id,
                       const std::vector<int> &species_ids,
                       PrsmViewMngPtr mng_ptr,
                       bool detail, bool add_ms) {
  xml_writer->write_str("<protein>");
  xml_writer->write_str("<sequence_id>" + std::to_string(prot_id) + "</sequence_id>");
  xml_writer->write_str("<sequence_name>" + prsm_ptrs[0]->getProteoformPtr()->getSeqName() + "</sequence_name>");
  xml_writer->write_str("<sequence_description>" + prsm_ptrs[0]->getProteoformPtr()->getSeqDesc() + "</sequence_description>");
  xml_writer->write_str("<compatible_proteoform_number>" + std::to_string(species_ids.size()) + "</compatible_proteoform_number>");
  for (size_t i = 0; i < species_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, species_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(), Prsm::cmpEValueInc);
    xml_writer->write(proteoformToXml(xml_writer->getDoc(), select_prsm_ptrs, mng_ptr, detail, add_ms));
  }
  xml_writer->write_str("</protein>");
}

}  // namespace prot
