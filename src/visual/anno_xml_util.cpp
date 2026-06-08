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

#include "visual/anno_xml_util.hpp"

#include <algorithm>

#include "prsm/prsm_util.hpp"
#include "visual/anno_prsm.hpp"

namespace toppic {

namespace anno_xml_util {

XmlDOMElement geneXmlForProteoform(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                          const PrsmPtrVec &prsm_ptrs,
                                          const PrsmViewMngPtr &mng_ptr,
                                          bool detail, bool add_ms) {
  XmlDOMElement proteoform_element = xml_doc->addElement(parent, "compatible_proteoform");
  std::string str = std::to_string(prsm_ptrs[0]->getProteoformPtr()->getProtId());
  xml_doc->addElement(proteoform_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(proteoform_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(proteoform_element, "sequence_description", str.c_str());
  str = std::to_string(prsm_ptrs[0]->getProteoformPtr()->getProteoClusterId());
  xml_doc->addElement(proteoform_element, "proteoform_id", str.c_str());
  int count = prsm_ptrs.size();
  str = std::to_string(count);
  xml_doc->addElement(proteoform_element, "prsm_number", str.c_str());
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    anno_prsm::geneAnnoPrsm(xml_doc, proteoform_element, prsm_ptrs[i], mng_ptr, detail, add_ms);
  }
  return proteoform_element;
}

void writeProteinToXml(const XmlWriterPtr &xml_writer,
                       const PrsmPtrVec &prsm_ptrs,
                       int prot_id,
                       const std::vector<int> &species_ids,
                       const PrsmViewMngPtr &mng_ptr,
                       bool detail, bool add_ms) {
  xml_writer->writeStr("<protein>");
  xml_writer->writeStr("<sequence_id>" + std::to_string(prot_id) + "</sequence_id>");
  xml_writer->writeStr("<sequence_name>" 
                        + prsm_ptrs[0]->getProteoformPtr()->getSeqName() + "</sequence_name>");
  xml_writer->writeStr("<sequence_description>" 
                        + prsm_ptrs[0]->getProteoformPtr()->getSeqDesc() 
                        + "</sequence_description>");
  xml_writer->writeStr("<compatible_proteoform_number>" 
                        + std::to_string(species_ids.size()) 
                        + "</compatible_proteoform_number>");
  for (size_t i = 0; i < species_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, species_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(),
              Prsm::cmpEValueIncProtInc);
    xml_writer->writeAndRelease(geneXmlForProteoform(xml_writer->getDoc(),
                                                     xml_writer->getDoc()->getDocumentElement(),
                                                     select_prsm_ptrs, mng_ptr, detail, add_ms));
  }
  xml_writer->writeStr("</protein>");
}

XmlDOMElement geneXmlForProteinList(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                           const PrsmPtrVec &prsm_ptrs,
                                           int prot_id,
                                           const std::vector<int> &cluster_ids,
                                           const PrsmViewMngPtr &mng_ptr,
                                           bool detail, bool add_ms) {
  XmlDOMElement prot_element = xml_doc->addElement(parent, "protein");
  std::string str = std::to_string(prot_id);
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  int count = cluster_ids.size();
  str = std::to_string(count);
  xml_doc->addElement(prot_element, "compatible_proteoform_number", str.c_str());
  for (size_t i = 0; i < cluster_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, cluster_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(),
              Prsm::cmpEValueIncProtInc);
    geneXmlForProteoform(xml_doc, prot_element, select_prsm_ptrs,
                         mng_ptr, detail, add_ms);
  }
  return prot_element;
}
XmlDOMElement geneXmlForPrsmList(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                           const PrsmPtrVec &prsm_ptrs,
                                           int prot_id,
                                           const std::vector<int> &cluster_ids,
                                           const PrsmViewMngPtr &mng_ptr,
                                           bool detail, bool add_ms) {
  XmlDOMElement prot_element = xml_doc->addElement(parent, "prsm");
  std::string str = std::to_string(prot_id);
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  int count = cluster_ids.size();
  str = std::to_string(count);
  xml_doc->addElement(prot_element, "compatible_proteoform_number", str.c_str());
  for (size_t i = 0; i < cluster_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, cluster_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(),
              Prsm::cmpEValueIncProtInc);
    geneXmlForProteoform(xml_doc, prot_element, select_prsm_ptrs,
                         mng_ptr, detail, add_ms);
  }
  return prot_element;
}
}  // namespace anno_xml_util

}  // namespace toppic
