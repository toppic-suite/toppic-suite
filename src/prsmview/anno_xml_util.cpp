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

#include "common/util/str_util.hpp"
#include "common/base/residue_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/proteoform_factory.hpp"
#include "spec/peak.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_util.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_prsm.hpp"
#include "prsmview/anno_xml_util.hpp"

namespace toppic {

namespace anno_xml_util {

xercesc::DOMElement* geneXmlForProteoform(XmlDOMDocument* xml_doc, 
                                          const PrsmPtrVec &prsm_ptrs,
                                          PrsmViewMngPtr mng_ptr, 
                                          bool detail, bool add_ms) {
  xercesc::DOMElement* proteoform_element = xml_doc->createElement("compatible_proteoform");
  std::string str = str_util::toString(prsm_ptrs[0]->getProteoformPtr()->getProtId());
  xml_doc->addElement(proteoform_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(proteoform_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(proteoform_element, "sequence_description", str.c_str());
  str = str_util::toString(prsm_ptrs[0]->getProteoformPtr()->getProteoClusterId());
  xml_doc->addElement(proteoform_element, "proteoform_id", str.c_str());
  int count = prsm_ptrs.size();
  str = str_util::toString(count);
  xml_doc->addElement(proteoform_element, "prsm_number", str.c_str());
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    proteoform_element->appendChild(
        anno_prsm::geneAnnoPrsm(xml_doc, prsm_ptrs[i], mng_ptr, detail, add_ms));
  }
  return proteoform_element;
}

void writeProteinToXml(XmlWriterPtr xml_writer,
                       const PrsmPtrVec &prsm_ptrs,
                       int prot_id,
                       const std::vector<int> &species_ids,
                       PrsmViewMngPtr mng_ptr,
                       bool detail, bool add_ms) {
  xml_writer->write_str("<protein>");
  xml_writer->write_str("<sequence_id>" + str_util::toString(prot_id) + "</sequence_id>");
  xml_writer->write_str("<sequence_name>" 
                        + prsm_ptrs[0]->getProteoformPtr()->getSeqName() + "</sequence_name>");
  xml_writer->write_str("<sequence_description>" 
                        + prsm_ptrs[0]->getProteoformPtr()->getSeqDesc() 
                        + "</sequence_description>");
  xml_writer->write_str("<compatible_proteoform_number>" 
                        + str_util::toString(species_ids.size()) 
                        + "</compatible_proteoform_number>");
  for (size_t i = 0; i < species_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, species_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(), Prsm::cmpEValueInc);
    xml_writer->write(geneXmlForProteoform(xml_writer->getDoc(), select_prsm_ptrs, mng_ptr, 
                                           detail, add_ms));
  }
  xml_writer->write_str("</protein>");
}

xercesc::DOMElement* geneXmlForProteinList(XmlDOMDocument* xml_doc,
                                           const PrsmPtrVec &prsm_ptrs,
                                           int prot_id,
                                           const std::vector<int> &cluster_ids,
                                           PrsmViewMngPtr mng_ptr,
                                           bool detail, bool add_ms) {
  xercesc::DOMElement* prot_element = xml_doc->createElement("protein");
  std::string str = str_util::toString(prot_id);
  xml_doc->addElement(prot_element, "sequence_id", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqName();
  xml_doc->addElement(prot_element, "sequence_name", str.c_str());
  str = prsm_ptrs[0]->getProteoformPtr()->getSeqDesc();
  xml_doc->addElement(prot_element, "sequence_description", str.c_str());
  int count = cluster_ids.size();
  str = str_util::toString(count);
  xml_doc->addElement(prot_element, "compatible_proteoform_number", str.c_str());
  for (size_t i = 0; i < cluster_ids.size(); i++) {
    PrsmPtrVec select_prsm_ptrs = prsm_util::selectClusterPrsms(prsm_ptrs, cluster_ids[i]);
    std::sort(select_prsm_ptrs.begin(), select_prsm_ptrs.end(), Prsm::cmpEValueInc);
    prot_element->appendChild(geneXmlForProteoform(xml_doc, select_prsm_ptrs, 
                                                   mng_ptr, detail, add_ms));
  }
  return prot_element;
}

}  // namespace anno_xml_util

}  // namespace toppic
