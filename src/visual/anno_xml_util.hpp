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

#ifndef TOPPIC_VISUAL_ANNO_XML_UTIL_HPP_
#define TOPPIC_VISUAL_ANNO_XML_UTIL_HPP_

#include <map>
#include <string>
#include <vector>

#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_writer.hpp"
#include "prsm/prsm.hpp"
#include "visual/prsm_view_mng.hpp"

namespace toppic {

namespace anno_xml_util {

// proteoform
XmlDOMElement geneXmlForProteoform(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                          const PrsmPtrVec &prsm_ptrs,
                                          const PrsmViewMngPtr &mng_ptr,
                                          bool detail = true,
                                          bool add_ms = true);
// single protein
void writeProteinToXml(const XmlWriterPtr &xml_writer,
                       const PrsmPtrVec &prsm_ptrs,
                       int prot_id,
                       const std::vector<int> &species_ids,
                       const PrsmViewMngPtr &mng_ptr, 
                       bool detail = true, 
                       bool add_ms = true);

// protein list
XmlDOMElement geneXmlForProteinList(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                           const PrsmPtrVec &prsm_ptrs,
                                           int prot_id,
                                           const std::vector<int> &cluster_ids,
                                           const PrsmViewMngPtr &mng_ptr,
                                           bool detail = true,
                                           bool add_ms = true);
// prsm list
XmlDOMElement geneXmlForPrsmList(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                           const PrsmPtrVec &prsm_ptrs,
                                           int prot_id,
                                           const std::vector<int> &cluster_ids,
                                           const PrsmViewMngPtr &mng_ptr, 
                                           bool detail = true, 
                                           bool add_ms = true);                                           
}  // namespace anno_xml_util

}  // namespace toppic

#endif /* TOPPIC_VISUAL_ANNO_XML_UTIL_HPP_ */
