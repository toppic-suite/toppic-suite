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


#ifndef PROT_ANNO_VIEW_HPP_
#define PROT_ANNO_VIEW_HPP_

#include <map>
#include <string>
#include <vector>

#include "prsm/extreme_value.hpp"
#include "seq/proteoform.hpp"
#include "base/xml_writer.hpp"
#include "spec/deconv_peak.hpp"
#include "spec/extend_peak.hpp"
#include "spec/sp_para.hpp"
#include "prsm/prsm.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_cleavage.hpp"
#include "prsmview/prsm_view_mng.hpp"

namespace toppic {

class AnnoView {
 public:
  std::vector<std::vector<std::string>> file_list_;
  xercesc::DOMElement* geneFileList(XmlDOMDocument* xml_doc);
};

typedef std::shared_ptr<AnnoView> AnnoViewPtr;


std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name);

xercesc::DOMElement* proteoformToXml(XmlDOMDocument* xml_doc, const PrsmPtrVec &prsm_ptrs,
                                     PrsmViewMngPtr mng_ptr, bool detail = true, bool add_ms = true);

xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  int prot_id,
                                  const std::vector<int> &cluster_ids,
                                  PrsmViewMngPtr mng_ptr, bool detail = true, bool add_ms = true);

void writeProteinToXml(XmlWriterPtr xml_writer,
                       const PrsmPtrVec &prsm_ptrs,
                       int prot_id,
                       const std::vector<int> &species_ids,
                       PrsmViewMngPtr mng_ptr, bool detail = true, bool add_ms = true);

}  // namespace toppic

#endif /* PROT_ANNO_VIEW_HPP_ */
