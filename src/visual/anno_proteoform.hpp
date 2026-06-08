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

#ifndef TOPPIC_VISUAL_ANNO_PROTEOFORM_HPP_
#define TOPPIC_VISUAL_ANNO_PROTEOFORM_HPP_

#include "common/xml/xml_dom_document.hpp"
#include "prsm/prsm.hpp"
#include "visual/prsm_view_mng.hpp"

namespace toppic {

namespace anno_proteoform {

XmlDOMElement geneAnnoProteoform(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                        const PrsmPtr &prsm_ptr,
                                        const PrsmViewMngPtr &mng_ptr);
XmlDOMElement geneAnnoProteoformBrief(XmlDOMDocument* xml_doc, XmlDOMElement parent,
                                        const PrsmPtr &prsm_ptr,
                                        const PrsmViewMngPtr &mng_ptr);
}
}

#endif /* TOPPIC_VISUAL_ANNO_PROTEOFORM_HPP_ */
