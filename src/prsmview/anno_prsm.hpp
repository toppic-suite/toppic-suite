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


#ifndef PROT_ANNO_PRSM_HPP_
#define PROT_ANNO_PRSM_HPP_

#include "base/xml_dom_util.hpp"
#include "prsm/prsm.hpp"
#include "prsmview/prsm_view_mng.hpp"

namespace prot{

xercesc::DOMElement* geneAnnoPrsm(XmlDOMDocument* xml_doc, PrsmPtr prsm_ptr,
                                  PrsmViewMngPtr mng_ptr, bool detail = true,
                                  bool add_ms_peaks = true);

}

#endif /* PROT_ANNO_PRSM_HPP_ */
