#ifndef PROT_ANNO_PRSM_HPP_
#define PROT_ANNO_PRSM_HPP_

#include "base/xml_dom_util.hpp"
#include "prsm/prsm.hpp"
#include "prsmview/prsm_view_mng.hpp"

namespace prot{

xercesc::DOMElement* geneAnnoPrsm(XmlDOMDocument* xml_doc, PrsmPtr prsm_ptr, 
                                  PrsmViewMngPtr mng_ptr, bool detail = true);

}

#endif /* PROT_ANNO_PRSM_HPP_ */
