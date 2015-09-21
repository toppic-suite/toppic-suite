#ifndef PROT_ANNO_VIEW_HPP_
#define PROT_ANNO_VIEW_HPP_

#include <map>

#include "base/extreme_value.hpp"
#include "base/proteoform.hpp"
#include "spec/deconv_peak.hpp"
#include "spec/extend_peak.hpp"
#include "spec/sp_para.hpp"
#include "prsm/prsm.hpp"
#include "prsmview/anno_residue.hpp"
#include "prsmview/anno_cleavage.hpp"
#include "prsmview/prsm_view_mng.hpp"

namespace prot{
class AnnoView {
 public:
  std::vector<std::vector<std::string>> file_list_;
  xercesc::DOMElement* geneFileList(XmlDOMDocument* xml_doc);
};

typedef std::shared_ptr<AnnoView> AnnoViewPtr;

std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name);

xercesc::DOMElement* genePrsmView(XmlDOMDocument* xml_doc, PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr);

xercesc::DOMElement* geneProteinView(XmlDOMDocument* xml_doc, PrsmPtr prsm_ptr,
                                     PrsmViewMngPtr mng_ptr, double err);

xercesc::DOMElement* proteoformToXml(XmlDOMDocument* xml_doc, const PrsmPtrVec &prsm_ptrs, 
                                     PrsmViewMngPtr mng_ptr);

xercesc::DOMElement* proteinToXml(XmlDOMDocument* xml_doc,
                                  const PrsmPtrVec &prsm_ptrs,
                                  ProteoformPtr proteo_ptr,
                                  const std::vector<int> &species_ids,
                                  PrsmViewMngPtr mng_ptr);


xercesc::DOMElement* allProteinToXml(XmlDOMDocument* xml_doc,
                                     const PrsmPtrVec &prsm_ptrs,
                                     PrsmViewMngPtr mng_ptr);

xercesc::DOMElement* allModificationToXml(XmlDOMDocument* xml_doc,
        PrsmPtrVec & prsm_ptrs, const PrsmViewMngPtr & mng_ptr);

std::vector<xercesc::DOMElement*> modificationToXml(XmlDOMDocument* xml_doc,
        PrsmPtrVec & prsm_ptrs, const PrsmViewMngPtr & mng_ptr);

}

#endif /* PROT_ANNO_VIEW_HPP_ */
