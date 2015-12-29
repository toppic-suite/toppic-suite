#ifndef PROT_ANNO_PRSM_VIEW_HPP_
#define PROT_ANNO_PRSM_VIEW_HPP_

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
class AnnoPrsmView {
 public:
  std::vector<std::vector<std::string>> file_list_;
  xercesc::DOMElement* geneFileList(XmlDOMDocument* xml_doc);
};

typedef std::shared_ptr<AnnoPrsmView> AnnoPrsmViewPtr;

std::vector<std::vector<std::string>> readViewXmlFiles(const std::string &file_name);

xercesc::DOMElement* genePrsmView(XmlDOMDocument* xml_doc, PrsmPtr prsm_ptr, PrsmViewMngPtr mng_ptr);

}

#endif /* PROT_ANNO_VIEW_HPP_ */
