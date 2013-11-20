#include "ptm_vec.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace proteomics {
      

const std::vector<PtmPtr> PtmVec::getInstance(std::vector<AcidPtr> &acid_ptrs,
                                              const char* file_name) {
  std::vector<PtmPtr> ptm_ptrs;
  ptm_ptrs.push_back(Ptm::getEmptyPtmPtr(acid_ptrs));
  proteomics::XmlDOMParser* parser = proteomics::getXmlDOMInstance();
  if (parser) {
    proteomics::XmlDOMDocument* doc = new proteomics::XmlDOMDocument(parser, file_name);
    if (doc) {
      int acid_num = doc->getChildCount("ptm_list", 0, "ptm");
      for (int i = 0; i < acid_num; i++) {
        xercesc::DOMElement* element = doc->getElement("ptm", i);
        std::string abbr_name = getChildValue(element, "abbreviation");
        std::vector<AcidPtr> valid_acid_ptrs;
        double mono_mass = getDoubleChildValue(element, "mono_mass");
        ptm_ptrs.push_back(PtmPtr(new Ptm(abbr_name, valid_acid_ptrs, mono_mass)));

      }
      delete doc;
    }
    delete parser;
  }
  return ptm_ptrs;
}

/**
 *   Returns a PTM based on the abbreviation name. Returns null if the
 *   abbreviation name does not exist.
 */
PtmPtr PtmVec::getPtmByAbbrName(std::vector<PtmPtr> &ptm_ptrs, 
                               const std::string &abbr_name) {
  for (unsigned int i = 0; i < ptm_ptrs.size(); i++) {
    std::string n = ptm_ptrs[i]->getAbbrName();
    if (n == abbr_name) {
      return ptm_ptrs[i];
    }
  }
  return PtmPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool PtmVec::containAbbrsName(std::vector<PtmPtr> &ptm_ptrs, 
                      const std::string &abbr_name) {
    if (getPtmByAbbrName(ptm_ptrs, abbr_name).get() == nullptr) {
        return false;
    }
    else {
        return true;
    }
}

PtmPtr PtmVec::getEmptyPtm(std::vector<PtmPtr> &ptm_ptrs) {
  for (unsigned int i = 0; i < ptm_ptrs.size(); i++) {
    if (ptm_ptrs[i]->isEmpty()) {
      return ptm_ptrs[i];
    }
  }
  throw "Empty ptm does not exist!";
}

std::vector<AcidPtr> PtmVec::getValidAcids(std::vector<PtmPtr> &ptm_ptrs) {
  std::vector<AcidPtr> acid_ptrs;
  for (unsigned int i = 0; i < ptm_ptrs.size(); i++) {
    std::vector<AcidPtr> cur_ptrs = ptm_ptrs[i]->getValidAcids();
    acid_ptrs.insert(acid_ptrs.end(), cur_ptrs.begin(), cur_ptrs.end());
  }
  return acid_ptrs;
}

	
};


