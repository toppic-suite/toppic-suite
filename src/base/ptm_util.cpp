
#include "ptm_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace proteomics {
      

PtmPtrVec getPtmPtrVecInstance(AcidPtrVec &acid_ptr_vec,
                               const char* file_name) {
  PtmPtrVec ptm_ptr_vec;
  ptm_ptr_vec.push_back(Ptm::getEmptyPtmPtr(acid_ptr_vec));
  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
    if (doc) {
      int ptm_num = doc->getChildCount("ptm_list", 0, "ptm");
      for (int i = 0; i < ptm_num; i++) {
        xercesc::DOMElement* element = doc->getElement("ptm", i);
        std::string abbr_name = getChildValue(element, "abbreviation");
        AcidPtrVec valid_acid_ptr_vec;
        double mono_mass = getDoubleChildValue(element, "mono_mass");
        ptm_ptr_vec.push_back(PtmPtr(new Ptm(abbr_name, valid_acid_ptr_vec, mono_mass)));

      }
      delete doc;
    }
    delete parser;
  }
  return ptm_ptr_vec;
}

/**
 *   Returns a PTM based on the abbreviation name. Returns null if the
 *   abbreviation name does not exist.
 */
PtmPtr getPtmPtrByAbbrName(PtmPtrVec &ptm_ptr_vec, 
                           const std::string &abbr_name) {
  for (unsigned int i = 0; i < ptm_ptr_vec.size(); i++) {
    std::string n = ptm_ptr_vec[i]->getAbbrName();
    if (n == abbr_name) {
      return ptm_ptr_vec[i];
    }
  }
  return PtmPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containAbbrsName(PtmPtrVec &ptm_ptr_vec, 
                      const std::string &abbr_name) {
    if (getPtmPtrByAbbrName(ptm_ptr_vec, abbr_name).get() == nullptr) {
        return false;
    }
    else {
        return true;
    }
}

PtmPtr findEmptyPtmPtr(PtmPtrVec &ptm_ptr_vec) {
  for (unsigned int i = 0; i < ptm_ptr_vec.size(); i++) {
    if (ptm_ptr_vec[i]->isEmpty()) {
      return ptm_ptr_vec[i];
    }
  }
  throw "Empty ptm does not exist!";
}

AcidPtrVec getValidAcidPtrVec(PtmPtrVec &ptm_ptr_vec) {
  AcidPtrVec acid_ptr_vec;
  for (unsigned int i = 0; i < ptm_ptr_vec.size(); i++) {
    AcidPtrVec cur_ptrs = ptm_ptr_vec[i]->getValidAcidPtrVec();
    acid_ptr_vec.insert(acid_ptr_vec.end(), cur_ptrs.begin(), cur_ptrs.end());
  }
  return acid_ptr_vec;
}

	
};


