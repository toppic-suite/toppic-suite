#include "residue.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

Residue::Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  acid_ptr_ = acid_ptr;
  ptm_ptr_ = ptm_ptr;
  mass_ = acid_ptr->getMonoMass() + ptm_ptr->getMonoMass();
}

Residue::Residue(AcidPtrVec acid_ptr_vec, PtmPtrVec ptm_ptr_vec,
          std::string acid_name, std::string ptm_abbr_name) {
  acid_ptr_ = getAcidPtrByName(acid_ptr_vec, acid_name);
  ptm_ptr_ = getPtmPtrByAbbrName(ptm_ptr_vec, ptm_abbr_name);
  mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
}

std::string Residue::toString(std::string delim_bgn, std::string delim_end) {
  if (ptm_ptr_->isEmpty()) {
    return acid_ptr_->getOneLetter();
  } else {
    return acid_ptr_->getOneLetter() + delim_bgn + ptm_ptr_->getAbbrName()
        + delim_end;
  }
}

ResiduePtr getResiduePtrByAcidPtm(ResiduePtrVec &residue_ptr_vec,
                                  AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  for (unsigned int i = 0; i < residue_ptr_vec.size(); i++) {
    if (residue_ptr_vec[i]->isSame(acid_ptr, ptm_ptr)) {
      return residue_ptr_vec[i];
    }
  }
  return ResiduePtr(nullptr);
}

ResiduePtrVec getResiduePtrVecInstance(AcidPtrVec &acid_ptr_vec, 
                                       PtmPtrVec &ptm_ptr_vec,
                                       const char* file_name) {
  ResiduePtrVec residue_ptr_vec;
  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
    if (doc) {
      int residue_num = doc->getChildCount("residue_list", 0, "residue");
      for (int i = 0; i < residue_num; i++) {
        xercesc::DOMElement* element = doc->getElement("residue", i);
        std::string acid_name = getChildValue(element, "acid");
        std::string ptm_abby_name = getChildValue(element, "ptm");
        residue_ptr_vec.push_back(ResiduePtr(
                new Residue(acid_ptr_vec, ptm_ptr_vec, acid_name, ptm_abby_name)));
      }
      delete doc;
    }
    delete parser;
  }
  return residue_ptr_vec;
}


ResiduePtr addResidue(ResiduePtrVec &residue_ptr_vec, AcidPtr acid_ptr,
                      PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = getResiduePtrByAcidPtm(residue_ptr_vec, acid_ptr, ptm_ptr);
  if (residue_ptr.get() == nullptr) {
    ResiduePtr new_ptr(new Residue(acid_ptr, ptm_ptr));
    residue_ptr_vec.push_back(new_ptr);
    return new_ptr;
  }
  else {
    return residue_ptr;
  }
}

}
