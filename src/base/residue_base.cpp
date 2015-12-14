#include "base/logger.hpp"
#include "base/acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/file_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ResiduePtrVec ResidueBase::residue_ptr_vec_;
ResiduePtr ResidueBase::empty_residue_ptr_;

void ResidueBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Residue::getXmlElementName();
    int residue_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG( "residue num " << residue_num);
    for (int i = 0; i < residue_num; i++) {
      xercesc::DOMElement* element
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ResiduePtr residue_ptr(new Residue(element));
      if (residue_ptr->getAcidPtr() == AcidBase::getEmptyAcidPtr() 
          && residue_ptr->getPtmPtr() == PtmBase::getEmptyPtmPtr()) {
        empty_residue_ptr_ = residue_ptr;
      }
      residue_ptr_vec_.push_back(residue_ptr);
    }
  }
}

// use residue in residue_base to remove duplications and reduce memory usage
// add the residue to base if it is a new one
ResiduePtr ResidueBase::getBaseResiduePtr(ResiduePtr residue_ptr) {
  for (size_t i = 0; i < residue_ptr_vec_.size(); i++) {
    if (residue_ptr_vec_[i]->isSame(residue_ptr)) {
      return residue_ptr_vec_[i];
    }
  }
  residue_ptr_vec_.push_back(residue_ptr);
  return residue_ptr;
}

ResiduePtrVec ResidueBase::getBaseNonePtmResiduePtrVec() {
  ResiduePtrVec result;
  for (size_t i = 0; i < residue_ptr_vec_.size(); i++) {
    if (residue_ptr_vec_[i]->getPtmPtr() == PtmBase::getEmptyPtmPtr()) {
      result.push_back(residue_ptr_vec_[i]);
    }
  }
  return result;
}

ResiduePtr ResidueBase::getBaseResiduePtr(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = ResiduePtr(new Residue(acid_ptr, ptm_ptr));
  return getBaseResiduePtr(residue_ptr);
}

ResiduePtr ResidueBase::getResiduePtrFromXml(xercesc::DOMElement * element) {
  ResiduePtr ptr(new Residue(element));
  return getBaseResiduePtr(ptr);
}

}
