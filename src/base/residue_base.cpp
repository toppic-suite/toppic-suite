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

ResiduePtr ResidueBase::getResiduePtrFromXml(xercesc::DOMElement * element) {
  ResiduePtr ptr(new Residue(element));
  return getBaseResiduePtr(ptr);
}

/*
ResiduePtrVec ResidueBase::getResiduePtrVecInstance(const std::string &file_name) {
  ResiduePtrVec new_list;
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int residue_num = getChildCount(parent, "residue");
    LOG_DEBUG( "residue num " << residue_num);
    for (int i = 0; i < residue_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "residue", i);
      std::string acid_name = getChildValue(element, "acid", 0);
      std::string ptm_abbr_name = getChildValue(element, "ptm", 0);
      LOG_DEBUG("acid " << acid_name << " ptm " << ptm_abbr_name);
      AcidPtr acid_ptr = AcidFactory::getBaseAcidPtrByName(acid_name);
      if (acid_ptr.get() == nullptr) {
        LOG_ERROR( "acid " << acid_name  << " not found ");
        throw("acid not found");
      }
      PtmPtr ptm_ptr = PtmFactory::getBasePtmPtrByAbbrName(ptm_abbr_name);
      if (ptm_ptr.get() == nullptr) {
        LOG_ERROR( "ptm " << ptm_abbr_name  << " not found ");
        throw("ptm not found");
      }
      ResiduePtr residue_ptr = getBaseResiduePtrByAcidPtm(acid_ptr, ptm_ptr);
      if (residue_ptr.get() == nullptr) {
        residue_ptr = addBaseResidue(acid_ptr, ptm_ptr);
      }
      new_list.push_back(residue_ptr);
    }
  }
  return new_list;
}
*/

}
