#include "base/logger.hpp"

#include "base/residue_base.hpp"
#include "base/file_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

ResiduePtrVec ResidueBase::residue_ptr_vec_;

void ResidueBase::initBase(const std::string &file_name) {
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
      residue_ptr_vec_.push_back(ResiduePtr(new Residue(acid_name, ptm_abbr_name)));
    }
  }
}

ResiduePtr ResidueBase::getResiduePtrByAcidPtm(AcidPtr acid_ptr, 
                                               PtmPtr ptm_ptr) {
  for (size_t i = 0; i < residue_ptr_vec_.size(); i++) {
    if (residue_ptr_vec_[i]->isSame(acid_ptr, ptm_ptr)) {
      return residue_ptr_vec_[i];
    }
  }
  return ResiduePtr(nullptr);
}

ResiduePtr ResidueFactory::addResidue(AcidPtr acid_ptr, 
                                      PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = getResiduePtrByAcidPtm(acid_ptr, ptm_ptr);
  if (residue_ptr.get() == nullptr) {
    ResiduePtr new_ptr(new Residue(acid_ptr, ptm_ptr));
    residue_ptr_vec_.push_back(new_ptr);
    return new_ptr;
  }
  else {
    return residue_ptr;
  }
}

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

}
