#include <base/logger.hpp>

#include "base/residue.hpp"
#include "base/file_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

ResiduePtrVec ResidueFactory::residue_ptr_vec_;

Residue::Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  acid_ptr_ = acid_ptr;
  ptm_ptr_ = ptm_ptr;
  mass_ = acid_ptr->getMonoMass() + ptm_ptr->getMonoMass();
}

Residue::Residue(const std::string &acid_name, 
                 const std::string &ptm_abbr_name) {
  acid_ptr_ = AcidFactory::getBaseAcidPtrByName(acid_name);
  ptm_ptr_ = PtmFactory::getBasePtmPtrByAbbrName(ptm_abbr_name);
  mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
}

std::string Residue::toString(const std::string &delim_bgn, 
                              const std::string &delim_end) {
  if (ptm_ptr_->isEmpty()) {
    return acid_ptr_->getOneLetter();
  } else {
    return acid_ptr_->getOneLetter() + delim_bgn + ptm_ptr_->getAbbrName()
        + delim_end;
  }
}

void Residue::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
    xercesc::DOMElement* element = xml_doc->createElement("residue");
    std::string str = convertToString(mass_);
    xml_doc->addElement(element, "mass", str.c_str());
    acid_ptr_->appendxml(xml_doc,element);
    ptm_ptr_->appendxml(xml_doc,element);
    parent->appendChild(element);
}

ResiduePtr getResiduePtrByAcid(const ResiduePtrVec &residue_ptr_vec,
                               AcidPtr acid_ptr) {
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    if (residue_ptr_vec[i]->getAcidPtr() == acid_ptr) {
      return residue_ptr_vec[i];
    }
  }
  return ResiduePtr(nullptr);
}

int findResidue(const ResiduePtrVec &residue_list, ResiduePtr residue_ptr) {
  for (size_t i = 0; i < residue_list.size(); i++) {
    if (residue_list[i] == residue_ptr) {
      return i;
    }
  }
  return -1;
}


ResiduePtrVec convertAcidToResidueSeq(const ResiduePtrVec &residue_list,
                                      const AcidPtrVec &acid_ptrs) {
  ResiduePtrVec result_seq;
  for (size_t i = 0; i < acid_ptrs.size(); i++) {
    ResiduePtr residue_ptr = getResiduePtrByAcid(residue_list, acid_ptrs[i]);
    result_seq.push_back(residue_ptr);
  }
  return result_seq;
}

void ResidueFactory::initFactory(const std::string &file_name) {
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

ResiduePtr ResidueFactory::getBaseResiduePtrByAcidPtm(AcidPtr acid_ptr, 
                                                      PtmPtr ptm_ptr) {
  for (size_t i = 0; i < residue_ptr_vec_.size(); i++) {
    if (residue_ptr_vec_[i]->isSame(acid_ptr, ptm_ptr)) {
      return residue_ptr_vec_[i];
    }
  }
  return ResiduePtr(nullptr);
}

ResiduePtr ResidueFactory::addBaseResidue(AcidPtr acid_ptr, 
                                          PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = getBaseResiduePtrByAcidPtm(acid_ptr, ptm_ptr);
  if (residue_ptr.get() == nullptr) {
    ResiduePtr new_ptr(new Residue(acid_ptr, ptm_ptr));
    residue_ptr_vec_.push_back(new_ptr);
    return new_ptr;
  }
  else {
    return residue_ptr;
  }
}

ResiduePtrVec ResidueFactory::getResiduePtrVecInstance(const std::string &file_name) {
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
      LOG_DEBUG(" acid " << acid_name << " ptm " << ptm_abbr_name);
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

/** FixResidueFactory */

std::map<std::string,ResiduePtrVec> FixResidueFactory::fix_res_list_map_;

void FixResidueFactory::initFactory(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  std::string file_dir = directory(file_name);
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int residue_list_num = getChildCount(parent, "fix_mod_residue_list");
    LOG_DEBUG( "residue list num " << residue_list_num);
    for (int i = 0; i < residue_list_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "fix_mod_residue_list", i);
      std::string id = getChildValue(element, "id", 0);
      std::string fix_mod_residue_file_name = getChildValue(element, "fix_mod_residue_file_name", 0);
      fix_mod_residue_file_name = file_dir + FILE_SEPARATOR + fix_mod_residue_file_name;
      ResiduePtrVec residue_ptr_vec = 
          ResidueFactory::getResiduePtrVecInstance(fix_mod_residue_file_name);
      std::pair<std::string,ResiduePtrVec> pair(id, residue_ptr_vec);
      fix_res_list_map_.insert(pair);
    }
  }
}

}
