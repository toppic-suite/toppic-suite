#include "residue.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

#include "log4cxx/logger.h"

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Residue"));

Residue::Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  acid_ptr_ = acid_ptr;
  ptm_ptr_ = ptm_ptr;
  mass_ = acid_ptr->getMonoMass() + ptm_ptr->getMonoMass();
}

Residue::Residue(AcidPtrVec acid_list, PtmPtrVec ptm_list,
          std::string acid_name, std::string ptm_abbr_name) {
  acid_ptr_ = getAcidPtrByName(acid_list, acid_name);
  ptm_ptr_ = getPtmPtrByAbbrName(ptm_list, ptm_abbr_name);
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

ResiduePtr getResiduePtrByAcid(ResiduePtrVec &residue_ptrs,
                                  AcidPtr acid_ptr) {
  for (unsigned int i = 0; i < residue_ptrs.size(); i++) {
    if (residue_ptrs[i]->getAcidPtr().get() == acid_ptr.get()) {
      return residue_ptrs[i];
    }
  }
  return ResiduePtr(nullptr);
}

ResiduePtr getResiduePtrByAcidPtm(ResiduePtrVec &residue_list,
                                  AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  for (unsigned int i = 0; i < residue_list.size(); i++) {
    if (residue_list[i]->isSame(acid_ptr, ptm_ptr)) {
      return residue_list[i];
    }
  }
  return ResiduePtr(nullptr);
}

ResiduePtrVec getResiduePtrVecInstance(AcidPtrVec &acid_list, 
                                       PtmPtrVec &ptm_list,
                                       const char* file_name) {
  ResiduePtrVec residue_list;
  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
    LOG4CXX_DEBUG(logger, "doc " << doc);
    xercesc::DOMElement* parent = doc->getDocumentElement();
    if (doc) {
      int residue_num = getChildCount(parent, "residue");
      LOG4CXX_DEBUG(logger, "residue num " << residue_num);
      for (int i = 0; i < residue_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "residue", i);
        std::string acid_name = getChildValue(element, "acid", 0);
        std::string ptm_abby_name = getChildValue(element, "ptm", 0);
        LOG4CXX_DEBUG(logger, "acid vec " << acid_list.size() << " ptm vec " << ptm_list.size() << " acid " << acid_name << " ptm " << ptm_abby_name);
        residue_list.push_back(ResiduePtr(
                new Residue(acid_list, ptm_list, acid_name, ptm_abby_name)));
      }
      delete doc;
    }
  }
  return residue_list;
}


ResiduePtr addResidue(ResiduePtrVec &residue_list, AcidPtr acid_ptr,
                      PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = getResiduePtrByAcidPtm(residue_list, acid_ptr, ptm_ptr);
  if (residue_ptr.get() == nullptr) {
    ResiduePtr new_ptr(new Residue(acid_ptr, ptm_ptr));
    residue_list.push_back(new_ptr);
    return new_ptr;
  }
  else {
    return residue_ptr;
  }
}

ResiduePtrVec convertAcidToResidueSeq(ResiduePtrVec residue_list,
                                      AcidPtrVec acid_ptrs) {
  ResiduePtrVec result_seq;
  for (unsigned int i = 0; i < acid_ptrs.size(); i++) {
    ResiduePtr residue_ptr = getResiduePtrByAcid(residue_list, acid_ptrs[i]);
    result_seq.push_back(residue_ptr);
  }
  return result_seq;
}

}
