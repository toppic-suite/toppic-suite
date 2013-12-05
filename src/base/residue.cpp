#include <log4cxx/logger.h>

#include "base/residue.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

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
                                       std::string file_name) {
  ResiduePtrVec residue_list;
  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name.c_str());
    LOG4CXX_DEBUG(logger, "doc " << doc);
    xercesc::DOMElement* parent = doc->getDocumentElement();
    if (doc) {
      int residue_num = getChildCount(parent, "residue");
      LOG4CXX_DEBUG(logger, "residue num " << residue_num);
      for (int i = 0; i < residue_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "residue", i);
        std::string acid_name = getChildValue(element, "acid", 0);
        std::string ptm_abbr_name = getChildValue(element, "ptm", 0);
        LOG4CXX_DEBUG(logger, "acid vec " << acid_list.size() << " ptm vec " 
                      << ptm_list.size() << " acid " << acid_name << " ptm " << ptm_abbr_name);
        residue_list.push_back(ResiduePtr(
                new Residue(acid_list, ptm_list, acid_name, ptm_abbr_name)));
      }
      delete doc;
    }
  }
  return residue_list;
}

ResiduePtrVec getResiduePtrVecInstance(AcidPtrVec &acid_list, 
                                       PtmPtrVec &ptm_list,
                                       ResiduePtrVec &residue_list,
                                       std::string file_name) {
  ResiduePtrVec new_list;
  XmlDOMParser* parser = getXmlDOMInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name.c_str());
    LOG4CXX_DEBUG(logger, "doc " << doc);
    xercesc::DOMElement* parent = doc->getDocumentElement();
    if (doc) {
      int residue_num = getChildCount(parent, "residue");
      LOG4CXX_DEBUG(logger, "residue num " << residue_num);
      for (int i = 0; i < residue_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "residue", i);
        std::string acid_name = getChildValue(element, "acid", 0);
        std::string ptm_abbr_name = getChildValue(element, "ptm", 0);
        LOG4CXX_DEBUG(logger, "acid vec " << acid_list.size() << " ptm vec " << ptm_list.size() 
                      << " acid " << acid_name << " ptm " << ptm_abbr_name);
        AcidPtr acid_ptr = getAcidPtrByName(acid_list, acid_name);
        if (acid_ptr.get() == nullptr) {
          LOG4CXX_ERROR(logger, "acid " << acid_name  << " not found ");
          throw("acid not found");
        }
        PtmPtr ptm_ptr = getPtmPtrByAbbrName(ptm_list, ptm_abbr_name);
        if (ptm_ptr.get() == nullptr) {
          LOG4CXX_ERROR(logger, "ptm " << ptm_abbr_name  << " not found ");
          throw("ptm not found");
        }
        ResiduePtr residue_ptr = getResiduePtrByAcidPtm(residue_list, acid_ptr, ptm_ptr);
        if (residue_ptr.get() == nullptr) {
          residue_ptr = addResidue(residue_list, acid_ptr, ptm_ptr);
        }
        new_list.push_back(residue_ptr);
      }
      delete doc;
    }
  }
  return new_list;
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
