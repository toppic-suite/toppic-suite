#include "residue_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace prot {

ResiduePtr getResiduePtrByAcid(ResiduePtrVec residue_ptr_vec, 
                               AcidPtr acid_ptr) {
  for (unsigned int i = 0; i < residue_ptr_vec.size(); i++) {
    if (residue_ptr_vec[i]->getAcidPtr().get() == acid_ptr.get()) {
      return residue_ptr_vec[i];
    }
  }
  return ResiduePtr(nullptr);
}


ResiduePtr getResiduePtrByAcidPtm(ResiduePtrVec residue_ptr_vec, 
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

}

