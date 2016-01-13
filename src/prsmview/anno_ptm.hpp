#ifndef PROT_ANNO_PTM_HPP_
#define PROT_ANNO_PTM_HPP_

#include "base/ptm.hpp"
#include "base/change_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoPtm;
typedef std::shared_ptr<AnnoPtm> AnnoPtmPtr;
typedef std::vector<AnnoPtmPtr> AnnoPtmPtrVec;

class AnnoPtm {
 public:

  AnnoPtm(PtmPtr ptm_ptr, ChangeTypePtr change_type_ptr);

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  ChangeTypePtr getChangeTypePtr() {return change_type_ptr_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  void addOccurence(int pos, const std::string &acid_letter);

  static AnnoPtmPtr findPtm(const AnnoPtmPtrVec &ptm_ptrs, PtmPtr ptm_ptr, 
                            ChangeTypePtr change_type_ptr);

 private:
  PtmPtr ptm_ptr_;
  ChangeTypePtr change_type_ptr_;
  std::vector<std::pair<int, std::string>> occurences_;
};

}
#endif

