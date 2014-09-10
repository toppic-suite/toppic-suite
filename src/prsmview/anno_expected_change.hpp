#ifndef PROT_EXPECTED_CHANGE_HPP_
#define PROT_EXPECTED_CHANGE_HPP_

#include "base/ptm.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoExpectedChange {
 public:
  AnnoExpectedChange(int change_type, PtmPtr ptm_ptr);

  void addOccurence(int pos, const std::string &acid_letter);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  int getChangeType() {return change_type_;}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

 private:
  int change_type_;
  PtmPtr ptm_ptr_;
  std::vector<std::pair<int, std::string>> occurences_;
};

typedef std::shared_ptr<AnnoExpectedChange> AnnoExpectedChangePtr;
typedef std::vector<AnnoExpectedChangePtr> AnnoExpectedChangePtrVec;

AnnoExpectedChangePtr findExpectedChange(const AnnoExpectedChangePtrVec &expected_change_ptrs, 
                                         int change_type, PtmPtr ptm_ptr);

}
#endif

