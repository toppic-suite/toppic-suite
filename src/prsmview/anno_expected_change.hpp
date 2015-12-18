#ifndef PROT_EXPECTED_CHANGE_HPP_
#define PROT_EXPECTED_CHANGE_HPP_

#include "base/ptm.hpp"
#include "base/change_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoExpectedChange {
 public:
  AnnoExpectedChange(ChangeTypePtr change_type_ptr, ModPtr mod_ptr);

  void addOccurence(int pos, const std::string &acid_letter);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  ChangeTypePtr getChangeTypePtr() {return change_type_ptr_;}

  ModPtr getModPtr() {return mod_ptr_;}

 private:
  ChangeTypePtr change_type_ptr_;
  ModPtr mod_ptr_;
  std::vector<std::pair<int, std::string>> occurences_;
};

typedef std::shared_ptr<AnnoExpectedChange> AnnoExpectedChangePtr;
typedef std::vector<AnnoExpectedChangePtr> AnnoExpectedChangePtrVec;

AnnoExpectedChangePtr findExpectedChange(const AnnoExpectedChangePtrVec &expected_change_ptrs, 
                                         ChangeTypePtr change_type_ptr, ModPtr mod_ptr);

}
#endif

