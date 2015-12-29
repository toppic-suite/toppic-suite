#ifndef PROT_ANNO_CHANGE_HPP_
#define PROT_ANNO_CHANGE_HPP_

#include "base/change.hpp"
#include "base/change_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoChange: public Change {
 public:

  void addOccurence(int pos, const std::string &acid_letter);

  void appendExpectedXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  void appendUnexpectedXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent, int precison);

 private:
  int color_;
  std::vector<std::pair<int, std::string>> occurences_;
};

typedef std::shared_ptr<AnnoChange> AnnoChangePtr;
typedef std::vector<AnnoChangePtr> AnnoChangePtrVec;

AnnoChangePtr findExpectedChange(const AnnoChangePtrVec &expected_change_ptrs, 
                                 ChangeTypePtr change_type_ptr, ModPtr mod_ptr);

}
#endif

