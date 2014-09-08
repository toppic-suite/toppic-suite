#ifndef PROT_ANNO_CHANGE_HPP_
#define PROT_ANNO_CHANGE_HPP_

#include <memory>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoChange {
 public:
  AnnoChange(int left_bp_pos, int right_bp_pos, 
             double mass_shift, int color, int type);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  int left_bp_pos_;
  int right_bp_pos_;
  double mass_shift_;
  int color_;
  int type_;
};

typedef std::shared_ptr<AnnoChange> AnnoChangePtr;
typedef std::vector<AnnoChangePtr> AnnoChangePtrVec;


}
#endif

