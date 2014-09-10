#ifndef PROT_ANNO_UNEXPECTED_CHANGE_HPP_
#define PROT_ANNO_UNEXPECTED_CHANGE_HPP_

#include <memory>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoUnexpectedChange {
 public:
  AnnoUnexpectedChange(int left_pos, int right_pos, 
                       double mass_shift, int color, const std::string &type);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                 int decimal_point_num);

 private:
  int left_pos_;
  int right_pos_;
  double mass_shift_;
  int color_;
  std::string type_;
};

typedef std::shared_ptr<AnnoUnexpectedChange> AnnoUnexpectedChangePtr;
typedef std::vector<AnnoUnexpectedChangePtr> AnnoUnexpectedChangePtrVec;

}
#endif

