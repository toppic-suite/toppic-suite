#ifndef PROT_ANNO_SEGMENT_HPP_
#define PROT_ANNO_SEGMENT_HPP_

#include "base/change.hpp"
#include "base/change_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoSegment {
 public:

  void addOccurence(int pos, const std::string &acid_letter);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent, int precison);

 private:
  std::string segment_type_;
  int left_bp_pos_;
  int right_bp_pos_;
  double mass_shift_;
  int color_;
  PtmPtr ptm_ptr_;
  std::vector<std::pair<int, std::string>> occurences_;
};

typedef std::shared_ptr<AnnoSegment> AnnoSegmentPtr;
typedef std::vector<AnnoSegmentPtr> AnnoSegmentPtrVec;

}
#endif

