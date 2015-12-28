#ifndef PROT_ANNO_UNEXPECTED_CHANGE_HPP_
#define PROT_ANNO_UNEXPECTED_CHANGE_HPP_

#include <memory>
#include <vector>

#include "base/ptm.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoUnexpectedChange {
 public:
  AnnoUnexpectedChange(int left_pos, int right_pos, 
                       double mass_shift, int color, const std::string &type);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                 int decimal_point_num);

  std::string getChangeType() {return type_;}

  int getRightPos() {return right_pos_;}

  void setRightPos(int right_pos) {right_pos_ = right_pos;}

  //void addOccurence(int pos, const std::string &acid_letter);

  //ModPtr getModPtr() {return mod_ptr_;}

  //void setModPtr(ModPtr m) {mod_ptr_ = m;}

 private:
  int left_pos_;
  int right_pos_;
  double mass_shift_;
  int color_;
  std::string type_;
  //ModPtr mod_ptr_;
  //std::vector<std::pair<int, std::string> > occurences_;
};

typedef std::shared_ptr<AnnoUnexpectedChange> AnnoUnexpectedChangePtr;
typedef std::vector<AnnoUnexpectedChangePtr> AnnoUnexpectedChangePtrVec;

}
#endif

