#ifndef ANNO_RESIDUE_HPP_
#define ANNO_RESIDUE_HPP_

#include "base/residue.hpp"

namespace prot {

class AnnoResidue : public Residue {
 public:
  AnnoResidue(ResiduePtr residue_ptr);

  void setPos(int pos) {pos_ = pos; }

  void setDisplayPos(int display_pos) {display_pos_ = display_pos;}

  void setType(const std::string &type) {type_ = type;}

  void setShiftStyle(const std::string &shift_type) {shift_style_ = shift_type;}

  std::string getType() {return type_;}

  void setIsModifyed(bool is_modified) {is_modified_ = is_modified;}

  bool isModified() {return is_modified_;}

  void setShift(double shift) {shift_ = shift;}

  double getShift() {return shift_;}

  void setExpected(bool e) {is_expected_ = e;}

  void setDisplayBg(int bg){display_bg_ = bg;}

  void appendViewXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  //residues pos
  int pos_ = 0;
  //residues type
  std::string type_ = "normal";
  //display pos shown in match seq
  int display_pos_ = 0;
  //display background color in html file
  int display_bg_ = -1;
  //is modified
  bool is_modified_ = false;
  //is expected
  bool is_expected_ = false;
  //expected shift type
  std::string shift_style_;
  //modify mass
  double shift_ = 0;
};

typedef std::shared_ptr<AnnoResidue> AnnoResiduePtr;
typedef std::vector<AnnoResiduePtr> AnnoResiduePtrVec;

}

#endif /* ANNO_RESIDUE_HPP_ */
