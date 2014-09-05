#ifndef ANNO_RESIDUE_HPP_
#define ANNO_RESIDUE_HPP_

#include "base/residue.hpp"


namespace prot {

#define ANNO_RESIDUE_TYPE_NORMAL "normal"
#define ANNO_RESIDUE_TYPE_N_TRUNCATION "n_truncation"
#define ANNO_RESIDUE_TYPE_C_TRUNCATION "c_truncation"
#define ANNO_RESIDUE_TYPE_KNOWN_CHANGE "known_change"

class AnnoResidue : public Residue {
 public:
  AnnoResidue(ResiduePtr residue_ptr, int display_pos);

  void setType(const std::string &type) {type_ = type;}

  void setUnexpectedChange(bool u) {is_unexpected_change_ = u;}

  void setUnexpectedChangeColor(int color) {unexpected_change_color_ = color;}

  void appendViewXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  //residues display pos
  int display_pos_ = 0;
  //residues type
  std::string type_ = ANNO_RESIDUE_TYPE_NORMAL;
  //is expected
  bool is_unexpected_change_ = false;
  //unexpected change color
  int unexpected_change_color_ = 0;
};

typedef std::shared_ptr<AnnoResidue> AnnoResiduePtr;
typedef std::vector<AnnoResiduePtr> AnnoResiduePtrVec;

}

#endif /* ANNO_RESIDUE_HPP_ */
