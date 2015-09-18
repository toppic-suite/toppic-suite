#ifndef PROT_CHANGE_HPP_
#define PROT_CHANGE_HPP_

#include "base/ptm.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

#define INPUT_CHANGE      0
#define FIXED_CHANGE      1
#define PROTEIN_VARIABLE_CHANGE   2
#define VARIABLE_CHANGE   3
#define UNEXPECTED_CHANGE 4

class Change;
typedef std::shared_ptr<Change> ChangePtr;

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, const PtmPtr & ptm_ptr);

  Change(xercesc::DOMElement* change_element);

  Change(ChangePtr ori_change, int shift);
  
  int getLeftBpPos() {return left_bp_pos_;}

  void setLeftBpPos(int i) {left_bp_pos_ = i;}

  int getRightBpPos() {return right_bp_pos_;}

  void setRightBpPos(int i) { right_bp_pos_ = i;}

  void setSplit(int i) {split_pos_ = i;}

  int getChangeType() {return change_type_;}

  void setChangeType(int i) { change_type_ = i;}

  double getMassShift() {return mass_shift_;}

  void setMassShift(double m) { mass_shift_ = m;}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  void setPtmPtr(const PtmPtr& ptm) { ptm_ptr_ = ptm;}

  void setScr(std::vector<double> s) { scr_ = s;}

  std::vector<double> getScr() { return scr_;}

  void setConf(double c) { conf_ = c;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  // left and right positions are based on break point positions 
  int left_bp_pos_;
  int right_bp_pos_;
  int split_pos_ = 0;
  int change_type_;
  double mass_shift_;
  PtmPtr ptm_ptr_;
  std::vector<double> scr_;
  double conf_;
};

typedef std::vector<ChangePtr> ChangePtrVec;

inline bool compareChangeTypeUpPosUp(const ChangePtr &a, const ChangePtr &b) {
  if (a->getChangeType() < b->getChangeType()) {
    return true;
  }
  else if (a->getChangeType() > b->getChangeType()) {
    return false;
  }
  else {
    if (a->getLeftBpPos() < b->getLeftBpPos()) {
      return true;
    }
    else if (a->getLeftBpPos() > b->getLeftBpPos()) {
      return false;
    }
    else {
      return a->getRightBpPos() < b->getRightBpPos();
    }
  }
}

inline bool compareChangeUp(ChangePtr c1, ChangePtr c2) {
  if (c1->getLeftBpPos() < c2->getLeftBpPos()) {
    return true;
  }
  else if (c1->getLeftBpPos() > c2->getLeftBpPos()) {
    return false;
  }
  else {
    return c1->getRightBpPos() > c2->getRightBpPos();
  }
}

}
#endif

