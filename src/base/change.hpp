#ifndef PROT_CHANGE_HPP_
#define PROT_CHANGE_HPP_

#include "base/ptm.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Change;

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos, int change_type,
         double mass_shift, const PtmPtr &ptm_ptr);
  Change(xercesc::DOMElement* change_element);

  int getLeftBpPos() {return left_bp_pos_;}

  int getRightBpPos() {return right_bp_pos_;}

  int getChangeType() {return change_type_;}

  double getMassShift() {return mass_shift_;}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  // left and right positions are based on break point positions 
  int left_bp_pos_;
  int right_bp_pos_;
  int change_type_;
  double mass_shift_;
  PtmPtr ptm_ptr_;
};

typedef std::shared_ptr<Change> ChangePtr;
typedef std::vector<ChangePtr> ChangePtrVec;

bool compareChangePosUp(ChangePtr c1, ChangePtr c2);


}
#endif

