#ifndef PROT_BASE_CHANGE_HPP_
#define PROT_BASE_CHANGE_HPP_

#include "base/change_type.hpp"
#include "base/ptm.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Change;
typedef std::shared_ptr<Change> ChangePtr;

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos, 
         ChangeTypePtr change_type_ptr,
         double mass_shift, PtmPtr ptm_ptr);

  Change(xercesc::DOMElement* change_element);

  int getLeftBpPos() {return left_bp_pos_;}

  int getRightBpPos() {return right_bp_pos_;}

  ChangeTypePtr getChangeTypePtr() {return change_type_ptr_;}

  double getMassShift() {return mass_shift_;}

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "change";}

  static bool cmpPosIncrease(const ChangePtr &a, const ChangePtr &b);

 private:
  // left and right positions are based on break point positions 
  int left_bp_pos_;
  int right_bp_pos_;
  ChangeTypePtr change_type_ptr_;
  double mass_shift_;
  PtmPtr ptm_ptr_;
};

typedef std::vector<ChangePtr> ChangePtrVec;

}

#endif

