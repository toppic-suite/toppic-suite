#ifndef PROT_BASE_CHANGE_HPP_
#define PROT_BASE_CHANGE_HPP_

#include "base/change_type.hpp"
#include "base/mod.hpp"
#include "base/xml_dom_document.hpp"
#include "local/local_anno.hpp"

namespace prot {

class Change;
typedef std::shared_ptr<Change> ChangePtr;

class Change {
 public:
  Change(int left_bp_pos, int right_bp_pos, 
         ChangeTypePtr change_type_ptr,
         double mass_shift, ModPtr mod_ptr);

  Change(xercesc::DOMElement* change_element);

  int getLeftBpPos() {return left_bp_pos_;}

  void setLeftBpPos(int p) {left_bp_pos_ = p;}

  int getRightBpPos() {return right_bp_pos_;}

  void setRightBpPos(int p) {right_bp_pos_ = p;}

  ChangeTypePtr getChangeTypePtr() {return change_type_ptr_;}

  double getMassShift() {return mass_shift_;}

  void setMassShift(double m) {mass_shift_ = m;}

  ModPtr getModPtr() {return mod_ptr_;}

  LocalAnnoPtr getLocalAnno() {return local_anno_ptr_;}

  void setLocalAnno(LocalAnnoPtr p);

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "change";}

  static bool cmpPosInc(const ChangePtr &a, const ChangePtr &b);

  static ChangePtr geneChangePtr(ChangePtr ori_change_ptr, int start_pos);

 protected:
  // left and right positions are based on break point positions 
  int left_bp_pos_;
  int right_bp_pos_;
  ChangeTypePtr change_type_ptr_;
  double mass_shift_;
  ModPtr mod_ptr_;
  LocalAnnoPtr local_anno_ptr_;
};

typedef std::vector<ChangePtr> ChangePtrVec;

}

#endif

