#ifndef PROT_PRSM_N_TERM_SHIFT_HPP_
#define PROT_PRSM_N_TERM_SHIFT_HPP_

#include <memory>
#include <vector>

#include <xercesc/dom/DOM.hpp>

#include "base/xml_dom_document.hpp"

namespace prot {

#define N_TERM_TRUNCATION 0
#define C_TERM_TRUNCATION 1

class NTermShift {
 public:
  NTermShift(double shift, int type);

  NTermShift(xercesc::DOMElement* element);

  double getShift() {return shift_;}
  int getType() {return type_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "n_term_shift";}

 private:
  double shift_;
  int type_;
};

typedef std::shared_ptr<NTermShift> NTermShiftPtr;
typedef std::vector<NTermShiftPtr> NTermShiftPtrVec;

} /* namespace prot */

#endif 
