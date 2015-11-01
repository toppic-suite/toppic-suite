#ifndef PROT_BASE_PROT_MOD_HPP_
#define PROT_BASE_PROT_MOD_HPP_

#include "base/ptm.hpp"
#include "base/trunc.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class ProtMod {
 public:
  ProtMod(const std::string &name, TruncPtr trunc_ptr,
          PtmPtr ptm_ptr);

  ProtMod(xercesc::DOMElement* element); 

  const std::string& getName() { return name_;};

  TruncPtr getTruncPtr() { return trunc_ptr_;}

  PtmPtr getPtmPtr() { return ptm_ptr_;}

  double getProtShift() { return prot_shift_;}

  double getPepShift() { return pep_shift_;}

  void appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "prot_mod";}

 private:
  std::string name_;
  TruncPtr trunc_ptr_;
  PtmPtr ptm_ptr_;
  double prot_shift_;
  double pep_shift_;
};

typedef std::shared_ptr<ProtMod> ProtModPtr;
typedef std::vector<ProtModPtr> ProtModPtrVec;

}
#endif
