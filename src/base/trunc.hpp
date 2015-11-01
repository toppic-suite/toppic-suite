#ifndef PROT_BASE_TRUNC_HPP_
#define PROT_BASE_TRUNC_HPP_

#include <string>
#include "base/acid.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Trunc {
 public:
  Trunc(const std::string &name, int trunc_len, 
        const std::string &acid_str);

  Trunc(xercesc::DOMElement* element); 

  const std::string& getName() {return name_;}

  int getTruncLen() {return trunc_len_;}

  const AcidPtrVec& getAcidPtrVec() {return acid_ptr_vec_;}

  double getShift() {return shift_;}

  static std::string getNameFromXml(xercesc::DOMElement * element);

  static std::string getXmlElementName() {return "truncation";}

 private:
  std::string name_;
  int trunc_len_;
  AcidPtrVec acid_ptr_vec_;
  double shift_;
};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

}

#endif
