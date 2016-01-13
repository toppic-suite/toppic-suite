#ifndef PROT_BASE_SUPPORT_PEAK_TYPE_HPP_
#define PROT_BASE_SUPPORT_PEAK_TYPE_HPP_

#include <memory>
#include <string>
#include <vector>

#include "base/xml_dom_document.hpp"

namespace prot {

class SupportPeakType {
 public:
  SupportPeakType(int id, const std::string &name);

  SupportPeakType(xercesc::DOMElement* element);

  int getId(){return id_;}

  std::string getName(){return name_;}

  static std::string getXmlElementName() {return "support_peak_type";}

 private:
  int id_;
  std::string name_;
};

typedef std::shared_ptr<SupportPeakType> SPTypePtr;
typedef std::vector<SPTypePtr> SPTypePtrVec;

} /* namespace prot */

#endif /* SUPPORT_PEAK_TYPE_HPP_ */
