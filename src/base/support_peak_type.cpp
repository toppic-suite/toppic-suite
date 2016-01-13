#include "base/support_peak_type.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

SupportPeakType::SupportPeakType(int id, const std::string &name):
    id_(id), 
    name_(name) {
    }

SupportPeakType::SupportPeakType(xercesc::DOMElement* element) {
  id_ = XmlDomUtil::getIntChildValue(element, "id", 0);
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
}

} /* namespace prot */
