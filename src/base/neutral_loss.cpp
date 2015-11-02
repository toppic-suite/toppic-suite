
#include "base/neutral_loss.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

NeutralLoss::NeutralLoss(const std::string &name, double mass): 
    name_(name),
    mass_(mass) {
    }

NeutralLoss::NeutralLoss(xercesc::DOMElement* element) {
  name_ = XmlDomUtil::getChildValue(element,"name", 0);
  mass_ = XmlDomUtil::getDoubleChildValue(element,"mass", 0);
}

} /* namespace prot */
