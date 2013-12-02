#include "spec/peak_tolerance.hpp"

namespace prot {

PeakTolerance::PeakTolerance(double ppo, bool use_min_tolerance,
                             double min_tolerance) {
  ppo_ = ppo;
  use_min_tolerance_ = use_min_tolerance;
  min_tolerance_ = min_tolerance;
}
	
double PeakTolerance::compStrictErrorTole(double mass) {
  double tolerance = mass * ppo_;
  if (use_min_tolerance_ && tolerance < min_tolerance_) {
    tolerance = min_tolerance_;
  }
  return tolerance;
}

xercesc::DOMElement* PeakTolerance::toXml(XmlDOMDocument* xml_doc) {
  xercesc::DOMElement* element = xml_doc->createElement("peak_tolerance");
  std::string str = convertToString(ppo_);
  xml_doc->addElement(element, "ppo", str.c_str());
  str = convertToString(use_min_tolerance_);
  xml_doc->addElement(element, "use_min_tolerance", str.c_str());
  str = convertToString(min_tolerance_);
  xml_doc->addElement(element, "min_tolerance", str.c_str());
  return element;
}

} /* namespace prot */
