#include "spec/peak_tolerance.hpp"

namespace prot {

PeakTolerance::PeakTolerance(double ppo, bool use_min_tolerance,
                             double min_tolerance) {
  ppo_ = ppo;
  use_min_tolerance_ = use_min_tolerance;
  min_tolerance_ = min_tolerance;
}

PeakTolerance::PeakTolerance(xercesc::DOMElement* element){
  ppo_ = getDoubleChildValue(element,"ppo",0);
  use_min_tolerance_ = getDoubleChildValue(element,"use_min_tolerance",0);
  min_tolerance_ = getDoubleChildValue(element,"min_tolerance",0);
}

double PeakTolerance::compStrictErrorTole(double mass) {
  double tolerance = mass * ppo_;
  if (use_min_tolerance_ && tolerance < min_tolerance_) {
    tolerance = min_tolerance_;
  }
  return tolerance;
}

void PeakTolerance::appendXml(XmlDOMDocument* xml_doc, 
                              xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("peak_tolerance");
  std::string str = convertToString(ppo_);
  xml_doc->addElement(element, "ppo", str.c_str());
  str = convertToString(use_min_tolerance_);
  xml_doc->addElement(element, "use_min_tolerance", str.c_str());
  str = convertToString(min_tolerance_);
  xml_doc->addElement(element, "min_tolerance", str.c_str());
  parent->appendChild(element);
}

} /* namespace prot */
