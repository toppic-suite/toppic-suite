#include "spec/deconv_peak.hpp"

namespace prot {

DeconvPeak::DeconvPeak (int id, double mono_mass, 
                        double intensity, int charge) 
    : Peak (mono_mass, intensity) {
      id_ = id;
      charge_ = charge;
      score_ = 1.0;
    }

DeconvPeak::DeconvPeak(xercesc::DOMElement* element):
    Peak (getDoubleChildValue(element,"position",0), getDoubleChildValue(element,"intensity",0)) {
      id_ = getIntChildValue(element,"id",0);
      charge_ = getIntChildValue(element,"charge",0);
      score_ = getDoubleChildValue(element,"score",0);
    }

void DeconvPeak::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("deconv_peak");
  std::string str = convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = convertToString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = convertToString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  parent->appendChild(element);
}

}

