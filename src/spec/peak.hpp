#ifndef PROT_PEAK_HPP_
#define PROT_PEAK_HPP_

#include <memory>
#include "base/xml_dom_document.hpp"
#include "base/mass_constant.hpp"

namespace prot {

class Peak {
 public:
  Peak(double position, double intensity) {
    position_ = position;
    intensity_ = intensity;
  }

  double getIntensity() {return intensity_;}

  double getPosition() {return position_;}

  void setIntensity(double intensity) {
    intensity_ = intensity;
  }

  void setPosition(double position) {
    position_ = position;
  }

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  	  xercesc::DOMElement* element = xml_doc->createElement("peak");
  	  std::string str = convertToString(getPosition());
  	  xml_doc->addElement(element, "position", str.c_str());
  	  str = convertToString(getIntensity());
  	  xml_doc->addElement(element, "intensity", str.c_str());
  	  parent->appendChild(element);
  }

 private:
  double position_;
  double intensity_;
};

inline double compPeakMass(double mono_mz, int charge) {
  return mono_mz * charge - charge * MassConstant::getProtonMass();
}

inline double compMonoMz(double mono_mass, int charge) {
  return mono_mass / charge + MassConstant::getProtonMass();
}

}
#endif
