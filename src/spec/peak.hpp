#ifndef PROT_SPEC_PEAK_HPP_
#define PROT_SPEC_PEAK_HPP_

#include <memory>

#include "base/mass_constant.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Peak {
 public:
  Peak(double position, double intensity);

  double getIntensity() {return intensity_;}

  double getPosition() {return position_;}

  void setIntensity(double intensity) {intensity_ = intensity;}

  void setPosition(double position) {position_ = position;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "peak";}

  static double compPeakMass(double mono_mz, int charge) {
    return mono_mz * charge - charge * MassConstant::getProtonMass();
  }

  static double compMonoMz(double mono_mass, int charge) {
    return mono_mass / charge + MassConstant::getProtonMass();
  }

 private:
  double position_;
  double intensity_;
};

typedef std::shared_ptr<Peak> PeakPtr;


}
#endif
