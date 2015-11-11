#ifndef PROT_SPEC_DECONV_PEAK_HPP_
#define PROT_SPEC_DECONV_PEAK_HPP_

#include <memory>
#include <vector>
#include "spec/peak.hpp"

namespace prot {

class DeconvPeak;
typedef std::shared_ptr<DeconvPeak> DeconvPeakPtr;

class DeconvPeak : public Peak {
 public:
  DeconvPeak (int id, double mono_mass, double intensity, int charge);

  DeconvPeak(xercesc::DOMElement* element);

  int getCharge() {return charge_;}

  int getId() {return id_;}

  double getMonoMass() {return getPosition();}

  double getMonoMz() {return compMonoMz(getMonoMass(), charge_);}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static bool cmpPosIncreasep(const DeconvPeakPtr &a, const DeconvPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

  static std::string getXmlElementName() {return "deconv_peak";}

 private:
  int id_;
  int charge_;
  double score_ = 1.0;
};

typedef std::vector<DeconvPeakPtr> DeconvPeakPtrVec;

}
#endif
