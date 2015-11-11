#ifndef PROT_EXTEND_PEAK_HPP_
#define PROT_EXTEND_PEAK_HPP_

#include <vector>
#include <memory>
#include <algorithm>

#include "spec/deconv_peak.hpp"

namespace prot {

class ExtendPeak;
typedef std::shared_ptr<ExtendPeak> ExtendPeakPtr;

class ExtendPeak : public Peak{
 public:
  ExtendPeak();

  ExtendPeak(DeconvPeakPtr base_peak_ptr, double mono_mass, double score);

  DeconvPeakPtr getBasePeakPtr(){return base_peak_ptr_;}

  double getMonoMass(){return mono_mass_;}

  double getScore(){return score_;}

  double getOrigTolerance(){return orig_tolerance_;}

  double getReverseTolerance(){return reverse_tolerance_;}

  void setOrigTolerance(double orig_tolerance) {
    orig_tolerance_ = orig_tolerance;}

  void setReverseTolerance(double reverse_tolerance) {
    reverse_tolerance_ = reverse_tolerance;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static bool cmpPosIncrease(const ExtendPeakPtr &a, const ExtendPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

  static std::string getXmlElementName() {return "extend_peak";}

 private:
  DeconvPeakPtr base_peak_ptr_;
  double mono_mass_;
  double score_;
  double orig_tolerance_;
  double reverse_tolerance_;
};

typedef std::vector<ExtendPeakPtr> ExtendPeakPtrVec;



} /* namespace prot */

#endif /* EXTEND_PEAK_HPP_ */
