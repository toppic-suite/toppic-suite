#ifndef PROT_SP_PARA_HPP_
#define PROT_SP_PARA_HPP_

#include <memory>
#include "spec/peak_tolerance.hpp"
#include "base/activation.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class SpPara {
 public:
  SpPara(int min_peak_num, double min_mass_, double min_extend_mass, 
         const std::vector<double> &ext_offsets,
         PeakTolerancePtr peak_tolerance_ptr,
         ActivationPtr activation_ptr);

  SpPara(xercesc::DOMElement* element);

  double getMinMass() {return min_mass_;}

  double getExtendMinMass() {return extend_min_mass_;}

  const std::vector<double>& getExtendOffsets() {return ext_offsets_;}

  PeakTolerancePtr getPeakTolerancePtr() {return peak_tolerance_ptr_;}

  void setPeakTolerancePtr(PeakTolerancePtr peak_tolerance_ptr){
    peak_tolerance_ptr_ = peak_tolerance_ptr;}

  ActivationPtr getActivationPtr() {return activation_ptr_;}

  void setActivationPtr(ActivationPtr activation_ptr) {
    activation_ptr_ = activation_ptr;}

  int getMinPeakNum() {return min_peak_num_;}

  void setMinPeakNum(int min_peak_num) {min_peak_num_=min_peak_num;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  int min_peak_num_;

  // if the mass if smaller than min_mass, the mass is removed. 
  double min_mass_;
  // if the mass is smaller than extend_min_mass, the peak is not extended 
  double extend_min_mass_;
  std::vector<double> ext_offsets_;

  PeakTolerancePtr peak_tolerance_ptr_;
  ActivationPtr activation_ptr_;
};

typedef std::shared_ptr<SpPara> SpParaPtr;

} /* namespace prot */

#endif /* SP_PARA_HPP_ */
