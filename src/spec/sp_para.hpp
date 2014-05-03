/*
 * sp_para.hpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SP_PARA_HPP_
#define PROT_SP_PARA_HPP_

#include <memory>
#include "spec/peak_tolerance.hpp"
#include "spec/extend_sp_para.hpp"
#include "base/activation.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class SpPara {
 public:
  SpPara(int min_peak_num,double min_mass, 
         const PeakTolerancePtr &peak_tolerance,
         const ExtendSpParaPtr &extend_sp_para,
         const ActivationPtr &activation);

  SpPara(xercesc::DOMElement* element);

  PeakTolerancePtr getPeakTolerancePtr(){return peak_tolerance_ptr_;}

  void setPeakTolerance(PeakTolerancePtr peak_tolerance_ptr){
    peak_tolerance_ptr_ = peak_tolerance_ptr;}

  ExtendSpParaPtr getExtendSpPara(){return extend_sp_para_;}

  void setExtendSpPara(const ExtendSpParaPtr &extend_sp_para){
    extend_sp_para_ = extend_sp_para;}

  ActivationPtr getActivation(){return activation_;}

  void setActivation(const ActivationPtr &activation){
    activation_ = activation;}

  int getMinPeakNum(){return min_peak_num_;}

  void setMinPeakNum(int min_peak_num){min_peak_num_=min_peak_num;}

  double getMinMass(){return min_mass_;}

  void setMinMass(double min_mass){min_mass_=min_mass;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

 private:
  int min_peak_num_;
  double min_mass_;
  PeakTolerancePtr peak_tolerance_ptr_;
  ExtendSpParaPtr extend_sp_para_;
  ActivationPtr activation_;
};

typedef std::shared_ptr<SpPara> SpParaPtr;

} /* namespace prot */

#endif /* SP_PARA_HPP_ */
