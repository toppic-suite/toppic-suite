//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


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
