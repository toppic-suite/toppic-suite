//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_SPEC_PEAK_TOLERANCE_HPP_
#define PROT_SPEC_PEAK_TOLERANCE_HPP_

#include <memory>
#include <string>

#include "base/xml_dom_document.hpp"

namespace toppic {

class PeakTolerance {
 public:
  PeakTolerance(double ppo, bool use_min_tolerance,
                double min_tolerance):
      ppo_(ppo), 
      use_min_tolerance_(use_min_tolerance),
      min_tolerance_(min_tolerance) {}

  explicit PeakTolerance(xercesc::DOMElement* element);

  double compStrictErrorTole(double mass);

  // consider zero ptm relaxed error
  double compRelaxErrorTole(double m1, double m2) {
    return compStrictErrorTole(m1 + m2);
  }

  double getPpo() {return ppo_;}

  bool isUseMinTolerance() {return use_min_tolerance_;}

  double getMinTolerance() {return min_tolerance_;}

  void setPpo(double ppo) {ppo_ = ppo;}

  void setUseMinTolerance(bool use_min_tolerance) {
    use_min_tolerance_ = use_min_tolerance;
  }

  void setMinTolerance(double min_tolerance) {
    min_tolerance_ = min_tolerance;
  }

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "peak_tolerance";}

 private:
  double ppo_;
  /* whether or not use minimum tolerance */
  bool use_min_tolerance_;
  double min_tolerance_;
};

typedef std::shared_ptr<PeakTolerance> PeakTolerancePtr;

} /* namespace toppic */

#endif    
