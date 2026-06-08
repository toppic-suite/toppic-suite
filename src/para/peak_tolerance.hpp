//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_PARA_PEAK_TOLERANCE_HPP_
#define TOPPIC_PARA_PEAK_TOLERANCE_HPP_

#include <memory>
#include <string>
#include <vector>

#include "common/xml/xml_dom_element.hpp"

namespace toppic {

class XmlDOMDocument;

class PeakTolerance {
 public:
  explicit PeakTolerance(double ppo);

  explicit PeakTolerance(XmlDOMElement element);

  double compStrictErrorTole(double mass) const;

  // consider zero ptm relaxed error
  double compRelaxErrorTole(double m1, double m2) const;

  double getPpo() const {return ppo_;}

  int getIntPpm() const;

  bool isUseMinTolerance() const {return use_min_tolerance_;}

  double getMinTolerance() const {return min_tolerance_;}

  void setPpo(double ppo) {ppo_ = ppo;}

  void setUseMinTolerance(bool use_min_tolerance) {
    use_min_tolerance_ = use_min_tolerance;}

  void setMinTolerance(double min_tolerance) {
    min_tolerance_ = min_tolerance;}

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement parent) const;

  static std::string getXmlElementName() {return "peak_tolerance";}

 private:
  double ppo_;
  /* whether or not use minimum tolerance */
  bool use_min_tolerance_ = true;
  double min_tolerance_ = 0.01;
};

using PeakTolerancePtr = std::shared_ptr<PeakTolerance>;
using PeakTolerancePtrVec = std::vector<PeakTolerancePtr>;

}  // namespace toppic

#endif
