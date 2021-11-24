//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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


#ifndef TOPPIC_MS_FEATURE_SINGLE_CHARGE_FEATURE_HPP_
#define TOPPIC_MS_FEATURE_SINGLE_CHARGE_FEATURE_HPP_

#include <memory>
#include <vector>

#include "common/xml/xml_dom_element.hpp"
#include "common/xml/xml_dom_document.hpp"

namespace toppic {

class SingleChargeFeature {
 public:

  SingleChargeFeature(int charge,
                      double time_begin, double time_end,
                      int scan_begin, int scan_end,
                      double intensity, int env_num);

  SingleChargeFeature(XmlDOMElement* element);

  int getCharge() {return charge_;}

  double getIntensity() {return intensity_;}

  double getTimeBegin() {return time_begin_;}

  double getTimeEnd() {return time_end_;}

  double getTimeMiddle() {return (time_begin_ + time_end_)/2;}

  int getScanBegin() {return scan_begin_;}

  int getScanEnd() {return scan_end_;}

  int getEnvNum() {return env_num_;}

  static std::string getXmlElementName() {return "single_charge_feature";}

  void appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement *parent);

 protected:
  int charge_;
  double time_begin_;
  double time_end_;
  int scan_begin_;
  int scan_end_;
  double intensity_;
  int env_num_ = 0;
};

typedef std::shared_ptr<SingleChargeFeature> SingleChargeFeaturePtr;
typedef std::vector<SingleChargeFeaturePtr> SingleChargeFeaturePtrVec;

}
#endif
