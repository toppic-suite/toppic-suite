//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_COMMON_BASE_SUPPORT_PEAK_TYPE_HPP_
#define TOPPIC_COMMON_BASE_SUPPORT_PEAK_TYPE_HPP_

#include <string>
#include <memory>
#include <vector>

#include "common/xml/xml_dom_element.hpp"

namespace toppic {

class SupportPeakType {
 public:
  SupportPeakType(int id, const std::string &name);

  explicit SupportPeakType(XmlDOMElement* element);

  int getId() {return id_;}

  std::string getName() {return name_;}

  static std::string getXmlElementName() {return "support_peak_type";}

 private:
  int id_;

  std::string name_;
};

typedef std::shared_ptr<SupportPeakType> SPTypePtr;
typedef std::vector<SPTypePtr> SPTypePtrVec;

}  // namespace toppic

#endif /* SUPPORT_PEAK_TYPE_HPP_ */
