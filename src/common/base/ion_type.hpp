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

#ifndef TOPPIC_COMMON_BASE_ION_TYPE_HPP_
#define TOPPIC_COMMON_BASE_ION_TYPE_HPP_

#include <string>
#include <vector>
#include <memory>

#include "common/xml/xml_dom_element.hpp"

namespace toppic {

class XmlDOMDocument;

class IonType {
 public:
  IonType(const std::string &name, bool n_term, double shift);

  explicit IonType(XmlDOMElement* element);

  std::string getName() {return name_;}

  bool isNTerm() {return n_term_;}

  double getShift() {return shift_;}

  double getBYShift() {return b_y_shift_;}

  void appendNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getXmlElementName() {return "ion_type";}

 private:
  // ion name
  std::string name_;
  // A B C are n-terminal ions and X Y Z are non-n-terminal ions
  bool n_term_;
  /**
   * Shift stands for the shift of the ion compared to residue mass. For example, the
   * shift for B ion is 0, and the shift for Y ion is 18 (chrg 0);
   */
  double shift_;

  // shifts compared to b or y-ions
  double b_y_shift_;
};

typedef std::shared_ptr<IonType> IonTypePtr;
typedef std::vector<IonTypePtr> IonTypePtrVec;

}  // namespace toppic

#endif
