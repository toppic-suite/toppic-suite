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


#ifndef TOPPIC_BASE_ACTIVATION_HPP_
#define TOPPIC_BASE_ACTIVATION_HPP_

#include <string>
#include <vector>

#include "base/ion_type.hpp"

namespace toppic {

class XmlDOMDocument;

class Activation {
 public:
  Activation(const std::string &name, IonTypePtr n_ion_type_ptr,
             IonTypePtr c_ion_type_ptr);

  explicit Activation(XmlDOMElement * element);

  std::string getName() {return name_;}

  double getNShift() {return n_ion_type_ptr_->getBYShift();}

  double getCShift() {return c_ion_type_ptr_->getBYShift();}

  IonTypePtr getNIonTypePtr() {return n_ion_type_ptr_;}

  IonTypePtr getCIonTypePtr() {return c_ion_type_ptr_;}

  void appendNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  static std::string getNameFromXml(XmlDOMElement * element);

  static std::string getXmlElementName() {return "activation";}

 private:
  std::string name_;
  // n terminal ion type
  IonTypePtr n_ion_type_ptr_;
  // c terminal ion type
  IonTypePtr c_ion_type_ptr_;
};


typedef std::shared_ptr<Activation> ActivationPtr;
typedef std::vector<ActivationPtr> ActivationPtrVec;

}  // namespace toppic

#endif
