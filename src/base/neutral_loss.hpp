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


#ifndef PROT_BASE_NEUTRAL_LOSS_HPP_
#define PROT_BASE_NEUTRAL_LOSS_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/xml_dom_document.hpp"

namespace prot {

class NeutralLoss {
 public:
  NeutralLoss(const std::string &name, double mass): name_(name), mass_(mass) { }

  explicit NeutralLoss(xercesc::DOMElement* element);

  std::string getName() {return name_;}

  double getMass() {return mass_;}

  static std::string getXmlElementName() {return "neutral_loss";}

 private:
  std::string name_;

  double mass_;
};

typedef std::shared_ptr<NeutralLoss> NeutralLossPtr;

typedef std::vector<NeutralLossPtr> NeutralLossPtrVec;

}  // namespace prot

#endif /* NEUTRALLOSS_HPP_ */
