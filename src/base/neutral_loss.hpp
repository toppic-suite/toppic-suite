#ifndef PROT_BASE_NEUTRAL_LOSS_HPP_
#define PROT_BASE_NEUTRAL_LOSS_HPP_

#include <string>
#include <vector>
#include <memory>
#include "base/xml_dom_document.hpp"

namespace prot {

class NeutralLoss {

 public:
  NeutralLoss(const std::string &name, double mass);
  NeutralLoss(xercesc::DOMElement* element);
  std::string getName(){return name_;}
  double getMass(){return mass_;}

  static std::string getXmlElementName() {return "neutral_loss";}

 private:
  std::string name_;
  double mass_;
};

typedef std::shared_ptr<NeutralLoss> NeutralLossPtr;
typedef std::vector<NeutralLossPtr> NeutralLossPtrVec;

} /* namespace prot */

#endif /* NEUTRALLOSS_HPP_ */
