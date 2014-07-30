
#include "base/neutral_loss.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

NeutralLossPtrVec NeutralLossFactory::neutral_loss_ptr_vec_;

NeutralLoss::NeutralLoss(const std::string &name, double mass){
  name_ = name;
  mass_ = mass;
}

void NeutralLossFactory::initFactory(const std::string &file_name){
  prot::XmlDOMParser* parser 
      = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int neutral_loss_num = getChildCount(parent, "neutral_loss");
    for (int i = 0; i < neutral_loss_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "neutral_loss", i);
      std::string name = getChildValue(element,"name", 0);
      double mass = getDoubleChildValue(element,"mass", 0);
      neutral_loss_ptr_vec_.push_back(
          NeutralLossPtr(new NeutralLoss(name, mass)));
    }
  }
}

NeutralLossPtr NeutralLossFactory::getBaseNeutralLossPtrByName(
    const std::string &name){
  for (size_t i = 0; i < neutral_loss_ptr_vec_.size(); i++) {
    std::string n = neutral_loss_ptr_vec_[i]->getName();
    if (n == name) {
      return neutral_loss_ptr_vec_[i];
    }
  }
  return NeutralLossPtr(nullptr);
}

} /* namespace prot */
