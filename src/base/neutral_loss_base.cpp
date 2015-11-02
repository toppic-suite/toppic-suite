
#include "base/neutral_loss_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

NeutralLossPtrVec NeutralLossBase::neutral_loss_ptr_vec_;
NeutralLossPtr NeutralLossBase::neutral_loss_ptr_NONE_;


void NeutralLossBase::initBase(const std::string &file_name){
  prot::XmlDOMParser* parser 
      = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = NeutralLoss::getXmlElementName();
    int neutral_loss_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < neutral_loss_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      NeutralLossPtr neutral_loss_ptr(new NeutralLoss(element));
      neutral_loss_ptr_vec_.push_back(neutral_loss_ptr);
      if (neutral_loss_ptr->getName() == getName_NONE()) {
        neutral_loss_ptr_NONE_ = neutral_loss_ptr;
      }
    }
  }
}

NeutralLossPtr NeutralLossBase::getNeutralLossPtrByName(
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
