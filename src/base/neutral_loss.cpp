/*
 * NeutralLoss.cpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#include "base/neutral_loss.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

NeutralLoss::NeutralLoss(std::string name,double mass){
	name_ = name;
	mass_ = mass;
}

NeutralLossPtrVec getNeutralLossPtrVecInstance(const char* file_name){
	NeutralLossPtrVec neutralLossPtrVec;
	prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMInstance();
  if (parser) {
    prot::XmlDOMDocument* doc = new prot::XmlDOMDocument(parser, file_name);
    if (doc) {
      xercesc::DOMElement* parent = doc->getDocumentElement();
      int neutral_loss_num = getChildCount(parent, "neutral_loss");
      for (int i = 0; i < neutral_loss_num; i++) {
        xercesc::DOMElement* element = getChildElement(parent, "neutral_loss", i);
        std::string name = getChildValue(element,"name", 0);
        double mass = getDoubleChildValue(element,"mass", 0);
        neutralLossPtrVec.push_back(NeutralLossPtr(new NeutralLoss(name, mass)));
      }
      delete doc;
    }
  }
  return neutralLossPtrVec;
}
NeutralLossPtr getNeutralLossPtrByName(NeutralLossPtrVec &neutralLoss_ptr_vec, const std::string &name){
	for (unsigned int i = 0; i < neutralLoss_ptr_vec.size(); i++) {
	    std::string n = neutralLoss_ptr_vec[i]->getName();
	    if (n == name) {
	      return neutralLoss_ptr_vec[i];
	    }
	  }
	  return NeutralLossPtr(nullptr);
}
} /* namespace prot */
