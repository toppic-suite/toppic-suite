#include <cmath>

#include "extreme_value.hpp"

namespace prot {

ExtremeValue::ExtremeValue (double one_prot_prob, double test_num, 
                            double adjust_factor) {
  one_prot_prob_ = one_prot_prob;
  test_num_ = test_num;
  adjust_factor_ = adjust_factor;
  init();
}

ExtremeValuePtr getMaxEvaluePtr() {
  ExtremeValuePtr evalue_ptr = ExtremeValuePtr(new ExtremeValue(1.0, 1.0, 1.0));
  return evalue_ptr;
}


ExtremeValue::ExtremeValue (xercesc::DOMElement* element){
  one_prot_prob_ = getDoubleChildValue(element, "one_protein_probability", 0);
  test_num_ = getDoubleChildValue(element, "test_number", 0);
  adjust_factor_ = getDoubleChildValue(element, "adjust_factor", 0);
  init();
}
    
void ExtremeValue::init() {
  e_value_ = one_prot_prob_ * test_num_ * adjust_factor_;
  if (one_prot_prob_ > 1) {
    p_value_  = 1.0;
  }
  else {
    double n = test_num_ * adjust_factor_;
    // approximation of 1 - (1- one_prot_prob)^n
    p_value_ =  n * one_prot_prob_ 
        - (n * (n-1))/2 * one_prot_prob_ * one_prot_prob_
        + (n*(n-1) * (n-2))/6 * std::pow(one_prot_prob_,3);
    if (p_value_ > 1.0) {
      p_value_ = 1.0;
    }
  }
}
    
void ExtremeValue::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("extreme_value");
  std::string str = convertToString(one_prot_prob_);
  xml_doc->addElement(element, "one_protein_probability", str.c_str());
  str = convertToString(test_num_);
  xml_doc->addElement(element, "test_number", str.c_str());
  str = convertToString(adjust_factor_);
  xml_doc->addElement(element, "adjust_factor", str.c_str());
  str = convertToString(p_value_);
  xml_doc->addElement(element, "p_value", str.c_str());
  str = convertToString(e_value_);
  xml_doc->addElement(element, "e_value", str.c_str());
  parent->appendChild(element);
}

}
