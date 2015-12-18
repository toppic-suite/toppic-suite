// Author: Xiaowen Liu
// Created on: Nov 25, 2013

#ifndef PROT_BASE_ACTIVATION_HPP_
#define PROT_BASE_ACTIVATION_HPP_

#include "base/ion_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Activation {
 public:
  Activation(const std::string &name, IonTypePtr n_ion_type_ptr, 
             IonTypePtr c_ion_type_ptr);

  Activation(xercesc::DOMElement * element);

  std::string getName() {return name_;}

  double getNShift() {return n_ion_type_ptr_->getBYShift();}
  
  double getCShift() {return c_ion_type_ptr_->getBYShift();}

  IonTypePtr getNIonTypePtr() {return n_ion_type_ptr_;}

  IonTypePtr getCIonTypePtr() {return c_ion_type_ptr_;}

  void appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getNameFromXml(xercesc::DOMElement * element);

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

}

#endif 
