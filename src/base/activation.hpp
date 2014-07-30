/*
 * Activation.hpp
 *
 *  Created on: Nov 25, 2013
 *  Author: Xiaowen Liu
 */

#ifndef PROT_ACTIVATION_HPP_
#define PROT_ACTIVATION_HPP_

#include "xercesc/dom/DOM.hpp"
#include "base/ion_type.hpp"

namespace prot {

class Activation {
 public:
  Activation(const std::string &name, IonTypePtr n_ion_type_ptr, 
             IonTypePtr c_ion_type_ptr);

  Activation(xercesc::DOMElement * element);

  std::string getName() {return name_;}

  double getNShift() {return n_ion_type_ptr_->getShift()-getBIonShift();}

  double getCShift() {return c_ion_type_ptr_->getShift()-getYIonShift();}

  IonTypePtr getNIonTypePtr() {return n_ion_type_ptr_;}

  IonTypePtr getCIonTypePtr() {return c_ion_type_ptr_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  static double getBIonShift() {return 0;}
  static double getYIonShift() {return 18.0106;}

  std::string name_;
  // n terminal ion type 
  IonTypePtr n_ion_type_ptr_;
  // c terminal ion type
  IonTypePtr c_ion_type_ptr_;
};


typedef std::shared_ptr<Activation> ActivationPtr;
typedef std::vector<ActivationPtr> ActivationPtrVec;

/* activation factory */
class ActivationFactory {
 private:
  static ActivationPtrVec activation_ptr_vec_;

 public:
  static void initFactory(const std::string &file_name);

  static const ActivationPtrVec& getBaseActivationPtrVec() {return activation_ptr_vec_;}

  static ActivationPtr getBaseActivationPtrByName(const std::string &name);
};

} /* namespace prot */

#endif /* ACTIVATION_HPP_ */
