/*
 * Activation.hpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#ifndef PROT_ACTIVATION_HPP_
#define PROT_ACTIVATION_HPP_

#include "xercesc/dom/DOM.hpp"
#include "base/ion_type.hpp"

namespace prot {

class Activation {
 public:
  Activation(std::string name, IonTypePtr n_ion_type_ptr, 
             IonTypePtr c_ion_type_ptr);

  Activation(IonTypePtrVec ion_type_list, xercesc::DOMElement * element);

  std::string getName(){return name_;}

  double getNShit(){return n_ion_type_->getShift()-getBIonShift();}

  double getCShit(){return c_ion_type_->getShift()-getYIonShift();}

  IonTypePtr getNIonType(){return n_ion_type_;}

  IonTypePtr getCIonType(){return c_ion_type_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  static double getBIonShift() {return 0;}
  static double getYIonShift() {return 18.0106;}

  std::string name_;
  IonTypePtr n_ion_type_;
  IonTypePtr c_ion_type_;
};


typedef std::shared_ptr<Activation> ActivationPtr;
typedef std::vector<ActivationPtr> ActivationPtrVec;

/* activation factory */
class ActivationFactory {
 private:
  static ActivationPtrVec activation_ptr_vec_;

 public:
  static void initFactory(const std::string file_name);

  static ActivationPtrVec& getActivationPtrVec() {return activation_ptr_vec_;}

  static ActivationPtr getActivationPtrByName(std::string name);
};

} /* namespace prot */

#endif /* ACTIVATION_HPP_ */
