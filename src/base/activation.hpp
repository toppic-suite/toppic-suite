/*
 * Activation.hpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#ifndef PROT_ACTIVATION_HPP_
#define PROT_ACTIVATION_HPP_

#include <xercesc/dom/DOM.hpp>
#include "ion_type.hpp"

namespace prot {

class Activation {
public:
	Activation(std::string name,IonTypePtr n_ion_type_ptr,IonTypePtr c_ion_type_ptr);
	Activation(IonTypePtrVec ion_type_list,xercesc::DOMElement * element);
	std::string getName(){return name_;}
	IonTypePtr getNIonType(){return n_ion_type_;}
	IonTypePtr getCIonType(){return c_ion_type_;}
private:
	std::string name_;
	IonTypePtr n_ion_type_;
	IonTypePtr c_ion_type_;
};

typedef std::shared_ptr<Activation> ActivationPtr;
typedef std::vector<ActivationPtr> ActivationPtrVec;

ActivationPtrVec getActivationPtrVectInst(const char* file_name);
ActivationPtr getActivationPtrByName(ActivationPtrVec activation_list,std::string name);

} /* namespace prot */

#endif /* ACTIVATION_HPP_ */
