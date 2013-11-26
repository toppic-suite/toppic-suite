/*
 * ion.hpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#ifndef PROT_ION_H_
#define PROT_ION_H_

#include <sstream>
#include "ion_type.hpp"
#include "neutral_loss.hpp"

namespace prot {

class Ion {
public:
	Ion(int charge,int pos,int display,IonTypePtr ion_type_ptr,NeutralLossPtr neutral_loss_ptr);
	int getCharge(){return charge_;}
	int getPos(){return pos_;}
	int getDisplayPos(){return display_pos_;}
	IonTypePtr getIonType(){return ion_type_;}
	std::string getDisplayName(){
		std::stringstream s;
		s<<display_pos_;
		return ion_type_ -> getName() + s.str();
	}
private:
	int charge_;
	int pos_;
	int display_pos_;
	IonTypePtr ion_type_;
	NeutralLossPtr neutral_loss_;
};

typedef std::shared_ptr<Ion> IonPtr;
typedef std::vector<IonPtr> IonPtrVec;

} /* namespace prot */

#endif /* ION_HPP_ */
