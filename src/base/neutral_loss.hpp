/*
 * NeutralLoss.hpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#ifndef PROT_NEUTRAL_LOSS_HPP_
#define PROT_NEUTRAL_LOSS_HPP_

#include <string>
#include <vector>
#include <memory>

namespace prot {

class NeutralLoss {

public:
	NeutralLoss(std::string name,double mass);
	std::string getName(){return name_;}
	double getMass(){return mass_;}
private:
	std::string name_;
	double mass_;
};

typedef std::shared_ptr<NeutralLoss> NeutralLossPtr;
typedef std::vector<NeutralLossPtr> NeutralLossPtrVec;

NeutralLossPtrVec getNeutralLossPtrVecInstance(const std::string &file_name);
NeutralLossPtr getNeutralLossPtrByName(NeutralLossPtrVec &neutralLoss_ptr_vec, 
                                       const std::string &name);

} /* namespace prot */

#endif /* NEUTRALLOSS_HPP_ */
