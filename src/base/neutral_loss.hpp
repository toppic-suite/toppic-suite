/*
 * NeutralLoss.hpp
 *
 *  Created on: Nov 25, 2013
 *      Author: xunlikun
 */

#ifndef PROT_NEUTRALLOSS_H_
#define PROT_NEUTRALLOSS_H_

#include <string>
#include <vector>
#include <memory>
#include <xercesc/dom/DOM.hpp>

namespace prot {

class NeutralLoss {

public:
	NeutralLoss(std::string name,double mass);
	NeutralLoss(xercesc::DOMElement * element);
	std::string getName(){return name_;}
	double getMass(){return mass_;}
private:
	std::string name_;
	double mass_;
};

typedef std::shared_ptr<NeutralLoss> NeutralLossPtr;
typedef std::vector<NeutralLossPtr> NeutralLossPtrVec;

NeutralLossPtrVec getNeutralLossPtrVecInstance(const char* file_name);
NeutralLossPtr getNeutralLossPtrByName(NeutralLossPtrVec &neutralLoss_ptr_vec, const std::string &name);

} /* namespace prot */

#endif /* NEUTRALLOSS_HPP_ */
