/*
 * BreakPoint.hpp
 *
 *  Created on: Nov 26, 2013
 *      Author: xunlikun
 */

#ifndef PROT_BREAK_POINT_HPP_
#define PROT_BREAK_POINT_HPP_

#include "ion_type.hpp"

namespace prot {

class BreakPoint {
public:
	BreakPoint(double prm,double srm);
	double getPrm(){return prm_;}
	double getSrm(){return srm_;}
	double getNTermMass(IonTypePtr ion_type){return prm_ + ion_type->getShift();}
	double getCTermMass(IonTypePtr ion_type){return srm_ + ion_type->getShift();}
private:
	double prm_;
	double srm_;
};

typedef std::shared_ptr<BreakPoint> BreakPointPtr;
typedef std::vector<BreakPointPtr> BreakPointPtrVec;

} /* namespace prot */

#endif /* BREAKPOINT_HPP_ */
