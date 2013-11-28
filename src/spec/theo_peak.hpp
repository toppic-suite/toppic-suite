/*
 * theo_peak.hpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#ifndef THEO_PEAK_HPP_
#define THEO_PEAK_HPP_

#include "ion.hpp"
#include "peak.hpp"

namespace prot {

class TheoPeak : public Peak {
public:
	TheoPeak(IonPtr ion,double unmode_mass,double shift);
	IonPtr getIon(){return ion_;}
	double getModMass(){return getPosition();}
	double getShift(){return shift_;}
private:
	double unmode_mass_;
	double shift_;

	IonPtr ion_;
};

} /* namespace prot */

#endif /* THEO_PEAK_HPP_ */
