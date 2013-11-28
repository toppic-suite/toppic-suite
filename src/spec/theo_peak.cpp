/*
 * theo_peak.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#include "theo_peak.hpp"

namespace prot {
TheoPeak::TheoPeak(IonPtr ion,double unmode_mass,double shift){
	ion_ = ion;
	unmode_mass_ = unmode_mass;
	shift_ = shift;
	setIntensity(1.0);
	setPosition(unmode_mass+shift);
}
} /* namespace prot */
