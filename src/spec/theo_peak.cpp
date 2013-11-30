/*
 * theo_peak.cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#include "theo_peak.hpp"

namespace prot {
TheoPeak::TheoPeak(IonPtr ion,double unmod_mass,double shift):
  Peak(unmod_mass + shift, 1.0) {
	ion_ = ion;
	unmod_mass_ = unmod_mass;
	shift_ = shift;
}
} /* namespace prot */
