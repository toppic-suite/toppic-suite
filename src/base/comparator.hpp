/*
 * mass_comparator.hpp
 *
 *  Created on: Nov 26, 2013
 *      Author: xunlikun
 */

#ifndef PROT_MASS_COMPARATOR_HPP_
#define PROT_MASS_COMPARATOR_HPP_

#include "base/theo_peak.hpp"

namespace prot {
class Comparator{
public:
 /*
	static bool mass_up(const double p, const double n){
		return p<n;
	}
	static bool mass_down(const double p, const double n){
		return p>n;
	}
  */

	static bool theopeak_up(const TheoPeakPtr p,TheoPeakPtr n){
		return p->getPosition() < n->getPosition();
	}
	static bool theopeak_down(const TheoPeakPtr p,TheoPeakPtr n){
		return p->getPosition() > n->getPosition();
	}

};
}

#endif /* MASS_COMPARATOR_HPP_ */
