/*
 * mass_comparator.hpp
 *
 *  Created on: Nov 26, 2013
 *      Author: xunlikun
 */

#ifndef PROT_MASS_COMPARATOR_HPP_
#define PROT_MASS_COMPARATOR_HPP_

namespace prot {
class Comparator{
public:
	static bool mass_up(const double p, const double n){
		return p<n;
	}
	static bool mass_down(const double p, const double n){
		return p>n;
	}
};
}

#endif /* MASS_COMPARATOR_HPP_ */
