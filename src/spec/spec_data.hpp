/*
 * spec_data.hpp
 *
 *  Created on: Dec 5, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SPEC_DATA_HPP_
#define PROT_SPEC_DATA_HPP_

#include "prm_peak_type.hpp"
#include "support_peak_type.hpp"

namespace prot {

class SpecData {
public:
	SpecData(){};
	static std::string prm_peak_type_original ;
	static std::string prm_peak_type_reversed ;
};

std::string SpecData::prm_peak_type_original = "ORIGINAL";
std::string SpecData::prm_peak_type_reversed = "REVERSED";

typedef std::shared_ptr<SpecData> SpecDataPtr;

} /* namespace prot */

#endif /* SPEC_DATA_HPP_ */
