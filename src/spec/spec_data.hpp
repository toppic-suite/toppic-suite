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
public :
	static PrmPeakTypePtrVec& getPrmPeakTypePtrVec(){return prm_peak_type_list_;}
	static SupportPeakTypePtrVec& getSupportPeakTypePtrVec(){return support_peak_type_list_;}
private:
	static PrmPeakTypePtrVec prm_peak_type_list_;
	static SupportPeakTypePtrVec support_peak_type_list_;
};

} /* namespace prot */

#endif /* SPEC_DATA_HPP_ */
