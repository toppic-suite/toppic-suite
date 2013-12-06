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
	SpecData(char const* config_file_name);
	PrmPeakTypePtrVec& getPrmPeakTypePtrVec(){return prm_peak_type_list_;}
	SupportPeakTypePtrVec& getSupportPeakTypePtrVec(){return support_peak_type_list_;}
private:
	PrmPeakTypePtrVec prm_peak_type_list_;
	SupportPeakTypePtrVec support_peak_type_list_;
};

typedef std::shared_ptr<SpecData> SpecDataPtr;

} /* namespace prot */

#endif /* SPEC_DATA_HPP_ */
