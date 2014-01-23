/*
 * ptm_fastfilter_mng.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include "ptm_fast_filter_mng.hpp"

namespace prot {
PtmFastFilterMng::PtmFastFilterMng(std::string config_file_name){
	base_data = BaseDataPtr(new BaseData("conf/configuration.xml"));
}

} /* namespace tools */
