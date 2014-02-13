/*
 * ptm_mng.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#include "ptm_mng.hpp"

namespace prot {
PtmMng::PtmMng(std::string config_file_name){
    base_data_=BaseDataPtr(new BaseData(config_file_name));
}
} /* namespace prot */
