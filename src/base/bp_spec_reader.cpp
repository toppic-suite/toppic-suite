/*
 * bp_spec_reader..cpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#include <bp_spec_reader.hpp>

namespace prot {
BpSpecPtrVec BpSpecReader::readDb(RSPtrVec rs_list){
	BpSpecPtrVec bpspec_ptr_list;
	for(unsigned int i =0; i<rs_list.size();i++){
		bpspec_ptr_list.push_back(BpSpecPtr(new BpSpec(rs_list[i])));
	}
	return bpspec_ptr_list;
}
} /* namespace prot */
