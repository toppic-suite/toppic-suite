/*
 * bp_spec_reader.hpp
 *
 *  Created on: Nov 28, 2013
 *      Author: xunlikun
 */

#ifndef BP_SPEC_READER_HPP_
#define BP_SPEC_READER_HPP_

#include <bp_spec.hpp>

namespace prot {

class BpSpecReader {
public:
	static BpSpecPtrVec readDb(RSPtrVec rs_list);
};

} /* namespace prot */

#endif /* BP_SPEC_READER_HPP_ */
