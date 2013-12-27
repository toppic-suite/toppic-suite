/*
 * anno_cleavage.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#ifndef ANNO_CLEAVAGE_HPP_
#define ANNO_CLEAVAGE_HPP_
#include <memory>
#include <vector>

namespace prot {

class AnnoCleavage {
};
typedef std::shared_ptr<AnnoCleavage> AnnoCleavagePtr;
typedef std::vector<AnnoCleavagePtr> AnnoCleavagePtrVec;

} /* namespace prot */

#endif /* ANNO_CLEAVAGE_HPP_ */
