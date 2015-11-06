#include "base/align_type.hpp"

namespace prot {

AlignTypePtr AlignType::COMPLETE = AlignTypePtr(new AlignType("COMPLETE", 0));
AlignTypePtr AlignType::PREFIX = AlignTypePtr(new AlignType("PREFIX", 1));
AlignTypePtr AlignType::SUFFIX  = AlignTypePtr(new AlignType("SUFFIX", 2));
AlignTypePtr AlignType::INTERNAL = AlignTypePtr(new AlignType("INTERNAL", 3));

AlignType::AlignType(const std::string &name, int id):
    name_(name), 
    id_(id) {
    }

}
