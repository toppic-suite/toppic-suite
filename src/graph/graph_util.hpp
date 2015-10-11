#ifndef PROT_GRAPH_UTIL_HPP_
#define PROT_GRAPH_UTIL_HPP_

#include <vector>

#include "graph/graph.hpp"

namespace prot {

void writeToDot(const std::string &file_name, MassGraphPtr graph_ptr);

}

#endif
