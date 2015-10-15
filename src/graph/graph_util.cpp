#include <fstream>

#include "base/logger.hpp"
#include "graph/graph_util.hpp"

namespace prot {

void writeToDot(const std::string &file_name, MassGraphPtr graph_ptr) {
  std::ofstream f(file_name);
  EdgeInfoWriter<MassGraph> w(*graph_ptr.get());
  std::map<std::string, std::string> graph_attr, vertex_attr, edge_attr;
  graph_attr["rankdir"] = "LR";
  LOG_DEBUG("start write graph");
  boost::write_graphviz(
      f, *graph_ptr.get(), boost::default_writer(), w,
      boost::make_graph_attributes_writer(graph_attr, vertex_attr, edge_attr));
  f.close();
  LOG_DEBUG("end write graph");
}

}

