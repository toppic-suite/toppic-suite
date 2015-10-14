#ifndef PROT_GRAPH_ALIGN_PROCESSOR_HPP_
#define PROT_GRAPH_ALIGN_PROCESSOR_HPP_

#include "graph/graph.hpp"
#include "graphalign/graph_align_mng.hpp"

namespace prot {

class GraphAlignProcessor {
 public:
  GraphAlignProcessor(GraphAlignMngPtr mng_ptr);
  void process();
 private:
  GraphAlignMngPtr mng_ptr_;
};

typedef std::shared_ptr<GraphAlignProcessor> GraphAlignProcessorPtr;

}

#endif

