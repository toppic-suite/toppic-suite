#ifndef __PEAK_NODE_HPP__
#define __PEAK_NODE_HPP__

#include "spec/spectrum_set.hpp"

namespace prot{

class PeakNode;

typedef std::shared_ptr<PeakNode> PeakNodePtr;

class PeakNode {
 public:
  PeakNode(PrmPeakPtr peak):peak_(peak),maxPrefix_(0){};

  PeakNode();

  PrmPeakPtr getPrmPeak() {return peak_;}

  int getComponentId() {return componentId_;}

  void setComponenetId(int id);

  void clearEdges() {
    next_.clear();
  };

  void addNext(PeakNodePtr peak);

  std::vector<PeakNodePtr> getNext() {return next_;}

  bool updateComponentId();

  int getMaxPrefix() {return maxPrefix_;}

  void setMaxPrefix(int m) {maxPrefix_ = m;}

  double getMonoMass();

 private:

  void doUpdateComponentId(PeakNodePtr peak);

  PrmPeakPtr peak_;	
  int componentId_;
  int maxPrefix_;
  std::vector<PeakNodePtr> next_;
};

}

#endif
