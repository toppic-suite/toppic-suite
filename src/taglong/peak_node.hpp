// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef PROT_TAG_LONG_PEAK_NODE_HPP
#define PROT_TAG_LONG_PEAK_NODE_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

namespace prot {

class PeakNode;

typedef std::shared_ptr<PeakNode> PeakNodePtr;

class PeakNode {
 public:
  explicit PeakNode(double m): mono_mass(m), maxPrefix(0) {}

  double getMass() {return mono_mass;}

  int getMaxPrefix() {return maxPrefix;}

  void setMaxPrefix(int m) {maxPrefix = m;}

  int getComponenetId() {return componentId;}

  void setComponenetId(int id) {
    if (componentId != id) {
      componentId = id;
      for (size_t i = 0; i < next.size(); i++) {
        next[i]->setComponenetId(id);
      }
    }
  }

  std::vector<PeakNodePtr> getNext() {return next;}

  void addNext(PeakNodePtr peak) {
    next.push_back(peak);
    if (peak->getComponenetId() != componentId) {
      doUpdateComponentId(peak);
    }
  }

  double diff(PeakNodePtr prev) {return mono_mass - prev->getMass();}

  void clearEdges() {
    next.clear();
    prev.clear();
  }

  bool updateComponentId() {
    for (size_t i = 0; i < next.size(); i++) {
      if (next[i]->getComponenetId() != componentId) {
        doUpdateComponentId(next[i]);
        return true;
      }
    }
    return false;
  }

 private:
  void doUpdateComponentId(PeakNodePtr peak) {
    int min_id = std::min(componentId, peak->getComponenetId());
    setComponenetId(min_id);
    peak->setComponenetId(min_id);
  }

  double mono_mass;

  int componentId;

  std::vector<PeakNodePtr> next, prev;

  int maxPrefix;
};

}  // namespace prot
#endif
