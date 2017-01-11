// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#ifndef PROT_PEAK_NODE_HPP
#define PROT_PEAK_NODE_HPP

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
  PeakNode(double m): mono_mass(m), maxPrefix(0){};

  double getMass() {return mono_mass;}

  //int getCharge() {return charge;}

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
  //int charge;
  int componentId;
  std::vector<PeakNodePtr> next, prev;
  int maxPrefix;
};

}
#endif
