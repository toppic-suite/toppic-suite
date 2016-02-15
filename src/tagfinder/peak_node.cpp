
#include "peak_node.hpp"

namespace prot{

PeakNode::PeakNode(){
  peak_ = std::make_shared<PrmPeak>(0, nullptr, nullptr, 0.0, 0.0);
  maxPrefix_ = 0;
}

void PeakNode::setComponenetId(int id){
  if (componentId_ != id) {
    componentId_ = id;
    for (size_t i = 0; i < next_.size(); i++) {
      next_[i]->setComponenetId(id);
    }
  }
}

void PeakNode::doUpdateComponentId(PeakNodePtr peak) {
  int min_id = std::min(componentId_, peak->getComponentId());
  setComponenetId(min_id);
  peak->setComponenetId(min_id);
}

void PeakNode::addNext(PeakNodePtr peak) {
  next_.push_back(peak);
  if (peak->getComponentId() != componentId_) {
    doUpdateComponentId(peak); 
  }
}

bool PeakNode::updateComponentId() {
  for (size_t i = 0; i < next_.size(); i++) {
    if (next_[i]->getComponentId() != componentId_) {
      doUpdateComponentId(next_[i]);
      return true;
    }
  }
  return false;
}

double PeakNode::getMonoMass() {
  if (peak_ == nullptr) {
    return 0.0;
  } else {
    return peak_->getMonoMass();
  }
}

}
