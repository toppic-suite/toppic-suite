#ifndef PROT_SPEC_MULTI_SPEC_PEAK_HPP_
#define PROT_SPEC_MULTI_SPEC_PEAK_HPP_

#include <memory>
#include <vector>
#include "spec/peak.hpp"

namespace prot {

class MultiSpecPeak;
typedef std::shared_ptr<MultiSpecPeak> MultiSpecPeakPtr;

class MultiSpecPeak : public Peak {
 public:
  MultiSpecPeak (double mono_mass, double intensity, int spec_id, 
                 int in_spec_peak_id, int base_type, double score);

  int getId() {return id_;}

  int getSpecId() {return spec_id_;}

  int getInSpecPeakId() {return in_spec_peak_id_;}

  int getBaseType() {return base_type_;}

  double getMonoMass() {return getPosition();}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

  static bool cmpPosIncrease(const MultiSpecPeakPtr &a, const MultiSpecPeakPtr &b){
    return a->getPosition() < b->getPosition();
  }

 private:
  int id_;
  int spec_id_;
  int in_spec_peak_id_;
  int base_type_;
  double score_;
};

typedef std::vector<MultiSpecPeakPtr> MultiSpecPeakPtrVec;

}
#endif
