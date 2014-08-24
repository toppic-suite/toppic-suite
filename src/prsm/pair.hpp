#ifndef PAIR_HPP_
#define PAIR_HPP_

#include "base/algorithm.hpp"
#include "spec/extend_peak.hpp"
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"

namespace prot {

class Pair;

typedef std::shared_ptr<Pair> PairPtr;

class Pair {
 public:
  Pair(int x,int y);

  int getX(){return x_;}

  int getY(){return y_;}

  void setX(int x){x_=x;}

  void setY(int y){y_=y;}

 protected:
  int x_=0;
  int y_=0;
};

inline bool comparePairUp(const PairPtr &a, const PairPtr &b) {
  if(a->getY() != b->getY()){
    return a->getY() < b->getY();
  }
  return a->getX() < b->getX();
}

std::vector<double> compPpoDeviation(ExtendMsPtr ms_ptr, const TheoPeakPtrVec &peak_ptrs,
                                     double ppo);

double compIonScore(ExtendMsPtr ms_ptr, const TheoPeakPtrVec &peak_ptrs,double recal,
                    double ppo);

// peak_ptrs are sorted with masses
PeakIonPairPtrVec findPairs(ExtendMsPtr ms_ptr, TheoPeakPtrVec peak_ptrs,
                            int bgn, int end, double add_tolerance);

std::vector<double> getNCScore(ExtendMsPtr ms_ptr, const TheoPeakPtrVec &peak_ptrs,
                               int bgn,int end,double delta,double ppo);
} /* namespace prot */

#endif /* PAIR_HPP_ */
