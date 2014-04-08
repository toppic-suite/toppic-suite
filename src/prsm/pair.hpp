/*
 * pair.hpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

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
  Pair(int x,int y){
    x_=x; y_=y;
  }

  int getX(){return x_;}

  int getY(){return y_;}

  void setX(int x){x_=x;}

  void setY(int y){y_=y;}

 protected:
  int x_=0;
  int y_=0;
};

inline bool comparePairUp(PairPtr c1, PairPtr c2) {
  if(c1->getY()!= c2->getY()){
    return c1->getY()<c2->getY();
  }
  return c1->getX() < c2->getX();
}

std::vector<double> compPpoDeviation(ExtendMsPtr ms,TheoPeakPtrVec ions,
                                     double ppo);
double compIonScore(ExtendMsPtr ms,TheoPeakPtrVec ions,double recal,
                    double ppo);

PeakIonPairPtrVec findPairs(ExtendMsPtr ms,TheoPeakPtrVec ions,int bgn,
                            int end);

std::vector<double> getNCScore(ExtendMsPtr ms,TheoPeakPtrVec ions,
                               int bgn,int end,double delta,double ppo);
} /* namespace prot */

#endif /* PAIR_HPP_ */
