#ifndef PROT_PRSM_PAIR_HPP_
#define PROT_PRSM_PAIR_HPP_

#include <memory>

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

  static bool cmpYIncXInc(const PairPtr &a, const PairPtr &b);

 protected:
  int x_=0;
  int y_=0;
};

} /* namespace prot */

#endif /* PAIR_HPP_ */
