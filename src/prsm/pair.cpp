#include "prsm/pair.hpp"

namespace prot {

Pair::Pair (int x, int y):
    x_(x),
    y_(y) {
    }

bool Pair::cmpYIncXInc(const PairPtr &a, const PairPtr &b) {
  if(a->getY() != b->getY()){
    return a->getY() < b->getY();
  }
  return a->getX() < b->getX();
}

}
