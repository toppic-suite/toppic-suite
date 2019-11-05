//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef TOPPIC_ONE_PTM_SEARCH_PTM_SEARCH_PAIR_HPP_
#define TOPPIC_ONE_PTM_SEARCH_PTM_SEARCH_PAIR_HPP_

#include <memory>

namespace toppic {

class Pair;

typedef std::shared_ptr<Pair> PairPtr;

class Pair {
 public:
  Pair(int x, int y): x_(x), y_(y) {}

  int getX() {return x_;}

  int getY() {return y_;}

  void setX(int x) {x_ = x;}

  void setY(int y) {y_ = y;}

  static bool cmpPosInc(const PairPtr &a, const PairPtr &b) {
    if (a->getY() != b->getY()) {
      return a->getY() < b->getY();
    }
    return a->getX() < b->getX();
  }

 protected:
  int x_ = 0;

  int y_ = 0;
};


} /* namespace toppic */

#endif /* PAIR_HPP_ */
