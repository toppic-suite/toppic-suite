//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEARCH_DIAG_DIAGONAL_HPP_
#define TOPPIC_SEARCH_DIAG_DIAGONAL_HPP_

#include "search/diag/diag_header.hpp"
#include "search/diag/diag_pair.hpp"

namespace toppic {

class Diagonal;

typedef std::shared_ptr<Diagonal> DiagonalPtr;
typedef std::vector<DiagonalPtr>  DiagonalPtrVec;

class Diagonal {
  public:
    Diagonal() {}

    // Add contruction method here because a template 
    // is used. 
    explicit Diagonal(DiagHeaderPtr header_ptr):
      header_ptr_(header_ptr) {}

    // need init pair_ptr_list after create
    explicit Diagonal(DiagHeaderPtr header_ptr, 
                      DiagPairPtrVec pair_ptr_list):
      header_ptr_(header_ptr),
      pair_ptr_list_(pair_ptr_list) {}

    size_t size() {return pair_ptr_list_.size(); }

    DiagHeaderPtr getHeader() {return header_ptr_;}

    const DiagPairPtrVec& getDiagPairPtrVec() {return pair_ptr_list_;}

    DiagPairPtr getDiagPair(int i) {return pair_ptr_list_[i];}

  private:
    DiagHeaderPtr header_ptr_;

    DiagPairPtrVec pair_ptr_list_;
};

} /* namespace toppic */

#endif /* DIAGONAL_HPP_ */
