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


#ifndef PROT_DIAGONAL_HPP_
#define PROT_DIAGONAL_HPP_

#include <memory>
#include "spec/theo_peak.hpp"
#include "spec/extend_ms.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "oneptmsearch/diagonal_header.hpp"

namespace prot {

template <class T>
class Diagonal{
 public:
  Diagonal(){};
  Diagonal(DiagonalHeaderPtr header_ptr){ header_ptr_ = header_ptr; };
  /**
   * need init pair_ptr_list after create
   */
  Diagonal(DiagonalHeaderPtr header_ptr,std::vector<T> pair_ptr_list) {
    header_ptr_ = header_ptr;
    pair_ptr_list_ = pair_ptr_list;
  };

  size_t size(){return pair_ptr_list_.size(); }

  DiagonalHeaderPtr getHeader(){return header_ptr_;}

  const std::vector<T>& getDiagPair(){return pair_ptr_list_;}

  T getDiagPair(int i){return pair_ptr_list_[i];}

 private:
  DiagonalHeaderPtr header_ptr_;
  std::vector<T> pair_ptr_list_;
};

double refinePrecursorAndHeaderShift(ProteoformPtr proteo_ptr,
                                     const ExtendMsPtrVec &ms_three_ptr_vec, 
                                     DiagonalHeaderPtrVec &header_ptrs,
                                     double ppo, double min_mass,
                                     double refine_prec_step_width);

DiagonalHeaderPtrVec refineHeadersBgnEnd(
    ProteoformPtr proteo_ptr,
    const ExtendMsPtrVec &ms_three_ptr_vec,
    const DiagonalHeaderPtrVec& heade_ptrs,
    double min_mass);

DiagonalHeaderPtrVec2D refineHeadersBgnEnd(
        ProteoformPtr proteo_ptr,
        const ExtendMsPtrVec &ms_three_ptr_vec,
        const DiagonalHeaderPtrVec2D& header_ptrs_2d,
        const DiagonalHeaderPtrVec& header_ptrs_1d,
        double min_mass);

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs);
int getNewEnd(const PeakIonPairPtrVec& pair_ptrs);

} /* namespace prot */

#endif /* DIAGONAL_HPP_ */
