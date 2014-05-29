/*
 * diagonal.hpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#ifndef PROT_DIAGONAL_HPP_
#define PROT_DIAGONAL_HPP_

#include <memory>
#include "spec/theo_peak.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/pair.hpp"
#include "ptmsearch/diagonal_header.hpp"

namespace prot {

template <class T>
class Diagonal{
 public:
  Diagonal(){};
  Diagonal(DiagonalHeaderPtr header){
    header_ = header;
  };
  /**
   * need init pair_ptr_list after create
   */
  Diagonal(DiagonalHeaderPtr header,std::vector<T> pair_ptr_list){
    header_ = header;
    pair_ptr_list_ = pair_ptr_list;
  };

  unsigned int size(){
    return pair_ptr_list_.size();
  }
  DiagonalHeaderPtr getHeader(){
    return header_;
  }
  std::vector<T> getDiagPair(){
    return pair_ptr_list_;
  }

  T getDiagPair(int i){
    return pair_ptr_list_[i];
  }

 private:
  DiagonalHeaderPtr header_;
  std::vector<T> pair_ptr_list_;
};

DiagonalHeaderPtrVec refineHeadersBgnEnd(
    int first_res_pos,
    int last_res_pos,
    ProteoformPtr seq,
    DeconvMsPtr deconv_ms,
    ExtendMsPtr ms_three,
    PtmMngPtr mng,
    DiagonalHeaderPtrVec headers);
int getNewBgn(PeakIonPairPtrVec pairs);
int getNewEnd(PeakIonPairPtrVec pairs);
} /* namespace prot */

#endif /* DIAGONAL_HPP_ */
