
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
                                     PtmMngPtr mng_ptr);

DiagonalHeaderPtrVec refineHeadersBgnEnd(
    ProteoformPtr proteo_ptr,
    const ExtendMsPtrVec &ms_three_ptr_vec,
    const DiagonalHeaderPtrVec& heade_ptrs,
    PtmMngPtr mng_ptr);

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs);
int getNewEnd(const PeakIonPairPtrVec& pair_ptrs);

} /* namespace prot */

#endif /* DIAGONAL_HPP_ */
