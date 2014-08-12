
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
    header_ptr = header_ptr;
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

DiagonalHeaderPtrVec refineHeadersBgnEnd(
    int first_res_pos,
    int last_res_pos,
    ProteoformPtr proteo_ptr,
    DeconvMsPtr deconv_ms_ptr,
    ExtendMsPtr ms_three_ptr,
    PtmMngPtr mng_ptr,
    const DiagonalHeaderPtrVec& heade_ptrs);

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs);
int getNewEnd(const PeakIonPairPtrVec& pair_ptrs);

} /* namespace prot */

#endif /* DIAGONAL_HPP_ */
