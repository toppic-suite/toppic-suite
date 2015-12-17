
#ifndef PROT_DIAGONAL_HPP_
#define PROT_DIAGONAL_HPP_

#include <memory>
#include "spec/theo_peak.hpp"
#include "spec/extend_ms.hpp"
#include "prsm/peak_ion_pair.hpp"
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
                                     double min_mass,
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
