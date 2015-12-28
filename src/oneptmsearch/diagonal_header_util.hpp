#ifndef PROT_DIAGONAL_HEADER_UTIL_HPP_
#define PROT_DIAGONAL_HEADER_UTIL_HPP_

#include <limits>
#include <cmath>

#include "oneptmsearch/diagonal_header.hpp"

namespace prot {

class DiagonalHeaderUtil {
 public:
  // get the header corresponding to the top left corner in the spectral grid
  static DiagonalHeaderPtr getTopLeftCornerHeader() {
    double shift = 0;
    // n_term strict; c_term nostrict; prot n_term match; prot c_term no_match
    // pep n_term no_match; pep c_term no_match
    return DiagonalHeaderPtr(
        new DiagonalHeader(shift, true, false, true, false, false, false));
  }

  static DiagonalHeaderPtr getBottomRightCornerHeader(double seq_mass,
                                                      double prec_mass) {
    double shift = prec_mass - seq_mass;
    // n term nostrict, c_term strict, prot n_term no match ; prot c_term match
    // pep n_term no match, pep c_term no match 
    return DiagonalHeaderPtr(
        new DiagonalHeader(shift, false, true, false, true, false, false));
  }

  static int findSimilarShiftPos(const std::vector<double> &shifts, double s) {
    int best_pos = -1;
    double best_diff = std::numeric_limits<double>::infinity();
    for(size_t i = 0; i < shifts.size();i++){
      if(std::abs(shifts[i] - s) < best_diff){
        best_pos = i;
        best_diff = std::abs(shifts[i] - s);
      }
    }
    return best_pos;
  }

  static bool isExistHeader(const DiagonalHeaderPtrVec &header_ptrs, double shift) {
    for(size_t i = 0; i < header_ptrs.size();i++){
      if(header_ptrs[i]->getProtNTermShift() == shift) {
        return true;
      }
    }
    return false;
  }
};

} 

#endif
