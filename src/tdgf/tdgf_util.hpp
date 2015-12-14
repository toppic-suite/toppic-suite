#ifndef PROT_TDGF_UTIL_HPP_
#define PROT_TDGF_UTIL_HPP_

#include "base/proteoform.hpp"
#include "base/residue_freq.hpp"
#include "tdgf/tdgf_mng.hpp"

namespace prot {

class TdgfUtil {
 public:
  static void updateNTermResidueCounts(ResiduePtrVec &residue_list, 
                                       std::vector<double> &counts,
                                       const ProteoformPtrVec &mod_proteo_ptrs);


  static void updateResidueCounts(const ResiduePtrVec &residue_list, 
                                  std::vector<double> &counts,
                                  ProteoformPtr prot_ptr);

  static ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list, 
                                       const std::vector<double> &counts);

  static int computeAvgLength(const ResFreqPtrVec &residue_ptrs, double convert_ratio);
};

}

#endif
