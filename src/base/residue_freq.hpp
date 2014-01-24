/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_RESIDUE_FREQ_HPP_
#define PROT_RESIDUE_FREQ_HPP_

#include "base/residue.hpp"

namespace prot {

class ResidueFreq: public Residue {
 public:
  ResidueFreq(std::string acid_name, std::string abbr_name, 
              double freq);

  double getFreq() {return freq_;}

 private:
  double freq_;
};

typedef std::shared_ptr<ResidueFreq> ResFreqPtr;
typedef std::vector<ResFreqPtr> ResFreqPtrVec;

ResFreqPtrVec getResidueFreqPtrVecInstance(AcidPtrVec &acid_list, PtmPtrVec &ptm_list,
                                           std::string &file_name); 

}
#endif
