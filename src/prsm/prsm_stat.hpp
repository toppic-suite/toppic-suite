#ifndef PROT_PRSM_PRSM_STAT_HPP_
#define PROT_PRSM_PRSM_STAT_HPP_

#include <string>

#include "base/string_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class PrsmStat {
 public:
  PrsmStat(PrsmParaPtr prsm_para_ptr, 
           const std::string &input_file_ext, 
           const std::string &output_file_ext);
  void process();

  void writePrsm(std::ofstream &file, PrsmPtr prsm_ptr);

 private:
  PrsmParaPtr prsm_para_ptr_;
  double min_mass_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  AcidPtrVec acid_ptr_vec_;
};

typedef std::shared_ptr<PrsmStat> PrsmStatPtr;


} /* namespace prot */

#endif /* PRSM_STAT_HPP_ */
