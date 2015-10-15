#ifndef PROT_PRSM_PROB_HPP_
#define PROT_PRSM_PROB_HPP_

#include <vector>
#include <string>
#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/proteoform_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrsmProb {
 public:
  PrsmProb(const std::string &db_file_name, 
           const std::string &spec_file_name, 
           const std::string &in_file_ext,
           const std::string &out_file_ext,
           double K1, double K2,
           double pref, double inte);

  void process();
 private:
  std::string spec_file_name_;
  std::string db_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  double K1_;
  double K2_;
  double pref_;
  double inte_;

};

typedef std::shared_ptr<PrsmProb> PrsmProbPtr;
} /* namespace prot */

#endif /* PROT_PRSM_PROB_HPP_ */
