#ifndef PRSM_SPECIES_HPP_
#define PRSM_SPECIES_HPP_

#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"

namespace prot {

class PrsmSpecies {
 public:
  PrsmSpecies(const std::string &db_file_name,
              const std::string &spec_file_name,
              const std::string &input_file_ext,
              const std::string &output_file_ext, 
              const ResiduePtrVec &residue_ptv_vec,
              double ppo);
  void process();
 private:
  std::string db_file_name_;
  std::string spec_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  ResiduePtrVec residue_ptr_vec_;
  double ppo_;
};

typedef std::shared_ptr<PrsmSpecies> PrsmSpeciesPtr;

void setSpeciesId(PrsmPtrVec& prsm_ptrs,double ppo);

} /* namespace prot */

#endif /* PRSM_SPECIES_HPP_ */
