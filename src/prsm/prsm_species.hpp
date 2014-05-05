
#include "base/species.hpp"
#include "prsm/prsm.hpp"
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
  PrsmSpecies(std::string db_file,
              std::string spec_file,
              std::string input_file,
              std::string output_file,
              double ppo);

  PrsmSpecies(std::map<std::string,std::string> arguments,
              std::string input_file,
              std::string output_file);
  void process();
 private:
  std::string db_file_;
  std::string spec_file_;
  std::string input_file_;
  std::string output_file_;
  double ppo_;
};

SpeciesPtrVec setSpeciesId(const PrsmPtrVec &prsms,double ppo);

} /* namespace prot */

#endif /* PRSM_SPECIES_HPP_ */
