#ifndef PROT_PRSM_PRSM_FEATURE_SPECIES_HPP_
#define PROT_PRSM_PRSM_FEATURE_SPECIES_HPP_

#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace prot {

class PrsmFeatureSpecies {
 public:
  PrsmFeatureSpecies(const std::string &db_file_name,
              const std::string &spec_file_name,
              const std::string &input_file_ext,
              const std::string &output_file_ext,
              const ModPtrVec &fix_mod_ptr_vec);
  void process();
 private:
  std::string db_file_name_;
  std::string spec_file_name_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  ModPtrVec fix_mod_ptr_vec_;

  void setProtId(PrsmPtrVec& prsm_ptrs);

  void setSpeciesId(const PrsmPtrVec& prsm_ptrs);
};

typedef std::shared_ptr<PrsmFeatureSpecies> PrsmFeatureSpeciesPtr;


} /* namespace prot */

#endif 
