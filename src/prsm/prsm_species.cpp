/*
 *  Created on: Feb 19, 2014
 */

#include "base/species.hpp"
#include "base/file_util.hpp"
#include "prsm/prsm_species.hpp"

namespace prot {
PrsmSpecies::PrsmSpecies(std::string db_file,
                 std::string spec_file,
                 std::string input_file,
                 std::string output_file,
                 double ppo) {
  db_file_= db_file;
  spec_file_ = spec_file;
  input_file_=input_file;
  output_file_ = output_file;
  ppo_ = ppo;
}

PrsmSpecies::PrsmSpecies(std::map<std::string,std::string> arguments,
                 std::string input_file,
                 std::string output_file) {
  db_file_= arguments["databaseFileName"];
  spec_file_ = arguments["spectrumFileName"];
  input_file_=input_file;
  output_file_ = output_file;
  ppo_ = atof(arguments["ppo"].c_str());
}

void PrsmSpecies::process(){
  std::string base_name = basename(spec_file_);
  std::string input_file_name = base_name+"."+input_file_;
  ProteoformPtrVec proteoforms_ = prot::readFastaToProteoform(db_file_,ResidueFactory::getBaseResiduePtrVec());
  PrSMPtrVec prsms = readPrsm(input_file_name,proteoforms_);
  sort(prsms.begin(),prsms.end(),prsmSpectrumIdUpMatchFragUp);
  setSpeciesId(prsms,ppo_);

  //output
  std::string output_file_name = base_name +"."+output_file_;
  PrSMWriter writer(output_file_name);
  writer.writeVector(prsms);
  //because the prsm_writer ~PrSMWriter changed and the fileclosing is an independant function
  writer.close();
}
} /* namespace prot */
