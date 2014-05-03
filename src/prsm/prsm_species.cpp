/*
 *  *  Created on: Feb 19, 2014
 *   */

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
  //because the prsm_writer ~PrSMWriter changed and the
  // fileclosing is an independant function
  writer.close();
}

ProteoformPtrVec2D groupProteins(const PrSMPtrVec &prsms){
  //get max shift number
  unsigned int max_shift_number = 0;
  for(unsigned int i=0;i<prsms.size();i++){
    if(max_shift_number < prsms[i]->getProteoformPtr()->getChangePtrVec().size()){
      max_shift_number = prsms[i]->getProteoformPtr()->getChangePtrVec().size();
    }
  }
  //get proteoform groups
  ProteoformPtrVec2D proteogroups;
  for(unsigned int shift =0;shift<=max_shift_number;shift++){
    ProteoformPtrVec proteoforms;
    for(unsigned int i=0;i<prsms.size();i++ ){
      if(shift == prsms[i]->getProteoformPtr()->getChangePtrVec().size()){
        proteoforms.push_back(prsms[i]->getProteoformPtr());
      }
    }
    proteogroups.push_back(proteoforms);
  }
  return proteogroups;
}

SpeciesPtrVec getZeroPtmList(const ProteoformPtrVec& proteoforms, double ppo){
  SpeciesPtrVec list;
  for(unsigned int i=0;i<proteoforms.size();i++){
    bool is_break = false;
    for(unsigned int j=0;j<list.size();j++){
      if(isSamePeptideAndMass(proteoforms[i],list[j]->getFistProteoform(),ppo)){
        list[j]->addProteoform(proteoforms[i]);
        is_break = true;
        break;
      }
    }
    if(!is_break){
      list.push_back(SpeciesPtr(new Species(proteoforms[i])));
    }
  }
  return list;
}

SpeciesPtrVec setSpeciesId(const PrSMPtrVec& prsms,double ppo){
  ProteoformPtrVec2D proteogroups = groupProteins(prsms);

  SpeciesPtrVec list = getZeroPtmList(proteogroups[0],ppo);

  for(unsigned int i=1;i<proteogroups.size();i++){
    for(unsigned int j=0;j<proteogroups[i].size();j++){
      bool is_break = false;
      for(unsigned int m = 0; m<list.size();m++){
        if(isStrictCompatiablePtmSpecies(proteogroups[i][j],
                                         list[m]->getFistProteoform(),
                                         ppo)){
          list[m]->addProteoform(proteogroups[i][j]);
          is_break = true;
          break;
        }
      }
      if(!is_break){
        list.push_back(SpeciesPtr(new Species(proteogroups[i][j])));
      }
    }
  }

  for(unsigned int i=0;i<list.size();i++){
    list[i]->setSpeciesId(i);
  }
  return list;
}



} /* namespace prot */


