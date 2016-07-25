#include "base/file_util.hpp"
#include "base/proteoform_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_species.hpp"

namespace prot {

PrsmSpecies::PrsmSpecies(const std::string &db_file_name,
                         const std::string &spec_file_name,
                         const std::string &input_file_ext,
                         const ModPtrVec &fix_mod_ptr_vec,
                         const std::string &output_file_ext,
                         double ppo): 
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    fix_mod_ptr_vec_(fix_mod_ptr_vec),
    output_file_ext_(output_file_ext),
    ppo_(ppo) {
    }

ProteoformPtrVec2D PrsmSpecies::groupProteins(const PrsmPtrVec &prsm_ptrs){
  //get max shift number
  int max_shift_number = 0;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    int cur_shift_number 
        = prsm_ptrs[i]->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED);
    if(max_shift_number < cur_shift_number) {
      max_shift_number = cur_shift_number;
    }
  }
  //get proteoform groups
  ProteoformPtrVec2D proteogroups;
  for(int shift =0;shift<=max_shift_number;shift++){
    ProteoformPtrVec proteo_ptrs;
    for(size_t i=0;i<prsm_ptrs.size();i++ ){
      if(shift == prsm_ptrs[i]->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED)) {
        proteo_ptrs.push_back(prsm_ptrs[i]->getProteoformPtr());
      }
    }
    proteogroups.push_back(proteo_ptrs);
  }
  return proteogroups;
}

ProteoformPtrVec2D PrsmSpecies::getZeroPtmList(const ProteoformPtrVec& proteo_ptrs, double ppo){
  ProteoformPtrVec2D species;
  for(size_t i=0;i<proteo_ptrs.size();i++){
    bool is_found = false;
    for(size_t j=0; j<species.size(); j++){
      if(ProteoformUtil::isSameSeqAndMass(proteo_ptrs[i], species[j][0],ppo)){
        species[j].push_back(proteo_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if(!is_found){
      ProteoformPtrVec new_species;
      new_species.push_back(proteo_ptrs[i]);
      species.push_back(new_species);
    }
  }
  return species;
}

void PrsmSpecies::setProtId(PrsmPtrVec& prsm_ptrs){
  PrsmPtrVec2D proteins;
  std::vector<std::string> protein_names;
  for(size_t i=0;i<prsm_ptrs.size();i++) {
    std::string name = prsm_ptrs[i]->getProteoformPtr()->getSeqName();
    bool is_found = false;
    for(size_t j=0; j<protein_names.size(); j++){
      if(protein_names[j] == name) {
        proteins[j].push_back(prsm_ptrs[i]);
        is_found = true;
        break;
      }
    }
    if(!is_found){
      PrsmPtrVec new_protein;
      new_protein.push_back(prsm_ptrs[i]);
      proteins.push_back(new_protein);
      protein_names.push_back(name);
    }
  }

  for(size_t i=0; i<proteins.size();i++){
    for (size_t j = 0; j < proteins[i].size(); j++) {
      proteins[i][j]->getProteoformPtr()->setProtId(i);
    }
  }
}

void PrsmSpecies::setSpeciesId(PrsmPtrVec& prsm_ptrs,double ppo){
  ProteoformPtrVec2D proteo_groups = groupProteins(prsm_ptrs);
  
  // find zero ptm species 
  ProteoformPtrVec2D species = getZeroPtmList(proteo_groups[0],ppo);

  for(size_t i=1; i<proteo_groups.size(); i++){
    for(size_t j=0; j<proteo_groups[i].size();j++){
      bool is_found = false;
      for(size_t m = 0; m< species.size(); m++){
        if(ProteoformUtil::isStrictCompatiablePtmSpecies(
                proteo_groups[i][j], species[m][0], ppo)){
          species[m].push_back(proteo_groups[i][j]);
          is_found = true;
          break;
        }
      }
      if(!is_found){
        ProteoformPtrVec new_species;
        new_species.push_back(proteo_groups[i][j]);
        species.push_back(new_species);
      }
    }
  }

  for(size_t i=0; i<species.size();i++){
    for (size_t j = 0; j < species[i].size(); j++) {
      species[i][j]->setSpeciesId(i);
    }
  }
}

void PrsmSpecies::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;
  
  PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name_,
                                                  fix_mod_ptr_vec_);
  //sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);
  sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpMatchFragmentDecMatchPeakDec);
  setProtId(prsm_ptrs);
  setSpeciesId(prsm_ptrs,ppo_);
  sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);
  //output
  std::string output_file_name = base_name +"."+output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}


} /* namespace prot */


