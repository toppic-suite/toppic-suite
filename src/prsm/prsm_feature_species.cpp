#include "base/file_util.hpp"
#include "base/proteoform_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_feature_species.hpp"

namespace prot {

PrsmFeatureSpecies::PrsmFeatureSpecies(const std::string &db_file_name,
                                       const std::string &spec_file_name,
                                       const std::string &input_file_ext,
                                       const std::string &output_file_ext,
                                       const ModPtrVec &fix_mod_ptr_vec): 
    db_file_name_(db_file_name),
    spec_file_name_(spec_file_name),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext),
    fix_mod_ptr_vec_(fix_mod_ptr_vec) {
    }

void PrsmFeatureSpecies::setProtId(PrsmPtrVec& prsm_ptrs){
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


void PrsmFeatureSpecies::setSpeciesId(const PrsmPtrVec& prsm_ptrs){
  PrsmPtrVec2D species;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    bool is_found = false;
    PrsmPtr cur_ptr = prsm_ptrs[i];
    for(size_t j=0; j < species.size(); j++){
      PrsmPtr ref_ptr = species[j][0];
      if (cur_ptr->getProteoformPtr()->getProtId() == ref_ptr->getProteoformPtr()->getProtId()) {
        if (cur_ptr->getPrecFeatureId() == ref_ptr->getPrecFeatureId()) {
          species[j].push_back(cur_ptr);
          is_found = true;
          break;
        }
      }
    }
    if(!is_found){
      PrsmPtrVec new_species;
      new_species.push_back(prsm_ptrs[i]);
      species.push_back(new_species);
    }
  }
  for(size_t i=0; i<species.size();i++){
    for (size_t j = 0; j < species[i].size(); j++) {
      species[i][j]->getProteoformPtr()->setSpeciesId(i);
    }
  }
}

void PrsmFeatureSpecies::process(){
  std::string base_name = FileUtil::basename(spec_file_name_);
  std::string input_file_name = base_name+"."+input_file_ext_;
  
  PrsmPtrVec prsm_ptrs = PrsmReader::readAllPrsms(input_file_name, db_file_name_,
                                                  fix_mod_ptr_vec_);
  sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpMatchFragmentDecMatchPeakDec);
  setProtId(prsm_ptrs);
  setSpeciesId(prsm_ptrs);
  sort(prsm_ptrs.begin(),prsm_ptrs.end(),Prsm::cmpSpectrumIdIncPrecursorIdInc);
  //output
  std::string output_file_name = base_name +"."+output_file_ext_;
  PrsmXmlWriter writer(output_file_name);
  writer.writeVector(prsm_ptrs);
  writer.close();
}


} /* namespace prot */


