/*
 * species.cpp
 *
 *  Created on: Feb 20, 2014
 *      Author: xunlikun
 */

#include <base/species.hpp>

namespace prot {
Species::Species(ProteoformPtr proteoform){
  id_=0;
  addProteoform(proteoform);
}
 void Species::addProteoform(ProteoformPtr proteoform){
   proteoforms_.push_back(proteoform);
 }
void Species::setSpeciesId(int id){
  id_=id;
  for(unsigned int i=0;i<proteoforms_.size();i++){
    proteoforms_[i]->setSpeciesId(id);
  }
}
ProteoformPtr Species::getFistProteoform(){
  return proteoforms_[0];
}

ProteoformPtrVec2D groupProteins(PrSMPtrVec& prsms){
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

SpeciesPtrVec getZeroPtmList(ProteoformPtrVec& proteoforms,double ppo){
  SpeciesPtrVec list;
  for(unsigned int i=0;i<proteoforms.size();i++){
    for(unsigned int j=0;j<list.size();j++){
      if(isSamePeptideAndMass(proteoforms[i],list[j]->getFistProteoform(),ppo)){
        list[j]->addProteoform(proteoforms[i]);
        break;
      }
    }
    list.push_back(SpeciesPtr(new Species(proteoforms[i])));
  }
  return list;
}

SpeciesPtrVec setSpeciesId(PrSMPtrVec& prsms,double ppo){
  ProteoformPtrVec2D proteogroups = groupProteins(prsms);

  SpeciesPtrVec list = getZeroPtmList(proteogroups[0],ppo);

  for(unsigned int i=1;i<proteogroups.size();i++){
    for(unsigned int j=0;j<proteogroups[i].size();j++){
      for(unsigned int m = 0; m<list.size();m++){
        if(isStrictCompatiablePtmSpecies(proteogroups[i][j],list[m]->getFistProteoform(),ppo)){
          list[m]->addProteoform(proteogroups[i][j]);
          break;
        }
      }
      list.push_back(SpeciesPtr(new Species(proteogroups[i][j])));
    }
  }

  for(unsigned int i=0;i<list.size();i++){
    list[i]->setSpeciesId(i);
  }
  return list;
}
} /* namespace prot */
