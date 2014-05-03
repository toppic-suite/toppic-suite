/*
 * species.cpp
 *
 *  Created on: Feb 20, 2014
 *      Author: xunlikun
 */

#include "base/species.hpp"

namespace prot {

Species::Species(const ProteoformPtr &proteoform){
  id_=0;
  addProteoform(proteoform);
}

void Species::addProteoform(const ProteoformPtr &proteoform){
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

} /* namespace prot */
