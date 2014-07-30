
#include "base/species.hpp"

namespace prot {

Species::Species(ProteoformPtr proteoform_ptr){
  id_=0;
  addProteoform(proteoform_ptr);
}

void Species::addProteoform(ProteoformPtr proteoform_ptr){
  proteoform_ptr_vec_.push_back(proteoform_ptr);
}

void Species::setSpeciesId(int id){
  id_=id;
  for(size_t i=0;i< proteoform_ptr_vec_.size();i++){
    proteoform_ptr_vec_[i]->setSpeciesId(id);
  }
}

ProteoformPtr Species::getFirstProteoform(){
  return proteoform_ptr_vec_[0];
}

} /* namespace prot */
