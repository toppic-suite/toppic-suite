
#ifndef SPECIES_HPP_
#define SPECIES_HPP_

#include "base/proteoform.hpp"

namespace prot {

class Species {
 public:
  Species(ProteoformPtr proteoform_ptr);
  void addProteoform(ProteoformPtr proteoform_ptr);
  void setSpeciesId(int id);
  ProteoformPtr getFirstProteoform();
 private:
  int id_;
  ProteoformPtrVec proteoform_ptr_vec_;
};

typedef std::shared_ptr<Species> SpeciesPtr;
typedef std::vector<SpeciesPtr> SpeciesPtrVec;

} /* namespace prot */

#endif /* SPECIES_HPP_ */
