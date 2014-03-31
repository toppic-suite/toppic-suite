/*
 * species.hpp
 *
 *  Created on: Feb 20, 2014
 *      Author: xunlikun
 */

#ifndef SPECIES_HPP_
#define SPECIES_HPP_

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"

namespace prot {

class Species {
 public:
  Species(const ProteoformPtr &proteoform);
  void addProteoform(const ProteoformPtr &proteoform);
  void setSpeciesId(int id);
  ProteoformPtr getFistProteoform();
 private:
  int id_;
  ProteoformPtrVec proteoforms_;
};

typedef std::shared_ptr<Species> SpeciesPtr;
typedef std::vector<SpeciesPtr> SpeciesPtrVec;

SpeciesPtrVec setSpeciesId(const PrSMPtrVec &prsms,double ppo);
} /* namespace prot */

#endif /* SPECIES_HPP_ */
