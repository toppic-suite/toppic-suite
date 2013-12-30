#ifndef PROT_DB_RESIDUE_SEQ_HPP_
#define PROT_DB_RESIDUE_SEQ_HPP_

#include <string>

#include "base/residue_seq.hpp"

namespace prot {

class DbResidueSeq: public ResidueSeq {
 public:
  DbResidueSeq(ResiduePtrVec residues, int id, std::string name): ResidueSeq(residues) {
    id_ = id;
    name_ = name;
  }

  int getId() {return id_;}
  std::string getName() {return name_;}

 private:
  int id_;
  std::string name_;
};

typedef std::shared_ptr<DbResidueSeq> DbResSeqPtr;

}
#endif
