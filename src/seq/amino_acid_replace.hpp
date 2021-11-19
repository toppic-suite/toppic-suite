//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_SEQ_AMINO_ACID_REPLACE_HPP_
#define TOPPIC_SEQ_AMICO_ACID_REPLACE_HPP_

#include "common/base/amino_acid.hpp"

namespace toppic {

class AminoAcidReplace {
 public:
  AminoAcidReplace(std::string ori_letter, AminoAcidPtr amino_acid_ptr, int pos); 

  std::string getOriLetter() {return ori_letter_;}

  AminoAcidPtr getAminoAcidPtr() {return amino_acid_ptr_;} 

  int getPos() {return pos_;}

 private:
  std::string ori_letter_;
  AminoAcidPtr amino_acid_ptr_;
  int pos_;
};

typedef std::shared_ptr<AminoAcidReplace> AminoAcidReplacePtr;
typedef std::vector<AminoAcidReplacePtr> AminoAcidReplacePtrVec;

}  // namespace toppic

#endif
