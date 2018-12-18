//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef TOPPIC_BASE_AMINO_ACID_BASE_HPP_
#define TOPPIC_BASE_AMINO_ACID_BASE_HPP_

#include <string>
#include <unordered_map>

#include "base/amino_acid.hpp"

namespace toppic {

class AminoAcidBase {
 public:
  static void initBase(const std::string &file_name);

  static AminoAcidPtr getEmptyAminoAcidPtr() {return empty_amino_acid_ptr_;}

  // Returns an amino acid based on the the name. Returns null if the amino
  // acid name does not exist.
  static AminoAcidPtr getAminoAcidPtrByName(const std::string &name);

  // Returns an amino acid based on the one letter representation. Returns
  // null if the one letter representation does not exist.
  static AminoAcidPtr getAminoAcidPtrByOneLetter(const std::string &one_letter);

  // Returns an amino acid based on the three letter representation. Returns
  // null if the three letter representation does not exist.
  static AminoAcidPtr getAminoAcidPtrByThreeLetter(const std::string &three_letter);

  // Checks if the list contains an amino acid with the specific name.
  static bool containsName(const std::string &name);

  // Checks if the list contains an amino acid with the specific one letter
  // representation.
  static bool containsOneLetter(const std::string &one_letter);

  // Checks if the list contains an amino acid with the specific three letter
  // representation.
  static bool containsThreeLetter(const std::string &three_letter);

  static AminoAcidPtr getAminoAcidPtrFromXml(xercesc::DOMElement * element);

 private:
  static AminoAcidPtrVec amino_acid_ptr_vec_;

  static std::unordered_map<std::string, AminoAcidPtr> amino_acid_one_letter_map_;

  static std::unordered_map<std::string, AminoAcidPtr> amino_acid_three_letter_map_;

  static std::unordered_map<std::string, AminoAcidPtr> amino_acid_name_map_;

  static AminoAcidPtr empty_amino_acid_ptr_;
};

}  // namespace toppic
#endif
