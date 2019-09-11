//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_SEQ_RESIDUE_SEQ_HPP_
#define TOPPIC_SEQ_RESIDUE_SEQ_HPP_

#include "common/base/mass_constant.hpp"
#include "common/base/residue.hpp"

namespace toppic {

class ResidueSeq;

typedef std::shared_ptr<ResidueSeq> ResSeqPtr;
typedef std::vector<ResSeqPtr> ResSeqPtrVec;

class ResidueSeq {
 public:
  explicit ResidueSeq(const ResiduePtrVec &residues);

  /**
   * Returns a sub-peptide of the original peptide.
   **/
  ResSeqPtr getSubResidueSeq(int bgn, int end);

  /** Gets length */
  int getLen() {return residues_.size();}

  /** Gets residue at position i */
  ResiduePtr getResiduePtr(int i) {return residues_[i];}

  /** Gets all residues */
  const ResiduePtrVec& getResidues() {return residues_;}

  /** Gets sequence molecular mass */
  double getSeqMass() {
    return residue_mass_sum_ + mass_constant::getWaterMass();
  }

  /** Gets the sum of residue masses */
  double getResMassSum() {return residue_mass_sum_;}

  std::string toString();

  std::string toAcidString();

  static std::string getXmlElementName() {return "residue_seq";}

  //void appendXml(XmlDOMDocument* xml_doc,XmlDOMElement* parent);

 private:
  /** residue list */
  ResiduePtrVec residues_;
  /** the sum of residue mass */
  double residue_mass_sum_;

  static ResSeqPtr getEmptyResidueSeq();
};

}  // namespace toppic
#endif
