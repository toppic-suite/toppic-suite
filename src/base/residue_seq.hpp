// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_BASE_RESIDUE_SEQ_HPP_
#define PROT_BASE_RESIDUE_SEQ_HPP_

#include <string>

#include "base/mass_constant.hpp"
#include "base/residue.hpp"

namespace prot {

class ResidueSeq;

typedef std::shared_ptr<ResidueSeq> ResSeqPtr;
typedef std::vector<ResSeqPtr> ResSeqPtrVec;

class ResidueSeq {
 public:
  ResidueSeq(const ResiduePtrVec &residues);

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
    return residue_mass_sum_ + MassConstant::getWaterMass();
  }

  /** Gets the sum of residue masses */
  double getResMassSum() {return residue_mass_sum_;}

  std::string toString();

  std::string toAcidString();

  static std::string getXmlElementName() {return "residue_seq";}

  //void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  /** residue list */
  ResiduePtrVec residues_;
  /** the sum of residue mass */
  double residue_mass_sum_;

  static ResSeqPtr getEmptyResidueSeq();
};


}
#endif
