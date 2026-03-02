//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_COMMON_BASE_RESIDUE_HPP_
#define TOPPIC_COMMON_BASE_RESIDUE_HPP_

#include "common/base/amino_acid.hpp"
#include "common/base/ptm.hpp"

namespace toppic {

class Residue;
typedef std::shared_ptr<Residue> ResiduePtr;

class Residue {
 public:
  Residue(AminoAcidPtr acid_ptr, PtmPtr ptm_ptr);

  explicit Residue(XmlDOMElement* element);
  /** Get amino acid. */
  AminoAcidPtr getAminoAcidPtr() const {return acid_ptr_; }
  /** Get residue mass. */
  double getMass() const { return mass_; }
  /** Get post-translational modification. */
  PtmPtr getPtmPtr() const { return ptm_ptr_; }
  /** Checks if the residue contains the same amino acid and ptm.  */
  bool isSame(ResiduePtr residue_ptr) const;
  /** Get string representation */
  std::string toString(const std::string &delim_bgn, const std::string &delim_end) const;

  std::string toString() const {return toString("[", "]");}

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent,
                 const std::string &element_name) const;

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) const;

  static std::string getXmlElementName() {return "residue";}

 private:
  /** amino acid */
  AminoAcidPtr acid_ptr_;
  /** post-translational modification */
  PtmPtr ptm_ptr_;
  /** residue mass */
  double mass_;
};

typedef std::vector<ResiduePtr> ResiduePtrVec;
typedef std::vector<ResiduePtrVec> ResiduePtrVec2D;

}  // namespace toppic

#endif
