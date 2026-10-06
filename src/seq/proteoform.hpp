// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_SEQ_PROTEOFORM_HPP_
#define TOPPIC_SEQ_PROTEOFORM_HPP_

#include <memory>
#include <string>
#include <vector>

#include "common/base/mod.hpp"
#include "common/base/prot_mod.hpp"
#include "common/base/proteoform_type.hpp"
#include "common/base/ptm.hpp"
#include "common/base/residue_freq.hpp"
#include "seq/bp_spec.hpp"
#include "seq/fasta_index_reader.hpp"
#include "seq/fasta_seq.hpp"
#include "seq/mass_shift.hpp"
#include "seq/seq_segment.hpp"

namespace toppic {

class Proteoform;

using ProteoformPtr = std::shared_ptr<Proteoform>;

class Proteoform {
 public:
  Proteoform(const FastaSeqPtr& fasta_seq_ptr, const ProtModPtr& prot_mod_ptr,
             int start_pos, int end_pos, const ResSeqPtr& res_seq_ptr,
             const MassShiftPtrVec& mass_shift_ptr_vec);

  Proteoform(XmlDOMElement element, const FastaIndexReaderPtr& reader_ptr,
             const ModPtrVec& fix_mod_list);

  FastaSeqPtr getFastaSeqPtr() const { return fasta_seq_ptr_; }

  const std::string& getSeqName() const { return fasta_seq_ptr_->getName(); }

  const std::string& getSeqDesc() const { return fasta_seq_ptr_->getDesc(); }

  int getStartPos() const { return start_pos_; }

  int getEndPos() const { return end_pos_; }

  ProtModPtr getProtModPtr() const { return prot_mod_ptr_; }

  // The variable PTM of the proteoform N-terminus (the alpha-amine of its
  // first residue), for a truncated proteoform; ModBase::getNTermNoneModPtr()
  // without one.
  const ModPtr& getTruncNTermModPtr() const { return trunc_n_term_mod_ptr_; }

  bool hasTruncNTermMod() const;

  // Sets the truncation N-terminal modification, which also becomes the N-
  // terminal modification of the residue sequence (its mass is in the
  // sequence mass and the prefix masses). Only a proteoform without a
  // protein N-terminal form (NONE) that starts at the third residue of the
  // protein or later can have one, and only on its first residue's amino
  // acid; otherwise nothing is changed and false is returned. The N-terminal
  // none modification removes it.
  bool setTruncNTermModPtr(const ModPtr& mod_ptr);

  ResSeqPtr getResSeqPtr() const { return residue_seq_ptr_; }

  BpSpecPtr getBpSpecPtr() const { return bp_spec_ptr_; }

  int getLen() const { return end_pos_ - start_pos_ + 1; }

  int getMassShiftNum() const {
    return static_cast<int>(mass_shift_list_.size());
  }

  int getAlterNum(const AlterTypePtr& type_ptr);

  // the number of variable PTMs: the variable mass shifts' and the
  // truncation N-terminal modification
  int getVarPtmNum();

  const MassShiftPtrVec& getMassShiftPtrVec() const { return mass_shift_list_; }

  MassShiftPtrVec getMassShiftPtrVec(const AlterTypePtr& type_ptr);

  int getProteoClusterId() const { return proteo_cluster_id_; }

  void setProteoClusterId(int id) { proteo_cluster_id_ = id; }

  int getProtId() const { return prot_id_; }

  void setProtId(int id) { prot_id_ = id; }

  // protein cluster (see prsm/prsm_prot_cluster); -1 until assigned
  int getProtClusterId() const { return prot_cluster_id_; }

  void setProtClusterId(int id) { prot_cluster_id_ = id; }

  double getProteoInte() const { return proteo_inte_; }

  void setProteoInte(double inte) { proteo_inte_ = inte; }

  double getMass();

  double getMinusWaterMass();

  ProteoformTypePtr getProteoformType();

  void addMassShiftPtrVec(const MassShiftPtrVec& shift_ptr_vec);

  SeqSegmentPtrVec getSeqSegmentPtrVec();

  std::string getProteoformMatchSeq();

  std::string getPrevAminoAcid();

  std::string getNextAminoAcid();

  std::string getAlterStr(const AlterTypePtr& type_ptr);

  // the variable PTMs as name:[position] (the Variable PTMs column): the
  // truncation N-terminal modification first, as name:[N-term], then the
  // variable mass shifts
  std::string getVarPtmStr();

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement parent);

  void parseXml(XmlDOMElement element, const ProteoformPtr& db_proteoform);

  static std::string getXmlElementName() { return "proteoform"; }

  PtmPtrVec getPtmVec(const AlterTypePtr& type);

  std::string getMIScore();

  // The level of the proteoform in the five-level classification of
  // proteoform identifications (Smith et al., Nat. Methods 16:939-940, 2019),
  // from its unexpected mass shifts, prot_hit_num, the number of database
  // proteins containing its sequence (the # Protein hits column), and its
  // variable PTM sites (hasUndecidedVarPtmSite with var_mod_ptrs):
  //                                           one protein   several proteins
  //   no unexpected shift, sites decided      "1"           "2D"
  //   no unexpected shift, a site undecided   "2A"          "3"
  //   unexpected shift(s)                     "4"           "5"
  // Protein N-terminal forms and fixed PTMs leave the level as it is; 0 hits
  // count as one protein.
  std::string getClassification(int prot_hit_num,
                                const ModPtrVec& var_mod_ptrs);

  // Whether the site of a variable PTM is undecided: a variable mass shift
  // covers more than one residue, and for one of its PTMs the residues of
  // its range that a modification of the PTM in var_mod_ptrs modifies are
  // not as many as the shift's alterations of the PTM, or one residue could
  // carry two of its PTMs. A shift on one residue is decided. toppic writes
  // a shift whose sites are decided one PTM per residue
  // (VarPtmAlign::splitDeterminedShift); TopMG's shifts span the gaps
  // between matched fragments, which may hold one candidate residue.
  bool hasUndecidedVarPtmSite(const ModPtrVec& var_mod_ptrs);

  // The level stored by setClassification, "" until set: the report stage
  // sets it from the protein hits (prsm_cutoff_selector), and the XML has it
  // as <classification>.
  const std::string& getClassification() const { return classification_; }

  void setClassification(const std::string& classification) {
    classification_ = classification;
  }

  void setStartPos(int start_pos) { start_pos_ = start_pos; }

  void setEndPos(int end_pos) { end_pos_ = end_pos; }

  void setFastaSeqPtr(const FastaSeqPtr& fasta_seq_ptr) {
    fasta_seq_ptr_ = fasta_seq_ptr;
  }

 private:
  FastaSeqPtr fasta_seq_ptr_;

  ProtModPtr prot_mod_ptr_;

  // the variable PTM of the N-terminus of a truncated proteoform (see
  // setTruncNTermModPtr); ModBase::getNTermNoneModPtr() without one
  ModPtr trunc_n_term_mod_ptr_;

  // start and end positions are relative to the
  // database sequence
  int start_pos_;
  int end_pos_;

  // residue_seq starts from start_pos_ and ends at end_pos_, and contains
  // fixed and variable modifications
  ResSeqPtr residue_seq_ptr_;

  // bp_spec is generated from residue_seq
  BpSpecPtr bp_spec_ptr_;

  int proteo_cluster_id_ = -1;

  int prot_id_ = -1;
  int prot_cluster_id_ = -1;

  double proteo_inte_ = -1;

  // the stored level of the proteoform classification, "" until set
  std::string classification_;

  MassShiftPtrVec mass_shift_list_;
};

using ProteoformPtrVec = std::vector<ProteoformPtr>;
using ProteoformPtrVec2D = std::vector<ProteoformPtrVec>;

}  // namespace toppic

#endif
