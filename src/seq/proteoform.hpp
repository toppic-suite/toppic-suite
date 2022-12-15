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

#ifndef TOPPIC_SEQ_PROTEOFORM_HPP_
#define TOPPIC_SEQ_PROTEOFORM_HPP_

#include "common/base/ptm.hpp"
#include "common/base/prot_mod.hpp"
#include "common/base/residue_freq.hpp"
#include "seq/fasta_seq.hpp"
#include "seq/fasta_index_reader.hpp"
#include "seq/bp_spec.hpp"
#include "seq/mass_shift.hpp"
#include "seq/seq_segment.hpp"
#include "seq/proteoform_type.hpp"

namespace toppic {

class Proteoform;

typedef std::shared_ptr<Proteoform> ProteoformPtr;

class Proteoform {
 public:
  Proteoform(FastaSeqPtr fasta_seq_ptr,
             ProtModPtr prot_mod_ptr,
             int start_pos, int end_pos,
             ResSeqPtr res_seq_ptr,
             const MassShiftPtrVec &mass_shift_ptr_vec);

  Proteoform(XmlDOMElement* element, FastaIndexReaderPtr reader_ptr,
             const ModPtrVec &fix_mod_list);

  FastaSeqPtr getFastaSeqPtr() {return fasta_seq_ptr_;}

  std::string getSeqName() { return fasta_seq_ptr_->getName();}

  std::string getSeqDesc() { return fasta_seq_ptr_->getDesc();}

  int getStartPos() { return start_pos_;}

  int getEndPos() { return end_pos_;}

  ProtModPtr getProtModPtr() { return prot_mod_ptr_;}

  ResSeqPtr getResSeqPtr() { return residue_seq_ptr_;}

  BpSpecPtr getBpSpecPtr() { return bp_spec_ptr_;}

  int getLen() { return end_pos_ - start_pos_ + 1; }

  int getMassShiftNum() {return static_cast<int>(mass_shift_list_.size());}

  int getAlterNum(AlterTypePtr type_ptr);

  int getVarPtmNum(); 

  MassShiftPtrVec getMassShiftPtrVec() {return mass_shift_list_;}

  MassShiftPtrVec getMassShiftPtrVec(AlterTypePtr type_ptr);

  int getProteoClusterId() {return proteo_cluster_id_;}

  void setProteoClusterId(int id) {proteo_cluster_id_ = id;}

  int getProtId() {return prot_id_;}

  void setProtId(int id) {prot_id_ = id;}

  double getMass();

  double getMinusWaterMass();

  ProteoformTypePtr getProteoformType();

  void addMassShiftPtrVec(const MassShiftPtrVec & shift_ptr_vec);

  SeqSegmentPtrVec getSeqSegmentPtrVec();

  std::string getProteoformMatchSeq();

  std::string getAlterStr(AlterTypePtr type_ptr);

  void appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent);

  void parseXml(XmlDOMElement* element, ProteoformPtr db_proteoform);

  static std::string getXmlElementName() {return "proteoform";}

  PtmPtrVec getPtmVec(AlterTypePtr type);

  std::string getMIScore();

  void setStartPos(int start_pos) {start_pos_ = start_pos;}

  void setEndPos(int end_pos) {end_pos_ = end_pos;}

  void setFastaSeqPtr(FastaSeqPtr fasta_seq_ptr) {fasta_seq_ptr_ = fasta_seq_ptr;}

 private:
  FastaSeqPtr fasta_seq_ptr_;

  ProtModPtr prot_mod_ptr_;

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

  MassShiftPtrVec mass_shift_list_;
};

typedef std::vector<ProteoformPtr> ProteoformPtrVec;
typedef std::vector<ProteoformPtrVec> ProteoformPtrVec2D;

} /* namespace toppic */

#endif /* PROTEOFORM_HPP_ */
