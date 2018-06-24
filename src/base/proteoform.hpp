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


#ifndef PROT_BASE_PROTEOFORM_HPP_
#define PROT_BASE_PROTEOFORM_HPP_

#include <string>
#include <vector>

#include "base/residue_freq.hpp"
#include "base/fasta_seq.hpp"
#include "base/fasta_index_reader.hpp"
#include "base/bp_spec.hpp"
#include "base/mass_shift.hpp"
#include "base/segment.hpp"
#include "base/prot_mod.hpp"
#include "base/align_type.hpp"
#include "base/ptm.hpp"

namespace prot {

class Proteoform;

typedef std::shared_ptr<Proteoform> ProteoformPtr;

class Proteoform {
 public:
  Proteoform(FastaSeqPtr fasta_seq_ptr,
             ProtModPtr prot_mod_ptr,
             int start_pos, int end_pos,
             ResSeqPtr res_seq_ptr,
             const MassShiftPtrVec & mass_shift_ptr_vec);

  Proteoform(xercesc::DOMElement* element, FastaIndexReaderPtr reader_ptr,
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

  int getMassShiftNum(MassShiftTypePtr ct_ptr);

  MassShiftPtrVec getMassShiftPtrVec() {return mass_shift_list_;}

  MassShiftPtrVec getMassShiftPtrVec(MassShiftTypePtr ct_ptr);

  int getProteoClusterId() {return proteo_cluster_id_;}

  void setProteoClusterId(int id) {proteo_cluster_id_ = id;}

  int getProtId() {return prot_id_;}

  void setProtId(int id) {prot_id_ = id;}

  double getMass();

  AlignTypePtr getAlignType();

  void addMassShiftPtrVec(const MassShiftPtrVec & shift_ptr_vec);

  SegmentPtrVec getSegmentPtrVec();

  std::string getProteinMatchSeq();

  std::string toString();

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  void parseXml(xercesc::DOMElement* element, ProteoformPtr db_proteoform);

  static std::string getXmlElementName() {return "proteoform";}

  void setVariablePtmNum(int n) {variable_ptm_num_ = n;}

  int getVariablePtmNum() {return variable_ptm_num_;}

  std::string getMIScore();

  PtmPtrVec getPtmVec();

  PtmPtrVec getPtmVec(MassShiftTypePtr type);

 private:
  FastaSeqPtr fasta_seq_ptr_;

  ProtModPtr prot_mod_ptr_;

  // start and end positions are relative to the
  // database sequence
  int start_pos_;
  int end_pos_;

  /* residue_seq starts from start_pos_ and ends at end_pos_, and contains
   * fixed and variable modifications */
  ResSeqPtr residue_seq_ptr_;

  /* bp_spec is generated from residue_seq */
  BpSpecPtr bp_spec_ptr_;

  int proteo_cluster_id_ = -1;

  int prot_id_ = -1;

  MassShiftPtrVec mass_shift_list_;

  // Number of variable ptms is used for the test of the mass graph approach
  int variable_ptm_num_ = 0;

  std::string mi_score_ = "";
};

typedef std::vector<ProteoformPtr> ProteoformPtrVec;
typedef std::vector<ProteoformPtrVec> ProteoformPtrVec2D;

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
