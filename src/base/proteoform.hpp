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


#ifndef PROT_BASE_PROTEOFORM_HPP_
#define PROT_BASE_PROTEOFORM_HPP_

#include "base/residue_freq.hpp"
#include "base/fasta_seq.hpp"
#include "base/fasta_index_reader.hpp"
#include "base/bp_spec.hpp"
#include "base/change.hpp"
#include "base/segment.hpp"
#include "base/prot_mod.hpp"
#include "base/align_type.hpp"

namespace prot {

class Proteoform;

typedef std::shared_ptr<Proteoform> ProteoformPtr;

class Proteoform {
 public:
  Proteoform(FastaSeqPtr fasta_seq_ptr, ProtModPtr prot_mod_ptr, 
             int start_pos, int end_pos, ResSeqPtr res_seq_ptr, 
             const ChangePtrVec &change_ptr_vec);

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

  int getChangeNum() {return change_list_.size();}

  ChangePtrVec getChangePtrVec() {return change_list_;}

  int getSpeciesId() {return species_id_;}

  void setSpeciesId(int id) {species_id_ = id;}

  int getProtId() {return prot_id_;}

  void setProtId(int id) {prot_id_ = id;}

  double getMass();

  AlignTypePtr getAlignType();

  int getChangeNum(ChangeTypePtr ct_ptr);

  ChangePtrVec getChangePtrVec(ChangeTypePtr ct_ptr);

  void addChangePtrVec(ChangePtrVec &change_ptr_vec);

  void addChangePtr(ChangePtr &change_ptr);

  void rmChangePtr(ChangePtr &change_ptr);

  SegmentPtrVec getSegmentPtrVec();

  std::string getProteinMatchSeq();

  std::string toString();

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  void parseXml(xercesc::DOMElement* element, ProteoformPtr db_proteoform,
                int adjust = 0);

  static std::string getXmlElementName() {return "proteoform";}

  void setVariablePtmNum(int n) {variable_ptm_num_ = n;}

  int getVariablePtmNum() {return variable_ptm_num_;}

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

  int species_id_ = 0;

  int prot_id_ = 0;

  ChangePtrVec change_list_;

  // Number of variable ptms is used for the test of the mass graph approach
  int variable_ptm_num_ = 0;
};

typedef std::vector<ProteoformPtr> ProteoformPtrVec;
typedef std::vector<ProteoformPtrVec> ProteoformPtrVec2D;

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
