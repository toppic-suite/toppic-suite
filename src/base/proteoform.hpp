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

   double getMass();

   AlignTypePtr getAlignType();

   int getChangeNum(ChangeTypePtr ct_ptr);

   ChangePtrVec getChangePtrVec(ChangeTypePtr ct_ptr);

   SegmentPtrVec getSegmentPtrVec();

   std::string getProteinMatchSeq();

   std::string toString();

   void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

   void parseXml(xercesc::DOMElement* element, ProteoformPtr db_proteoform);

   static std::string getXmlElementName() {return "proteoform";}

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

   ChangePtrVec change_list_;
};

typedef std::vector<ProteoformPtr> ProteoformPtrVec;
typedef std::vector<ProteoformPtrVec> ProteoformPtrVec2D;

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
