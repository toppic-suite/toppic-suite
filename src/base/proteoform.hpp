#ifndef PROT_PROTEOFORM_HPP_
#define PROT_PROTEOFORM_HPP_

#include "base/residue_freq.hpp"
#include "base/db_residue_seq.hpp"
#include "base/bp_spec.hpp"
#include "base/change.hpp"
#include "base/segment.hpp"
#include "base/prot_mod.hpp"
#include "base/semi_align_type.hpp"
#include "base/base_data.hpp"


namespace prot {

class Proteoform;
typedef std::shared_ptr<Proteoform> ProteoformPtr;
typedef std::vector<ProteoformPtr> ProteoformPtrVec;
typedef std::vector<ProteoformPtrVec> ProteoformPtrVec2D;

class Proteoform {
 public:
  Proteoform(DbResSeqPtr db_res_seq_ptr, ProtModPtr prot_mod_ptr,  
             ResSeqPtr res_seq_ptr, int start_pos, int end_pos, 
             ChangePtrVec change_list);

  Proteoform(xercesc::DOMElement* element,ProteoformPtrVec proteoforms);

  DbResSeqPtr getDbResSeqPtr() {return db_residue_seq_ptr_;}

  ProtModPtr getProtModPtr() {return prot_mod_ptr_;}

  ResSeqPtr getResSeqPtr() {return residue_seq_ptr_;}

  void setResSeqPtr(ResSeqPtr residue_seq_ptr){ residue_seq_ptr_ =residue_seq_ptr;}

  BpSpecPtr getBpSpecPtr() {return bp_spec_ptr_;}

  int getStartPos() {return start_pos_;}

  int getEndPos() {return end_pos_;}

  int getLen() {return end_pos_ - start_pos_ + 1;}

  int getSeqId() {return db_residue_seq_ptr_->getId();}

  std::string getName() {return db_residue_seq_ptr_->getName();}

  ChangePtrVec getChangePtrVec() {return change_list_;}

  void addUnexpectedChangePtrVec(ChangePtrVec &changes);

  SegmentPtrVec getSegmentPtrVec();

  std::string toString();

  int getUnexpectedChangeNum();

  ChangePtrVec getUnexpectedChangePtrVec() ;

  SemiAlignTypePtr getSemiAlignType();

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  int getSpeciesId(){return species_id_;}

  void setSpeciesId(int id){species_id_ = id;}

  double getMass();

  std::string getProteinMatchSeq();

 private:

  DbResSeqPtr db_residue_seq_ptr_;

  ProtModPtr prot_mod_ptr_;

  ResSeqPtr residue_seq_ptr_;

  BpSpecPtr bp_spec_ptr_;

  // start and end positions are relative to the 
  // database residue sequence
  int start_pos_;

  int end_pos_;

  int species_id_;

  ChangePtrVec change_list_;
};

/* get db proteoform */
ProteoformPtr getDbProteoformPtr(DbResSeqPtr db_res_seq_ptr, 
                                 ProtModPtr prot_mod_ptr);
/* generate a proteoform with protein mod */ 
ProteoformPtr getProtModProteoform(ProteoformPtr raw_form_ptr, 
                                   ProtModPtr prot_mod_ptr); 

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr, int start, int end);

/* generate a proteoform vector with protein mod */ 
ProteoformPtrVec generateProtModProteoform(ProteoformPtrVec &ori_forms,
                                           ProtModPtrVec &prot_mods);

/* calculate frequencies for n_terminal_residues */
ResFreqPtrVec compNTermResidueFreq(ProteoformPtrVec &prot_mod_forms);

/* calculater frequences for all residues */
ResFreqPtrVec compResidueFreq(ResiduePtrVec &residue_list, 
                              ProteoformPtrVec &raw_mods);

bool isSamePeptideAndMass(ProteoformPtr proteoform,ProteoformPtr another_proteoform,double ppo);
bool isStrictCompatiablePtmSpecies(ProteoformPtr a,ProteoformPtr b,double ppo);

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
