#ifndef PROT_PROTEOFORM_HPP_
#define PROT_PROTEOFORM_HPP_

#include "base/residue_seq.hpp"
#include "base/bp_spec.hpp"
#include "base/change.hpp"
#include "base/segment.hpp"
#include "base/prot_mod.hpp"

namespace prot {

class Proteoform;
typedef std::shared_ptr<Proteoform> ProteoformPtr;

class Proteoform {
public:
	Proteoform(ProteoformPtr ori_form_ptr, ProtModPtr prot_mod_ptr, std::string name, 
             ResSeqPtr res_seq_ptr, int start_pos, int end_pos, ChangePtrVec change_list);

  ProteoformPtr getRawFormPtr() {return raw_form_ptr_;}

  ProtModPtr getProtModPtr() {return prot_mod_ptr_;}

  std::string getName() {return name_;}

	ResSeqPtr getResSeqPtr() {return residue_seq_ptr_;}

	BpSpecPtr getBpSpecPtr() {return bp_spec_ptr_;}

  int getStartPos() {return start_pos_;}

  int getEndPos() {return end_pos_;}

  int getLen() {return end_pos_ - start_pos_ + 1;}

  ChangePtrVec getChangePtrVec() {return change_list_;}

  SegmentPtrVec getSegmentPtrVec();

  std::string toString();


private:
  
  ProteoformPtr raw_form_ptr_;

  ProtModPtr prot_mod_ptr_;
  
  std::string name_;

	ResSeqPtr residue_seq_ptr_;

  BpSpecPtr bp_spec_ptr_;

  int start_pos_;

  int end_pos_;

  ChangePtrVec change_list_;
};

typedef std::vector<ProteoformPtr> ProteoformPtrVec;

ProteoformPtr getRawProteoformPtr(std::string name, ResSeqPtr res_seq_ptr);

ProteoformPtr getProtModProteoform(ProteoformPtr raw_form_ptr, 
                                   ResiduePtrVec &residue_list, ProtModPtr prot_mod_ptr); 

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr, int start, int end);

ProteoformPtrVec generateProtModProteoform(ProteoformPtrVec &ori_forms,
                                           ResiduePtrVec &residue_list,
                                           ProtModPtrVec &prot_mods);

ProteoformPtrVec getProtModNoneProteoform(ProteoformPtrVec &all_forms); 


} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
