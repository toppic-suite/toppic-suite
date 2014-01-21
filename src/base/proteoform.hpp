#ifndef PROT_PROTEOFORM_HPP_
#define PROT_PROTEOFORM_HPP_

#include "base/db_residue_seq.hpp"
#include "base/bp_spec.hpp"
#include "base/change.hpp"
#include "base/segment.hpp"
#include "base/prot_mod.hpp"


namespace prot {

#define SEMI_ALIGN_TYPE_COMPLETE 0
#define SEMI_ALIGN_TYPE_PREFIX   1
#define SEMI_ALIGN_TYPE_SUFFIX   2
#define SEMI_ALIGN_TYPE_INTERNAL 3


class Proteoform;
typedef std::shared_ptr<Proteoform> ProteoformPtr;

class Proteoform {
public:
	Proteoform(DbResSeqPtr db_res_seq_ptr, ProtModPtr prot_mod_ptr,  
             ResSeqPtr res_seq_ptr, int start_pos, int end_pos, 
             ChangePtrVec change_list);

  DbResSeqPtr getDbResSeqPtr() {return db_residue_seq_ptr_;}

  ProtModPtr getProtModPtr() {return prot_mod_ptr_;}

	ResSeqPtr getResSeqPtr() {return residue_seq_ptr_;}

	BpSpecPtr getBpSpecPtr() {
		std::cout<<bp_spec_ptr_<<std::endl;
		return bp_spec_ptr_;
	}

  int getStartPos() {return start_pos_;}

  int getEndPos() {return end_pos_;}

  int getLen() {return end_pos_ - start_pos_ + 1;}

  int getSeqId() {return db_residue_seq_ptr_->getId();}

  std::string getName() {return db_residue_seq_ptr_->getName();}

  ChangePtrVec getChangePtrVec() {return change_list_;}

  SegmentPtrVec getSegmentPtrVec();

  std::string toString();

  int getUnexpectedChangeNum();
  
  int getSemiAlignType();

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

private:
  
  DbResSeqPtr db_residue_seq_ptr_;

  ProtModPtr prot_mod_ptr_;
  
  ResSeqPtr residue_seq_ptr_;

  BpSpecPtr bp_spec_ptr_;

  int start_pos_;

  int end_pos_;

  ChangePtrVec change_list_;
};

typedef std::vector<ProteoformPtr> ProteoformPtrVec;

ProteoformPtr getDbProteoformPtr(DbResSeqPtr db_res_seq_ptr, 
                                 ProtModPtr prot_mod_ptr);

ProteoformPtr getProtModProteoform(ProteoformPtr raw_form_ptr, 
                                   ResiduePtrVec &residue_list, 
                                   ProtModPtr prot_mod_ptr); 

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr, int start, int end);

ProteoformPtrVec generateProtModProteoform(ProteoformPtrVec &ori_forms,
                                           ResiduePtrVec &residue_list,
                                           ProtModPtrVec &prot_mods);
std::string convertSemiAlignmentTypeToString(int i);

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
