#ifndef PROT_PROTEOFORM_HPP_
#define PROT_PROTEOFORM_HPP_

#include "htslib/faidx.h"

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
               const ChangePtrVec &change_ptr_vec);

    Proteoform(xercesc::DOMElement* element, const ProteoformPtrVec &db_proteoforms);

    Proteoform(xercesc::DOMElement* element, faidx_t *fai, const ResiduePtrVec &residue_ptr_vec);

    void parseXml(xercesc::DOMElement* element, ProteoformPtr db_proteoform);

    DbResSeqPtr getDbResSeqPtr() { return db_residue_seq_ptr_;}

    ProtModPtr getProtModPtr() { return prot_mod_ptr_;}

    ResSeqPtr getResSeqPtr() { return residue_seq_ptr_;}

    BpSpecPtr getBpSpecPtr() { return bp_spec_ptr_;}

    int getStartPos() { return start_pos_;}

    void setStartPos(int i);

    int getEndPos() { return end_pos_;}

    void setEndPos(int i);

    void setSplit(int i) { split_pos_ = i;}

    void setSplit() { split_pos_ = end_pos_;}

    int getLen() { return end_pos_ - start_pos_ + 1; }

    int getSeqId() { return db_residue_seq_ptr_->getId(); }

    const std::string& getSeqName() { return db_residue_seq_ptr_->getName(); }

    const std::string& getSeqDesc() { return db_residue_seq_ptr_->getDesc(); }

    ChangePtrVec getChangePtrVec() { return change_list_; }

    void rmChangePtr(int j) { change_list_.erase(change_list_.begin() + j); }

    void rmChangePtr(ChangePtr c);

    void addChangePtr(ChangePtr c) { change_list_.push_back(c); }

    int getSpeciesId() { return species_id_; }

    void setSpeciesId(int id) { species_id_ = id; }

    void setResSeqPtr(ResSeqPtr residue_seq_ptr);

    int getUnexpectedChangeNum();

    SemiAlignTypePtr getSemiAlignType();

    double getMass();

    ChangePtrVec getUnexpectedChangePtrVec() ;

    std::vector<int> getUnexpectedChangeId();

    SegmentPtrVec getSegmentPtrVec();

    std::string getProteinMatchSeq();

    std::string getProteinMatchSeq(double err);

    std::string toString();

    void addUnexpectedChangePtrVec(const ChangePtrVec &changes);

    void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

    bool isAcety();

  private:

    /* db_residue_seq contains fixed modifications */
    DbResSeqPtr db_residue_seq_ptr_;

    ProtModPtr prot_mod_ptr_;

    /* residue_seq starts from start_pos_ and ends at end_pos_, and contains
     * fixed and variable modifications */
    ResSeqPtr residue_seq_ptr_;

    /* bp_spec is generated from residue_seq */
    BpSpecPtr bp_spec_ptr_;

    // start and end positions are relative to the
    // database residue sequence
    int start_pos_;

    int end_pos_;

    int species_id_ = 0;

    ChangePtrVec change_list_;

    int split_pos_ = 0;
    // the splitting position for changes
    //   //if only one change, it's the same with end_pos_

};

/* get db proteoform */
ProteoformPtr getDbProteoformPtr(DbResSeqPtr db_res_seq_ptr);

/* generate a proteoform with protein mod */
ProteoformPtr getProtModProteoform(ProteoformPtr db_form_ptr,
                                   ProtModPtr prot_mod_ptr);

ProteoformPtrVec generateProtModProteoform(ProteoformPtr db_form_ptr,
        const ProtModPtrVec &prot_mod_ptrs);

ProteoformPtrVec2D generate2DProtModProteoform(const ProteoformPtrVec &db_form_ptrs,
        const ProtModPtrVec &prot_mod_ptrs);
/*
 * get subproteoform. local_start and local_end are relatively to
 * the start position in the original proteoform
 */
ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr,
                               int local_start, int local_end);

/* generate a proteoform vector with protein mod */
ProteoformPtrVec generateProtModProteoform(const ProteoformPtrVec &ori_forms,
        const ProtModPtrVec &prot_mods);

/* calculate frequencies for n_terminal_residues */
ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms);

/* calculater frequences for all residues */
ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                              const ProteoformPtrVec &raw_mods);

bool isSamePeptideAndMass(ProteoformPtr a, ProteoformPtr b, double ppo);

bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b, double ppo);

ProteoformPtrVec2D getProteoBlocks(const ProteoformPtrVec &proteo_ptrs, int db_block_size);

} /* namespace prot */

#endif /* PROTEOFORM_HPP_ */
