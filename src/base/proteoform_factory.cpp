#include <sstream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/proteoform_factory.hpp"

namespace prot {

ProteoformPtr ProteoformFactory::geneDbProteoformPtr(DbResSeqPtr db_res_seq_ptr) {
  int start_pos = 0;
  int end_pos = db_res_seq_ptr->getLen() - 1;
  ChangePtrVec change_list;
  for (int i = 0; i < db_res_seq_ptr->getLen(); i++) {
    PtmPtr ptm_ptr = db_res_seq_ptr->getResiduePtr(i)->getPtmPtr();
    if (ptm_ptr == PtmBase::getEmptyPtmPtr()) {
      ChangePtr change_ptr = ChangePtr(
          new Change(i, i+1, ChangeType::FIXED, ptm_ptr->getMonoMass(), ptm_ptr));
      change_list.push_back(change_ptr);
    }
  }
  ProtModPtr none_prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  return ProteoformPtr(new Proteoform(db_res_seq_ptr, none_prot_mod_ptr,
                                      db_res_seq_ptr, start_pos, end_pos,
                                      change_list));
}

ProteoformPtr getProtModProteoform(ProteoformPtr db_form_ptr,
                                   ProtModPtr prot_mod_ptr) {
    // check if the proteoform can be truncated
    TruncPtr trunc_ptr = prot_mod_ptr->getTruncPtr();
    DbResSeqPtr db_res_seq_ptr = db_form_ptr->getDbResSeqPtr();
    bool valid_trunc = trunc_ptr->isValidTrunc(db_res_seq_ptr->getResidues());
    if (!valid_trunc) {
        //LOG_DEBUG("NO valid trunc");
        return ProteoformPtr(nullptr);
    }

    int start_res = trunc_ptr->getTruncLen();
    /* last bp index */
    int end_res = db_form_ptr->getLen();
    ChangePtrVec ori_change_ptrs = db_form_ptr->getChangePtrVec();
    ChangePtrVec change_ptrs;
    for (size_t i = 0; i < ori_change_ptrs.size(); i++) {
        if (ori_change_ptrs[i]->getLeftBpPos() >= start_res
                && ori_change_ptrs[i]->getRightBpPos() <= end_res + 1) {
            ChangePtr change_ptr = ChangePtr(new Change(ori_change_ptrs[i], start_res));
            change_ptrs.push_back(change_ptr);
        }
    }

    // first residue might be acetylated
    ResiduePtrVec residue_ptrs;
    ResiduePtr first_residue_ptr = db_res_seq_ptr->getResiduePtr(start_res);

    PtmPtr ori_ptm_ptr = first_residue_ptr->getPtmPtr();
    PtmPtr prot_ptm_ptr = prot_mod_ptr->getPtmPtr();
    /* if there is a conflict */
    if (!ori_ptm_ptr ->isEmpty() && !prot_ptm_ptr->isEmpty()) {
        return ProteoformPtr(nullptr);
    }

    if (prot_ptm_ptr->isEmpty()) {
        residue_ptrs.push_back(first_residue_ptr);
    } else {
        /* add protein n-terminal mod */
        AcidPtr acid_ptr = first_residue_ptr->getAcidPtr();
        ResiduePtr mut_residue_ptr = ResidueFactory::getBaseResiduePtrByAcidPtm(acid_ptr, prot_ptm_ptr);
        if (mut_residue_ptr == nullptr) {
            LOG_ERROR( "Proteoform:: residue not found");
            throw("Residue not found");
        }
        residue_ptrs.push_back(mut_residue_ptr);
        change_ptrs.push_back(ChangePtr(new Change(0,1, PROTEIN_VARIABLE_CHANGE,
                                        prot_ptm_ptr->getMonoMass(), prot_ptm_ptr)));
    }

    // add all other residues
    for (int i = start_res + 1; i < db_res_seq_ptr->getLen(); i++) {
        residue_ptrs.push_back(db_res_seq_ptr->getResiduePtr(i));
    }
    ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residue_ptrs));

    //LOG_DEBUG("mod protein sequence name " << db_res_seq_ptr->getName()
    //<< " len " << db_res_seq_ptr->getLen());
    return ProteoformPtr(
               new Proteoform(db_res_seq_ptr, prot_mod_ptr, seq_ptr, start_res,
                              db_res_seq_ptr->getLen()-1, change_ptrs));
}

ProteoformPtr getSubProteoform(ProteoformPtr proteoform_ptr,
                               int local_start, int local_end) {
    ResiduePtrVec residues;
    ResSeqPtr res_seq_ptr = proteoform_ptr->getResSeqPtr();
    for (int i = local_start; i <= local_end; i++) {
        residues.push_back(res_seq_ptr->getResiduePtr(i));
    }
    ResSeqPtr seq_ptr = ResSeqPtr(new ResidueSeq(residues));
    ChangePtrVec change_list;
    ChangePtrVec ori_change_list = proteoform_ptr->getChangePtrVec();
    for (size_t i = 0; i < ori_change_list.size(); i++) {
        if (ori_change_list[i]->getLeftBpPos() >= local_start
                && ori_change_list[i]->getRightBpPos() <= local_end + 1) {
            ChangePtr change_ptr = ChangePtr(new Change(ori_change_list[i], local_start));
            change_list.push_back(change_ptr);
        }
    }
    DbResSeqPtr db_res_seq_ptr = proteoform_ptr->getDbResSeqPtr();
    ProtModPtr prot_mod_ptr = proteoform_ptr->getProtModPtr();
    return ProteoformPtr(
               new Proteoform(db_res_seq_ptr, prot_mod_ptr, seq_ptr,
                              local_start + proteoform_ptr->getStartPos(),
                              local_end + proteoform_ptr->getStartPos(), change_list));
}

ProteoformPtrVec generateProtModProteoform(ProteoformPtr proteo_ptr,
        const ProtModPtrVec &prot_mods) {
    ProteoformPtrVec new_forms;
    for (size_t j = 0; j < prot_mods.size(); j++) {
        ProteoformPtr ptr = getProtModProteoform(proteo_ptr, prot_mods[j]);
        if (ptr.get() != nullptr) {
            new_forms.push_back(ptr);
        }
    }
    return new_forms;
}

ProteoformPtrVec generateProtModProteoform(const ProteoformPtrVec &ori_forms,
        const ProtModPtrVec &prot_mods) {
    ProteoformPtrVec new_forms;
    for (size_t i = 0; i < ori_forms.size(); i++) {
        for (size_t j = 0; j < prot_mods.size(); j++) {
            ProteoformPtr ptr = getProtModProteoform(ori_forms[i], prot_mods[j]);
            if (ptr.get() != nullptr) {
                new_forms.push_back(ptr);
            }
        }
    }
    return new_forms;
}

ProteoformPtrVec2D generate2DProtModProteoform(const ProteoformPtrVec &ori_forms,
        const ProtModPtrVec &prot_mods) {
    ProteoformPtrVec2D new_forms;
    for (size_t i = 0; i < ori_forms.size(); i++) {
        ProteoformPtrVec mod_forms;
        for (size_t j = 0; j < prot_mods.size(); j++) {
            ProteoformPtr ptr = getProtModProteoform(ori_forms[i], prot_mods[j]);
            if (ptr.get() != nullptr) {
                mod_forms.push_back(ptr);
            }
        }
        new_forms.push_back(mod_forms);
    }
    return new_forms;
}

ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms) {
    std::vector<double> counts;
    ResiduePtrVec residue_list;
    for (size_t i = 0; i < prot_mod_forms.size(); i++) {
        ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
        if (seq_ptr->getLen() >= 1) {
            ResiduePtr res_ptr = seq_ptr->getResiduePtr(0);
            int pos = findResidue(residue_list, res_ptr);
            if (pos >= 0) {
                // found
                counts[pos] = counts[pos]+1;
            } else {
                residue_list.push_back(res_ptr);
                counts.push_back(1);
            }
        }
    }

    double sum = 0;
    for (size_t i = 0; i < counts.size(); i++) {
        sum = sum + counts[i];
    }
    ResFreqPtrVec res_freq_list;
    for (size_t i = 0; i < residue_list.size(); i++) {
        ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(),
                                                residue_list[i]->getPtmPtr(),
                                                counts[i]/sum));
        res_freq_list.push_back(res_freq_ptr);
    }
    return res_freq_list;
}


ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                              const ProteoformPtrVec &prot_mod_forms) {
    std::vector<double> counts(residue_list.size(), 0.0);
    for (size_t i = 0; i < prot_mod_forms.size(); i++) {
        ResSeqPtr seq_ptr = prot_mod_forms[i]->getResSeqPtr();
        for (int j = 0; j < seq_ptr->getLen(); j++) {
            ResiduePtr res_ptr = seq_ptr->getResiduePtr(j);
            int pos = findResidue(residue_list, res_ptr);
            if (pos >= 0) {
                // found
                counts[pos] = counts[pos]+1;
            }
        }
    }

    double sum = 0;
    for (size_t i = 0; i < counts.size(); i++) {
        sum = sum + counts[i];
    }
    ResFreqPtrVec res_freq_list;
    for (size_t i = 0; i < residue_list.size(); i++) {
        ResFreqPtr res_freq_ptr(new ResidueFreq(residue_list[i]->getAcidPtr(),
                                                residue_list[i]->getPtmPtr(),
                                                counts[i]/sum));
        res_freq_list.push_back(res_freq_ptr);
    }
    return res_freq_list;
}

void Proteoform::addUnexpectedChangePtrVec(const ChangePtrVec &changes) {
    for (size_t i = 0; i < changes.size(); i++) {
        if (changes[i]->getChangeType() == UNEXPECTED_CHANGE) {
            change_list_.push_back(changes[i]);
        }
    }
}

bool isSamePeptideAndMass(ProteoformPtr a, ProteoformPtr b, double ppo) {
    if(a->getDbResSeqPtr()->getId() != b->getDbResSeqPtr()->getId()) {
        return false;
    }
    if(a->getStartPos() != b->getStartPos()) {
        return false;
    }
    if(a->getEndPos() != b->getEndPos()) {
        return false;
    }
    double thresh = a->getResSeqPtr()->getSeqMass() * ppo;
    if(std::abs(a->getResSeqPtr()->getSeqMass()
                -b->getResSeqPtr()->getSeqMass())> thresh) {
        return false;
    }
    return true;
}

bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b,
                                   double ppo) {
    if(!isSamePeptideAndMass(a,b,ppo)) {
        return false;
    }
    if(a->getChangePtrVec().size() != b->getChangePtrVec().size()) {
        return false;
    }
    double shift_tolerance = a->getResSeqPtr()->getSeqMass()*ppo;
    // sort changes
    ChangePtrVec a_change_vec = a->getChangePtrVec();
    ChangePtrVec b_change_vec = b->getChangePtrVec();
    std::sort(a_change_vec.begin(),a_change_vec.end(),compareChangeUp);
    std::sort(b_change_vec.begin(),b_change_vec.end(),compareChangeUp);
    for(size_t i=0; i< a->getChangePtrVec().size(); i++) {
        ChangePtr ac = a_change_vec[i];
        ChangePtr bc = b_change_vec[i];
        if(ac->getRightBpPos() <= bc->getLeftBpPos() || bc->getRightBpPos() <= ac->getLeftBpPos()) {
            return false;
        }
        if(std::abs(ac->getMassShift()-bc->getMassShift()) > shift_tolerance) {
            return false;
        }
    }
    return true;
}

ProteoformPtrVec2D getProteoBlocks(const ProteoformPtrVec &proteo_ptrs, int db_block_size) {
    size_t start_idx = 0;
    size_t proteo_idx =0;
    int block_len=0;
    ProteoformPtrVec2D proteo_blocks;
    while(proteo_idx < proteo_ptrs.size()) {
        int proteo_len = proteo_ptrs[proteo_idx]->getResSeqPtr()->getLen();
        if(block_len + proteo_len < db_block_size) {
            block_len = block_len + proteo_len;
            proteo_idx++;
        } else {
            size_t end_idx = proteo_idx;
            ProteoformPtrVec proteo_in_block;
            for(size_t i = start_idx; i<=end_idx; i++) {
                proteo_in_block.push_back(proteo_ptrs[i]);
            }
            proteo_blocks.push_back(proteo_in_block);
            start_idx = end_idx +1;
            proteo_idx = end_idx +1;
            block_len = 0;
        }
    }
    /* last block */
    if (start_idx < proteo_ptrs.size()) {
        ProteoformPtrVec proteo_in_block;
        for(size_t i = start_idx; i < proteo_ptrs.size(); i++) {
            proteo_in_block.push_back(proteo_ptrs[i]);
        }
        proteo_blocks.push_back(proteo_in_block);
    }
    return proteo_blocks;
}

void Proteoform::rmChangePtr(ChangePtr c) {
    change_list_.erase(std::remove(change_list_.begin(), change_list_.end(), c),
                       change_list_.end());
}

void Proteoform::addChangePtrVec(const ChangePtrVec& changes) {
    for (size_t i = 0; i < changes.size(); i++) {
        change_list_.push_back(changes[i]);
    }
}


} /* namespace prot */

