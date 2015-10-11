#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "graph/proteo_anno.hpp"

namespace prot {

ProteoAnno::ProteoAnno(const ResiduePtrVec &fix_mod_res_ptr_vec,
                       const ProtModPtrVec &prot_mod_ptr_vec, 
                       const ResiduePtrVec &residue_mod_ptr_vec) {
  fix_mod_res_ptr_vec_ = fix_mod_res_ptr_vec;
  prot_mod_ptr_vec_ = prot_mod_ptr_vec;
  residue_mod_ptr_vec_ = residue_mod_ptr_vec;
  is_nme_ = false;
  for (size_t i = 0; i < residue_mod_ptr_vec_.size(); i++) {
    ResiduePtr res_ptr = residue_mod_ptr_vec[i];
    AcidPtr acid_ptr = res_ptr->getAcidPtr();
    if (ptm_map_.find(acid_ptr) == ptm_map_.end()) {
      ResiduePtrVec cur_vec;
      cur_vec.push_back(res_ptr);
      ptm_map_[acid_ptr]= cur_vec;
    }
    else {
      ptm_map_[acid_ptr].push_back(res_ptr);
    }
  }
}

void ProteoAnno::anno(const std::string &seq) {
  AcidPtrVec acid_ptr_vec = AcidFactory::convertSeqToAcidSeq(seq);
  ResiduePtrVec residue_ptr_vec = convertAcidToResidueSeq(fix_mod_res_ptr_vec_, acid_ptr_vec);
  res_vec_2d_.clear();
  change_vec_2d_.clear();
  // input and fixed mod
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    ResiduePtrVec cur_res_vec;
    std::vector<int> cur_change_vec;
    ResiduePtr res_ptr = residue_ptr_vec[i];
    cur_res_vec.push_back(res_ptr);
    if (res_ptr->getPtmPtr()->isEmpty()) {
      cur_change_vec.push_back(Change::getInputChange());
    }
    else {
      cur_change_vec.push_back(Change::getFixedChange());
    }
    res_vec_2d_.push_back(cur_res_vec);
    change_vec_2d_.push_back(cur_change_vec);
  }
  LOG_DEBUG("input complete");

  // protein mod
  for (size_t i = 0; i < prot_mod_ptr_vec_.size(); i++) {
    ProtModPtr mod_ptr = prot_mod_ptr_vec_[i];
    LOG_DEBUG("i " << i << " mod " << mod_ptr);
    if (mod_ptr == ProtModFactory::getProtModPtr_NME()) {
      LOG_DEBUG("NME");
      TruncPtr trunc_ptr = mod_ptr->getTruncPtr();
      bool valid_trunc = trunc_ptr->isValidTrunc(residue_ptr_vec);
      if (valid_trunc) {
        // add empty residue to the first methinine residue
        is_nme_ = true;
        AcidPtr acid_ptr = AcidFactory::findEmptyAcidPtr();
        PtmPtr ptm_ptr = PtmFactory::findEmptyPtmPtr();
        ResiduePtr empty_residue_ptr = ResidueFactory::addBaseResidue(acid_ptr, ptm_ptr);
        LOG_DEBUG("empty acid mass " << acid_ptr->getMonoMass());
        LOG_DEBUG("empty ptm mass " << ptm_ptr->getMonoMass());
        LOG_DEBUG("empty residue mass " << empty_residue_ptr->getMass());
        if (empty_residue_ptr == nullptr) {
          LOG_ERROR( "Proteoform:: residue not found");
          throw("Residue not found");
        }
        res_vec_2d_[0].push_back(empty_residue_ptr);
        change_vec_2d_[0].push_back(Change::getProteinVariableChange());
      }
    }
    else if (mod_ptr == ProtModFactory::getProtModPtr_NME_ACETYLATION()) {
      LOG_DEBUG("NME_ACETYLATION");
      TruncPtr trunc_ptr = mod_ptr->getTruncPtr();
      bool valid_trunc = trunc_ptr->isValidTrunc(residue_ptr_vec);
      if (valid_trunc && residue_ptr_vec.size() >= 2) {
        // add acetylation to the second residue
        is_nme_ = true;
        ResiduePtr second_residue_ptr = residue_ptr_vec[1];
        PtmPtr ori_ptm_ptr = second_residue_ptr->getPtmPtr();
        PtmPtr prot_ptm_ptr = mod_ptr->getPtmPtr();
        LOG_DEBUG("ptm ptr");
        /* if they are different */
        if (ori_ptm_ptr != prot_ptm_ptr) {
          /* add protein n-terminal mod */
          AcidPtr acid_ptr = second_residue_ptr->getAcidPtr();
          ResiduePtr mut_residue_ptr = ResidueFactory::addBaseResidue(acid_ptr, prot_ptm_ptr);
          if (mut_residue_ptr == nullptr) {
            LOG_ERROR( "Proteoform:: residue not found");
            throw("Residue not found");
          }
          res_vec_2d_[1].push_back(mut_residue_ptr);
          change_vec_2d_[1].push_back(Change::getProteinVariableChange());
        }
      }
    }
    LOG_DEBUG("round complete");
  }
  LOG_DEBUG("protein mod complete");

  //variable ptms
  for (size_t i = 0; i < acid_ptr_vec.size(); i++) {
    AcidPtr acid_ptr = acid_ptr_vec[i];
    // if exist modified residues
    if(ptm_map_.count(acid_ptr)) {
      ResiduePtrVec mod_res_vec = ptm_map_[acid_ptr];
      //remove duplications
      ResiduePtrVec exist_res_vec = res_vec_2d_[i];
      for (size_t j = 0; j < mod_res_vec.size(); j++) {
        bool found = false;
        for (size_t k = 0; k < exist_res_vec.size(); k++) {
          if (mod_res_vec[j] == exist_res_vec[k]) {
            found = true;
            break;
          }
        }
        if (!found) {
          res_vec_2d_[i].push_back(mod_res_vec[j]);
          change_vec_2d_[i].push_back(Change::getVariableChange());
        }

      }
    }
  }
  LOG_DEBUG("variable mod complete");
}

}
