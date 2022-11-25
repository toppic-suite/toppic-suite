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

#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "search/varptmsearch/var_ptm_slow_match.hpp"
#include "search/varptmsearch/var_ptm_search_processor.hpp"

namespace toppic {

ProtModPtr findMatchedProtModPtr(ProtModPtrVec &mod_ptr_vec, double n_term_shift) {
  for (size_t i = 0; i < mod_ptr_vec.size(); i++) {
    double mod_shift = mod_ptr_vec[i]->getProtShift();
    if (mod_shift == n_term_shift) {
      return mod_ptr_vec[i];
    }
  }
  LOG_ERROR("Failed to find protein modification!" << n_term_shift);
  exit(EXIT_FAILURE);
}

int findMatchedPos(ProteoformPtr proteo_ptr, double n_term_shift) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  for (size_t i = 0; i < prms.size(); i++) {
    if (prms[i] == n_term_shift) {
      return i;
    }
  }
  LOG_ERROR("Failed to find start posisiton!" << n_term_shift);
  exit(EXIT_FAILURE);
}

PrsmPtrVec VarPtmSearchProcessor::varPtmSearchOneSpec(SpectrumSetPtr spec_set_ptr,
                                                      const SimplePrsmPtrVec &simple_prsm_ptr_vec,
                                                      FastaIndexReaderPtr reader_ptr,
                                                      VarPtmSearchMngPtr mng_ptr,
                                                      ProteoformTypePtr type_ptr) {
  ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ProtModPtrVec prot_mod_ptr_vec = mng_ptr->prsm_para_ptr_->getProtModPtrVec();
  ProteoformPtrVec proteoform_ptr_vec;
  for (size_t i = 0; i < simple_prsm_ptr_vec.size(); i++) {
    SimplePrsmPtr simple_prsm_ptr = simple_prsm_ptr_vec[i];
    std::string seq_name = simple_prsm_ptr->getSeqName();
    std::string seq_desc = simple_prsm_ptr->getSeqDesc();
    ProteoformPtr db_proteo_ptr
        = proteoform_factory::readFastaToProteoformPtr(reader_ptr, seq_name, seq_desc, fix_mod_list);
    double n_term_shift = simple_prsm_ptr->getNTermShifts()[0];
    double c_term_shift = simple_prsm_ptr->getCTermShifts()[0];
    
    ProteoformPtr proteo_with_prot_mod_ptr = db_proteo_ptr;
    // start position is relative to the database sequence
    int start_pos = 0;
    if (type_ptr == ProteoformType::COMPLETE || type_ptr == ProteoformType::PREFIX) {
      if (n_term_shift != 0) {
        // Get the proteoform with N-terminal modificaiton
        ProtModPtr prot_mod_ptr = findMatchedProtModPtr(prot_mod_ptr_vec, n_term_shift);
        proteo_with_prot_mod_ptr = proteoform_factory::geneProtModProteoform(db_proteo_ptr, prot_mod_ptr);
        start_pos = proteo_with_prot_mod_ptr->getStartPos();
        LOG_DEBUG("Start pos with N terminal mod" << start_pos);
      }
    }
    else {
      //type_ptr == ProteoformType::SUFFIX || type_ptr == ProteoformType::INTERNAL
      start_pos = findMatchedPos(db_proteo_ptr, n_term_shift);
    }
    int end_pos = db_proteo_ptr->getEndPos();
    if (type_ptr == ProteoformType::PREFIX || type_ptr == ProteoformType::INTERNAL) {
      double proteo_minus_water_mass = db_proteo_ptr->getMinusWaterMass();
      double mass = proteo_minus_water_mass - c_term_shift;
      // the reported position is the position of matched breakpoint
      // the sequence end position is breakpoint position - 1
      end_pos = findMatchedPos(db_proteo_ptr, mass) - 1; 
    }
    
    // get the start and end positions relative to the proteoform with
    // N-terminal modification
    int local_start_pos = start_pos - proteo_with_prot_mod_ptr->getStartPos();
    int local_end_pos = end_pos - proteo_with_prot_mod_ptr->getStartPos();
    FastaSeqPtr fasta_seq_ptr = db_proteo_ptr->getFastaSeqPtr();
    ProteoformPtr sub_proteo_ptr = proteoform_factory::geneSubProteoform(proteo_with_prot_mod_ptr,
                                                                         fasta_seq_ptr,
                                                                         local_start_pos, 
                                                                         local_end_pos);

    proteoform_ptr_vec.push_back(sub_proteo_ptr);
  }

  PrsmPtrVec prsms;
  for (size_t i = 0; i < proteoform_ptr_vec.size(); i++) {
    VarPtmSlowMatch slow_match(proteoform_ptr_vec[i], spec_set_ptr, mng_ptr);
    PrsmPtr tmp = slow_match.compute();
    if (tmp != nullptr) {
      prsms.push_back(tmp);
    }
  }
  std::sort(prsms.begin(), prsms.end(), Prsm::cmpMatchFragmentDecMatchPeakDec);
  if (prsms.size() > 0) {
    prsms.erase(prsms.begin() + mng_ptr->n_report_, prsms.end());
  }
  return prsms;
}

void VarPtmSearchProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  SimplePrsmReader comp_prsm_reader(input_file_name + "_" + ProteoformType::COMPLETE->getName());
  SimplePrsmReader pref_prsm_reader(input_file_name + "_" + ProteoformType::PREFIX->getName());
  SimplePrsmReader suff_prsm_reader(input_file_name + "_" + ProteoformType::SUFFIX->getName());
  SimplePrsmReader internal_prsm_reader(input_file_name + "_" + ProteoformType::INTERNAL->getName());
  SimplePrsmPtr comp_prsm_ptr = comp_prsm_reader.readOnePrsm();
  SimplePrsmPtr pref_prsm_ptr = pref_prsm_reader.readOnePrsm();
  SimplePrsmPtr suff_prsm_ptr = suff_prsm_reader.readOnePrsm();
  SimplePrsmPtr internal_prsm_ptr = internal_prsm_reader.readOnePrsm();

  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_;
  PrsmXmlWriter comp_writer(output_file_name + "_" + ProteoformType::COMPLETE->getName());
  PrsmXmlWriter pref_writer(output_file_name + "_" + ProteoformType::PREFIX->getName());
  PrsmXmlWriter suff_writer(output_file_name + "_" + ProteoformType::SUFFIX->getName());
  PrsmXmlWriter internal_writer(output_file_name + "_" + ProteoformType::INTERNAL->getName());

  // init variables
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();

  FastaIndexReaderPtr reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);
  int spectrum_num = msalign_util::getSpNum(sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  SimpleMsAlignReaderPtr msalign_reader_ptr = std::make_shared<SimpleMsAlignReader>(sp_file_name, 
                                                                                    group_spec_num,
                                                                                    sp_para_ptr->getActivationPtr());
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;
  // LOG_DEBUG("Start search");
  while ((spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(msalign_reader_ptr,sp_para_ptr))!= nullptr) {
    cnt+= group_spec_num;
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      // complete
      SimplePrsmPtrVec comp_selected_prsm_ptrs;
      while (comp_prsm_ptr != nullptr && comp_prsm_ptr->getSpectrumId() == spec_id) {
        comp_selected_prsm_ptrs.push_back(comp_prsm_ptr);
        comp_prsm_ptr = comp_prsm_reader.readOnePrsm();
      }
      if (comp_selected_prsm_ptrs.size() > 0) {
        // LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms = varPtmSearchOneSpec(spec_set_ptr, comp_selected_prsm_ptrs,
                                               reader_ptr, mng_ptr_, ProteoformType::COMPLETE);
        comp_writer.writeVector(prsms);
      }

      // prefix
      SimplePrsmPtrVec pref_selected_prsm_ptrs;
      while (pref_prsm_ptr != nullptr && pref_prsm_ptr->getSpectrumId() == spec_id) {
        pref_selected_prsm_ptrs.push_back(pref_prsm_ptr);
        pref_prsm_ptr = pref_prsm_reader.readOnePrsm();
      }
      if (pref_selected_prsm_ptrs.size() > 0) {
        // LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms = varPtmSearchOneSpec(spec_set_ptr, pref_selected_prsm_ptrs,
                                               reader_ptr, mng_ptr_, ProteoformType::PREFIX);
        pref_writer.writeVector(prsms);
      }

      // suffix
      SimplePrsmPtrVec suff_selected_prsm_ptrs;
      while (suff_prsm_ptr != nullptr && suff_prsm_ptr->getSpectrumId() == spec_id) {
        suff_selected_prsm_ptrs.push_back(suff_prsm_ptr);
        suff_prsm_ptr = suff_prsm_reader.readOnePrsm();
      }
      if (suff_selected_prsm_ptrs.size() > 0) {
        // LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms = varPtmSearchOneSpec(spec_set_ptr, suff_selected_prsm_ptrs,
                                               reader_ptr, mng_ptr_, ProteoformType::SUFFIX);
        suff_writer.writeVector(prsms);
      }

      // internal
      SimplePrsmPtrVec internal_selected_prsm_ptrs;
      while (internal_prsm_ptr != nullptr && internal_prsm_ptr->getSpectrumId() == spec_id) {
        internal_selected_prsm_ptrs.push_back(internal_prsm_ptr);
        internal_prsm_ptr = internal_prsm_reader.readOnePrsm();
      }
      if (internal_selected_prsm_ptrs.size() > 0) {
        // LOG_DEBUG("start processing one spectrum.");
        PrsmPtrVec prsms = varPtmSearchOneSpec(spec_set_ptr, internal_selected_prsm_ptrs,
                                               reader_ptr, mng_ptr_, ProteoformType::INTERNAL);
        internal_writer.writeVector(prsms);
      }
    }
    std::cout << std::flush <<  "Variable PTM search - processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
  }
  int remainder = spectrum_num - cnt;
  if (prsm_para_ptr->getGroupSpecNum() > remainder && remainder > 0){
      //if there are spectrum remaining because they were not combined due to not having enough pairs
      //fix the message as the processing is completed.
      //this code avoids error when no combined spectra is used but a scan is remaining unprocessed
      //because then it will not satisfy the first condition
      std::cout << std::flush <<  "Variable PTM search - processing " << spectrum_num
        << " of " << spectrum_num << " spectra.\r";
  } 
  comp_prsm_reader.close();
  pref_prsm_reader.close();
  suff_prsm_reader.close();
  internal_prsm_reader.close();
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
  std::cout << std::endl;
}

}  // namespace toppic
