#include "base/logger.hpp"
#include "base/file_util.hpp"
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/db_block.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

namespace prot {

void zeroPtmSearch(SpectrumSetPtr spec_set_ptr, 
                   SemiAlignTypePtr type_ptr,
                   ProteoformPtrVec &proteoform_ptr_vec, 
                   ZeroPtmMngPtr mng_ptr, 
                   PrsmPtrVec &prsms) {
  ExtendMsPtrVec ms_three_vec = spec_set_ptr->getMsThreePtrVec();
  //LOG_DEBUG("ms three vector size " << ms_three_vec.size());

  ZpFastMatchPtrVec fast_matches 
      = zeroPtmFastFilter(type_ptr, ms_three_vec, proteoform_ptr_vec, 
                          mng_ptr->zero_ptm_filter_result_num_);
  //LOG_DEBUG("fast_match ended size " << fast_matches.size());
  DeconvMsPtrVec deconv_ms_vec = spec_set_ptr->getDeconvMsPtrVec();
  ZpSlowMatchPtrVec slow_matches 
      = zeroPtmSlowFilter(deconv_ms_vec, fast_matches, mng_ptr); 

  //LOG_DEBUG("slow_match ended size " << slow_matches.size());
  for (size_t i = 0; i < slow_matches.size(); i++) {
      prsms.push_back(slow_matches[i]->geneResult());
  }
  //LOG_DEBUG("prsm generation ended size " << prsms.size());

  std::sort(prsms.begin(), prsms.end(), prsmMatchFragmentDown);
  if (prsms.size() > 0) {
    prsms.erase(prsms.begin() + 1, prsms.end());
  }
}

void zeroPtmSearchProcessBlock(ZeroPtmMngPtr mng_ptr, DbBlockPtr block_ptr, 
                               int total_block_num) { 
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName() 
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(db_block_file_name, 
                              prsm_para_ptr->getFixModResiduePtrVec(),
                              block_ptr->getSeqIdx());
  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, prsm_para_ptr->getAllowProtModPtrVec());

  int spectrum_num = getSpNum (prsm_para_ptr->getSpectrumFileName());

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(), group_spec_num);

  std::string output_file_name = basename(prsm_para_ptr->getSpectrumFileName())
                                          + "." + mng_ptr->output_file_ext_;  
  std::string block_str = "_" + std::to_string(block_ptr->getBlockIdx());

  PrsmWriter comp_writer(output_file_name + "_" 
                         + SemiAlignTypeFactory::getCompletePtr()->getName() 
                         + block_str);
  PrsmWriter pref_writer(output_file_name + "_"
                         + SemiAlignTypeFactory::getPrefixPtr()->getName()
                         + block_str);
  PrsmWriter suff_writer(output_file_name + "_"
                         + SemiAlignTypeFactory::getSuffixPtr()->getName()
                         + block_str);
  PrsmWriter internal_writer(output_file_name + "_"
                         + SemiAlignTypeFactory::getInternalPtr()->getName()
                         + block_str);
  PrsmWriter all_writer(output_file_name);

  //LOG_DEBUG("start reading");
  int n = 0;
  SpectrumSetPtr spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr);
  LOG_DEBUG("init ms_ptr");
  //double delta = 0;
  while (spec_set_ptr != nullptr) {
    n = n + group_spec_num;
    if (spec_set_ptr->isValid()) {
      PrsmPtrVec comp_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getCompletePtr(), 
                    prot_mod_forms, mng_ptr, comp_prsms);
      comp_writer.writeVector(comp_prsms);
      all_writer.writeVector(comp_prsms);
      PrsmPtrVec pref_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getPrefixPtr(), 
                    prot_mod_forms, mng_ptr, pref_prsms);
      pref_writer.writeVector(pref_prsms);
      all_writer.writeVector(pref_prsms);
      PrsmPtrVec suff_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getSuffixPtr(), 
                    raw_forms, mng_ptr, suff_prsms);
      suff_writer.writeVector(suff_prsms);
      all_writer.writeVector(suff_prsms);
      PrsmPtrVec internal_prsms;
      zeroPtmSearch(spec_set_ptr, SemiAlignTypeFactory::getInternalPtr(), 
                    raw_forms, mng_ptr, internal_prsms);
      internal_writer.writeVector(internal_prsms);
      all_writer.writeVector(internal_prsms);

      WebLog::percent_log( (double)block_ptr->getBlockIdx() / total_block_num * 0.03
              + (double) n / spectrum_num / total_block_num * 0.03);
          
      std::cout << std::flush << "Zero PTM search is processing " << n << " of " 
          << spectrum_num << " spectra.\r";
    }
    spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr);
    //LOG_DEBUG("spectrum " << n);
  }
  std::cout << std::endl;

  reader.close();

  //because the prsm_writer ~PrsmWriter changed and the fileclosing is an independant function
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
  all_writer.close();
}

void zeroPtmSearchProcess(ZeroPtmMngPtr mng_ptr) {
  std::string db_file_name = mng_ptr->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = readDbBlockIndex(db_file_name);

  for(size_t i=0; i< db_block_ptr_vec.size(); i++){
    std::cout << "Zero PTM search block " << (i+1) << " out of " 
        << db_block_ptr_vec.size() << " started." << std::endl; 
    zeroPtmSearchProcessBlock(mng_ptr, db_block_ptr_vec[i], db_block_ptr_vec.size());
    std::cout << "Zero PTM search block " << (i +1) 
        << " finished. " << std::endl;
  }

  std::cout << "Zero PTM search: combining blocks started." << std::endl; 
  std::string sp_file_name = mng_ptr->prsm_para_ptr_->getSpectrumFileName();
  std::string output_ext = mng_ptr->output_file_ext_;
  int block_num = db_block_ptr_vec.size();
  PrsmStrCombinePtr comp_combine_ptr(new PrsmStrCombine(sp_file_name, output_ext + "_COMPLETE",
                                                        block_num, output_ext + "_COMPLETE", 
                                                        mng_ptr->report_num_));
  comp_combine_ptr->process();
  comp_combine_ptr = nullptr;

  PrsmStrCombinePtr pref_combine_ptr(new PrsmStrCombine(sp_file_name, output_ext + "_PREFIX",
                                                        block_num, output_ext + "_PREFIX", 
                                                        mng_ptr->report_num_));
  pref_combine_ptr->process();
  pref_combine_ptr = nullptr;
  PrsmStrCombinePtr suff_combine_ptr(new PrsmStrCombine(sp_file_name, output_ext + "_SUFFIX",
                                                        block_num, output_ext + "_SUFFIX", 
                                                        mng_ptr->report_num_));
  suff_combine_ptr->process();
  suff_combine_ptr = nullptr;
  PrsmStrCombinePtr inter_combine_ptr(new PrsmStrCombine(sp_file_name, output_ext + "_INTERNAL",
                                                        block_num, output_ext + "_INTERNAL", 
                                                        mng_ptr->report_num_));
  inter_combine_ptr->process();
  inter_combine_ptr = nullptr;
  std::cout << "Zero PTM serach: combining blocks finished." << std::endl; 
}

} // end namespace
