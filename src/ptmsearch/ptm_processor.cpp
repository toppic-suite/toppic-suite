/*
 * ptm_processor.cpp
 *
 *  Created on: Dec 20, 2013
 *      Author: xunlikun
 */

#include "base/prot_mod.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {

PtmProcessor::PtmProcessor(PtmMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string output_file_name = basename(sp_file_name)+"."+mng_ptr_->output_file_ext_;

  all_writer_ = PrsmWriterPtr (new PrsmWriter (output_file_name));
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
    std::string file_name = output_file_name+"_"+ convertToString(s)
              +"_"+ SemiAlignTypeFactory::getCompletePtr()->getName();
    PrsmWriterPtr complete_writer= PrsmWriterPtr (new PrsmWriter (file_name));
    complete_writers_.push_back(complete_writer);
    file_name = output_file_name+"_"+ convertToString(s)
              +"_"+ SemiAlignTypeFactory::getPrefixPtr()->getName();
    PrsmWriterPtr prefix_writer= PrsmWriterPtr (new PrsmWriter (file_name));
    prefix_writers_.push_back(prefix_writer);
    file_name = output_file_name+"_"+ convertToString(s)
              +"_"+ SemiAlignTypeFactory::getSuffixPtr()->getName();
    PrsmWriterPtr suffix_writer= PrsmWriterPtr (new PrsmWriter (file_name));
    suffix_writers_.push_back(suffix_writer);
    file_name = output_file_name+"_"+ convertToString(s)
              +"_"+ SemiAlignTypeFactory::getInternalPtr()->getName();
    PrsmWriterPtr internal_writer= PrsmWriterPtr (new PrsmWriter (file_name));
    internal_writers_.push_back(internal_writer);
 }

  init();
}

void PtmProcessor::init(){
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  proteoforms_ = prot::readFastaToProteoform(
      prsm_para_ptr->getSearchDbFileName(), 
      prsm_para_ptr->getFixModResiduePtrVec());
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string simplePrsmFileName = basename(sp_file_name)
      + "." + mng_ptr_->input_file_ext_;
  simplePrsms_  = prot::readSimplePrsm(simplePrsmFileName.c_str());
  prsmFindSeq(simplePrsms_,proteoforms_);
  comp_shift_ = CompShiftLowMemPtr(new CompShiftLowMem());
}

void PtmProcessor::prsmFindSeq(SimplePrsmPtrVec simple_prsms,ProteoformPtrVec seqs){
  for(unsigned int i =0;i<simple_prsms.size();i++){
    simple_prsms[i]->findSeq(seqs);
  }
}

void PtmProcessor::process(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string output_file_name = basename(sp_file_name)+"."+mng_ptr_->output_file_ext_;

  MsAlignReader spReader(sp_file_name);
  DeconvMsPtr deconv_sp;
  int spectra_num = countSpNum (sp_file_name);
  int cnt = 0;
  while((deconv_sp = spReader.getNextMs())!= nullptr){
    cnt++;
    SpectrumSetPtr spectrum_set_ptr = getSpectrumSet(deconv_sp,0,
                                                mng_ptr_->prsm_para_ptr_->getSpParaPtr());
    if(spectrum_set_ptr != nullptr){
      SimplePrsmPtrVec selected_prsms 
          = findSimplePrsms(simplePrsms_,deconv_sp->getHeaderPtr());
      search(spectrum_set_ptr, selected_prsms);
    }
    std::cout << std::flush << "Ptm searching is processing " << cnt << " of " << spectra_num << " spectra.\r";
  }
  spReader.close();
  all_writer_->close();
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
      complete_writers_[s-1]->close();
      prefix_writers_[s-1]->close();
      suffix_writers_[s-1]->close();
      internal_writers_[s-1]->close();
  }
  std::cout << std::endl << "Ptm searching finished." << std::endl;
}

void PtmProcessor::chooseCompPrePrsms(PrsmPtrVec &all_prsms, PrsmPtrVec &sele_prsms) {
  int match_size = all_prsms.size();
  std::sort(all_prsms.begin(), all_prsms.end(), prsmCompPreMatchFragmentDown);
  if(all_prsms.size()!=0){
    for(int r=0;r<mng_ptr_->n_report_;r++){
      if(r >= match_size){
        break;
      }
      if(all_prsms[r]->getMatchFragNum() > 0){
        sele_prsms.push_back(all_prsms[r]);
      }
    }
  }
  std::sort(sele_prsms.begin(), sele_prsms.end(), prsmProteoformIdUp);
}

void PtmProcessor::chooseSuffIntPrsms(PrsmPtrVec &all_prsms, PrsmPtrVec &sele_prsms) {
  int match_size = all_prsms.size();
  std::sort(all_prsms.begin(), all_prsms.end(), prsmMatchFragmentDown);
  if(all_prsms.size()!=0){
    for(int r=0;r<mng_ptr_->n_report_;r++){
      if(r >= match_size){
        break;
      }
      if(all_prsms[r]->getMatchFragNum() > 0){
        sele_prsms.push_back(all_prsms[r]);
      }
    }
  }
  std::sort(sele_prsms.begin(), sele_prsms.end(), prsmProteoformIdUp);
}


void PtmProcessor::search(SpectrumSetPtr spectrum_set_ptr, 
                          SimplePrsmPtrVec matches) {
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string output_file_name = basename(sp_file_name)+"."+mng_ptr_->output_file_ext_;

  PtmSlowFilterPtr slow_filter = PtmSlowFilterPtr(
      new PtmSlowFilter(spectrum_set_ptr,matches,comp_shift_,mng_ptr_));
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
    PrsmPtrVec complete_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getCompletePtr());
    PrsmPtrVec sele_complete_prsms;
    chooseCompPrePrsms(complete_prsms, sele_complete_prsms);
    complete_writers_[s-1]->writeVector(sele_complete_prsms);
    all_writer_->writeVector(sele_complete_prsms);

    PrsmPtrVec prefix_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getPrefixPtr());
    PrsmPtrVec sele_prefix_prsms;
    chooseCompPrePrsms(prefix_prsms, sele_prefix_prsms);
    prefix_writers_[s-1]->writeVector(sele_prefix_prsms);
    all_writer_->writeVector(sele_prefix_prsms);

    PrsmPtrVec suffix_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getSuffixPtr());
    PrsmPtrVec sele_suffix_prsms;
    chooseSuffIntPrsms(suffix_prsms, sele_suffix_prsms);
    suffix_writers_[s-1]->writeVector(sele_suffix_prsms);
    all_writer_->writeVector(sele_suffix_prsms);

    PrsmPtrVec internal_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getInternalPtr());
    PrsmPtrVec sele_internal_prsms;
    chooseSuffIntPrsms(internal_prsms, sele_internal_prsms);
    internal_writers_[s-1]->writeVector(sele_internal_prsms);
    all_writer_->writeVector(sele_internal_prsms);


//    PrsmWriterPtr all_writer_plus= PrsmWriterPtr(new PrsmWriter("in/sp_bak.msalign_ALL_RESULT"));
//    all_writer_plus->writeVector(complete_prsms);
//    all_writer_plus->writeVector(prefix_prsms);
//    all_writer_plus->writeVector(suffix_prsms);
//    all_writer_plus->writeVector(internal_prsms);
  }
}

} /* namespace prot */
