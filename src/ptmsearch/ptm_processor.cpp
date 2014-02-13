/*
 * ptm_processor.cpp
 *
 *  Created on: Dec 20, 2013
 *      Author: xunlikun
 */

#include "base/prot_mod.hpp"
#include "base/fasta_reader.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {

PtmProcessor::PtmProcessor(PtmMngPtr mng){
  mng_ = mng;
  init();
}

void PtmProcessor::init(){
  seqs_ = prot::readFastaToProteoform(
      mng_->search_db_file_name_,
      mng_->base_data_->getFixModResiduePtrVec());
  std::string sp_file_name = mng_->spectrum_file_name_;
  std::string simplePrsmFileName = mng_->spectrum_file_name_
      + ".FILTER" + mng_->input_file_ext_;
  simplePrsms_  = prot::readSimplePrSM(simplePrsmFileName.c_str());
  prsmFindSeq(simplePrsms_,seqs_);
  comp_shift_ = CompShiftLowMemPtr(new CompShiftLowMem());
}

void PtmProcessor::prsmFindSeq(SimplePrSMPtrVec simple_prsms,ProteoformPtrVec seqs){
  for(unsigned int i =0;i<simple_prsms.size();i++){
    simple_prsms[i]->findSeq(seqs);
  }
}

void PtmProcessor::process(){
  std::string sp_file_name = mng_->spectrum_file_name_;
  std::string output_file_name = sp_file_name+"."+mng_->output_file_ext_;

  MsAlignReader spReader(sp_file_name);

  all_writer_= PrSMWriterPtr(new PrSMWriter(output_file_name));

  for(int i=0;i<mng_->n_unknown_shift_;i++){
    std::string file_name = output_file_name+"_"+ convertToString(i+1)
        +"_"+ SemiAlignTypeFactory::getCompletePtr()->getName();
    complete_writers_.push_back(PrSMWriterPtr(new PrSMWriter(file_name)));
    file_name = output_file_name+"_"+ convertToString(i+1)
        +"_"+ SemiAlignTypeFactory::getPrefixPtr()->getName();
    prefix_writers_.push_back(PrSMWriterPtr(new PrSMWriter(file_name)));
    file_name = output_file_name+"_"+ convertToString(i+1)
        +"_"+ SemiAlignTypeFactory::getSuffixPtr()->getName();
    suffix_writers_.push_back(PrSMWriterPtr(new PrSMWriter(file_name)));
    file_name = output_file_name+"_"+ convertToString(i+1)
        +"_"+ SemiAlignTypeFactory::getInternalPtr()->getName();
    internal_writers_.push_back(PrSMWriterPtr(new PrSMWriter(file_name)));
  }

  DeconvMsPtr deconv_sp;

  int cnt = 0;
  while((deconv_sp = spReader.getNextMs())!= nullptr){
    cnt++;
    SpectrumSetPtr spectrum_set_ptr = getSpectrumSet(deconv_sp,0,
                                                mng_->sp_para_);
    //set deconvMsPtr errorTolerance ??
    deconv_sp->getHeaderPtr()->setErrorToleranceByPpo(mng_->ppo_);
    if(spectrum_set_ptr != nullptr){
      SimplePrSMPtrVec selected_prsms 
          = prot::findSimplePrsms(simplePrsms_,deconv_sp->getHeaderPtr());
      search(spectrum_set_ptr, selected_prsms);
    }
  }
  spReader.close();
}

void PtmProcessor::choosePrsms(PrSMPtrVec &all_prsms, PrSMPtrVec &sele_prsms) {
  // sort prms TO DO likun
  if(all_prsms.size()!=0){
    for(int r=0;r<mng_->n_report_;r++){
      int match_size = all_prsms.size();
      if(r >= match_size){
        break;
      }
      if(all_prsms[r]->getMatchFragNum() > 0){
        sele_prsms.push_back(all_prsms[r]);
      }
    }
  }
  //sort again
}


void PtmProcessor::search(SpectrumSetPtr spectrum_set_ptr, 
                          SimplePrSMPtrVec matches) {

  PtmSlowFilterPtr slow_filter = PtmSlowFilterPtr(
      new PtmSlowFilter(spectrum_set_ptr,matches,comp_shift_,mng_));
  for (int s = 1; s <= mng_->n_unknown_shift_; s++) {
    PrSMPtrVec complete_prsms = slow_filter->getPrSMs(
        s, SemiAlignTypeFactory::getCompletePtr());
    PrSMPtrVec sele_complete_prsms;
    choosePrsms(complete_prsms, sele_complete_prsms);
    complete_writers_[s-1]->writeVector(sele_complete_prsms);
    all_writer_->writeVector(sele_complete_prsms);

    PrSMPtrVec prefix_prsms = slow_filter->getPrSMs(
        s, SemiAlignTypeFactory::getPrefixPtr());
    PrSMPtrVec sele_prefix_prsms;
    choosePrsms(prefix_prsms, sele_prefix_prsms);
    complete_writers_[s - 1]->writeVector(sele_prefix_prsms);
    all_writer_->writeVector(sele_prefix_prsms);

    PrSMPtrVec suffix_prsms = slow_filter->getPrSMs(
        s, SemiAlignTypeFactory::getSuffixPtr());
    PrSMPtrVec sele_suffix_prsms;
    choosePrsms(suffix_prsms, sele_suffix_prsms);
    complete_writers_[s - 1]->writeVector(sele_suffix_prsms);
    all_writer_->writeVector(sele_suffix_prsms);

    PrSMPtrVec internal_prsms = slow_filter->getPrSMs(
        s, SemiAlignTypeFactory::getInternalPtr());
    PrSMPtrVec sele_internal_prsms;
    choosePrsms(internal_prsms, sele_internal_prsms);
    complete_writers_[s - 1]->writeVector(sele_internal_prsms);
    all_writer_->writeVector(sele_internal_prsms);

    /* to do
    PrSMPtrVec prefix_prsms = slow_filter->getBestMatch(s, 
        SemiAlignTypeFactory::getPrefixPtr());
    choosePrsms(prefix_slow_matches, prefix_prsms, s, 1);
    prefix_writes[s-1]->writeVector(prefix_prsms);
    all_writers->writeVector(prefix_prsms);

    PrSMPtrVec suffix_prsms;
    PtmSlowMatchPtrVec suffix_slow_matches = slow_filter->getBestMatch(s, 2);
    choosePrsms(suffix_slow_matches, suffix_prsms, s, 2);
    suffix_writes[s-1]->writeVector(suffix_prsms);
    all_writers->writeVector(suffix_prsms);

    PrSMPtrVec internal_prsms;
    PtmSlowMatchPtrVec internal_slow_matches = slow_filter->getBestMatch(s, 3);
    choosePrsms(internal_slow_matches, internal_prsms, s, 3);
    internal_writes[s-1]->writeVector(internal_prsms);
    all_writers->writeVector(internal_prsms);
    */
  }
}

} /* namespace prot */
