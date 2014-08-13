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
  initWriters();
  initData();
}

// initialize writers 
void PtmProcessor::initWriters(){
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
}

// close writers
void PtmProcessor::closeWriters() {
  all_writer_->close();
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
      complete_writers_[s-1]->close();
      prefix_writers_[s-1]->close();
      suffix_writers_[s-1]->close();
      internal_writers_[s-1]->close();
  }
}

// initialize data 
void PtmProcessor::initData() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  proteoforms_ = readFastaToProteoform(
      prsm_para_ptr->getSearchDbFileName(), 
      prsm_para_ptr->getFixModResiduePtrVec());
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string simplePrsmFileName = basename(sp_file_name)
      + "." + mng_ptr_->input_file_ext_;
  simple_prsms_  = readSimplePrsms(simplePrsmFileName.c_str());
  // find sequences for simple prsms
  for(size_t i =0;i< simple_prsms_.size();i++){
    simple_prsms_[i]->assignProteoformPtr(proteoforms_);
  }
  comp_shift_ = CompShiftLowMemPtr(new CompShiftLowMem());
}

// process ptm search
void PtmProcessor::process(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  MsAlignReader spReader(sp_file_name);
  DeconvMsPtr deconv_sp;
  int spectra_num = countSpNum (sp_file_name);
  int cnt = 0;
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  while((deconv_sp = spReader.getNextMs())!= nullptr){
    cnt++;
    SpectrumSetPtr spectrum_set_ptr = getSpectrumSet(deconv_sp, 0, sp_para_ptr);
    if(spectrum_set_ptr != nullptr){
      SimplePrsmPtrVec selected_prsms 
          = getMatchedSimplePrsms(simple_prsms_,deconv_sp->getHeaderPtr());
      processOneSpectrum(spectrum_set_ptr, selected_prsms);
    }
    std::cout << std::flush << "Ptm searching is processing " << cnt 
        << " of " << spectra_num << " spectra.\r";
  }
  spReader.close();
  closeWriters();
  std::cout << std::endl << "Ptm searching finished." << std::endl;
}

void seleTopPrsms(PrsmPtrVec &all_prsms, PrsmPtrVec &sele_prsms, int n_report) {
  int match_size = all_prsms.size();
  if(all_prsms.size()!=0){
    for(int r=0;r< n_report;r++){
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


void PtmProcessor::processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                                      SimplePrsmPtrVec matches) {
  PtmSlowFilterPtr slow_filter = PtmSlowFilterPtr(
      new PtmSlowFilter(spectrum_set_ptr,matches,comp_shift_,mng_ptr_));
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
    PrsmPtrVec complete_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getCompletePtr());
    std::sort(complete_prsms.begin(), complete_prsms.end(), 
              prsmCompPreMatchFragDown);
    PrsmPtrVec sele_complete_prsms;
    seleTopPrsms(complete_prsms, sele_complete_prsms, mng_ptr_->n_report_);
    complete_writers_[s-1]->writeVector(sele_complete_prsms);
    all_writer_->writeVector(sele_complete_prsms);

    PrsmPtrVec prefix_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getPrefixPtr());
    std::sort(prefix_prsms.begin(), prefix_prsms.end(), 
              prsmCompPreMatchFragDown);
    PrsmPtrVec sele_prefix_prsms;
    seleTopPrsms(prefix_prsms, sele_prefix_prsms, mng_ptr_->n_report_);
    prefix_writers_[s-1]->writeVector(sele_prefix_prsms);
    all_writer_->writeVector(sele_prefix_prsms);

    PrsmPtrVec suffix_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getSuffixPtr());
    std::sort(suffix_prsms.begin(), suffix_prsms.end(), prsmMatchFragmentDown);
    PrsmPtrVec sele_suffix_prsms;
    seleTopPrsms(suffix_prsms, sele_suffix_prsms, mng_ptr_->n_report_);
    suffix_writers_[s-1]->writeVector(sele_suffix_prsms);
    all_writer_->writeVector(sele_suffix_prsms);

    PrsmPtrVec internal_prsms = slow_filter->getPrsms(
        s-1, SemiAlignTypeFactory::getInternalPtr());
    std::sort(internal_prsms.begin(), internal_prsms.end(), prsmMatchFragmentDown);
    PrsmPtrVec sele_internal_prsms;
    seleTopPrsms(internal_prsms, sele_internal_prsms, mng_ptr_->n_report_);
    internal_writers_[s-1]->writeVector(sele_internal_prsms);
    all_writer_->writeVector(sele_internal_prsms);
  }
}


} /* namespace prot */
