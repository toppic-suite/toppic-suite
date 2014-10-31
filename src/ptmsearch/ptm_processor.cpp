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

  all_writer_ptr_ = PrsmWriterPtr (new PrsmWriter (output_file_name));
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
    std::string file_name = output_file_name+"_"+ convertToString(s)
        +"_"+ SemiAlignTypeFactory::getCompletePtr()->getName();
    PrsmWriterPtr complete_writer_ptr = PrsmWriterPtr (new PrsmWriter (file_name));
    complete_writer_ptrs_.push_back(complete_writer_ptr);
    file_name = output_file_name+"_"+ convertToString(s)
        +"_"+ SemiAlignTypeFactory::getPrefixPtr()->getName();
    PrsmWriterPtr prefix_writer_ptr = PrsmWriterPtr (new PrsmWriter (file_name));
    prefix_writer_ptrs_.push_back(prefix_writer_ptr);
    file_name = output_file_name+"_"+ convertToString(s)
        +"_"+ SemiAlignTypeFactory::getSuffixPtr()->getName();
    PrsmWriterPtr suffix_writer_ptr = PrsmWriterPtr (new PrsmWriter (file_name));
    suffix_writer_ptrs_.push_back(suffix_writer_ptr);
    file_name = output_file_name+"_"+ convertToString(s)
        +"_"+ SemiAlignTypeFactory::getInternalPtr()->getName();
    PrsmWriterPtr internal_writer_ptr = PrsmWriterPtr (new PrsmWriter (file_name));
    internal_writer_ptrs_.push_back(internal_writer_ptr);
  }
}

// close writers
void PtmProcessor::closeWriters() {
  all_writer_ptr_->close();
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
      complete_writer_ptrs_[s-1]->close();
      prefix_writer_ptrs_[s-1]->close();
      suffix_writer_ptrs_[s-1]->close();
      internal_writer_ptrs_[s-1]->close();
  }
}

// initialize data 
void PtmProcessor::initData() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  proteo_ptrs_ = readFastaToProteoform(
      prsm_para_ptr->getSearchDbFileName(), 
      prsm_para_ptr->getFixModResiduePtrVec());
  mod_proteo_2d_ptrs_ =  
      generate2DProtModProteoform(proteo_ptrs_, prsm_para_ptr->getAllowProtModPtrVec());

  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string simplePrsmFileName = basename(sp_file_name)
      + "." + mng_ptr_->input_file_ext_;
  simple_prsm_ptrs_  = readSimplePrsms(simplePrsmFileName);
  // find sequences for simple prsms
  for(size_t i =0;i< simple_prsm_ptrs_.size();i++){
    simple_prsm_ptrs_[i]->assignProteoformPtr(proteo_ptrs_, mod_proteo_2d_ptrs_);
  }
  comp_shift_ptr_ = CompShiftLowMemPtr(new CompShiftLowMem());
}

// process ptm search
void PtmProcessor::process(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  
  std::string log_file_name = mng_ptr_->prsm_para_ptr_->getLogFileName();

  std::ofstream logfile;

  if (log_file_name.length() != 0){
    logfile.open(log_file_name, std::ios::out | std::ios::app);
  }    
  
  int spectra_num = countSpNum (sp_file_name);
  MsAlignReader sp_reader(sp_file_name);
  DeconvMsPtr deconv_sp;
  int cnt = 0;
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  double total_time = 0;
  while((deconv_sp = sp_reader.getNextMs())!= nullptr){
    long start_t = clock();
    cnt++;
    SpectrumSetPtr spectrum_set_ptr = getSpectrumSet(deconv_sp, 0, sp_para_ptr);
    if(spectrum_set_ptr != nullptr){

      SimplePrsmPtrVec selected_prsm_ptrs 
          = getMatchedSimplePrsms(simple_prsm_ptrs_,deconv_sp->getHeaderPtr());
      processOneSpectrum(spectrum_set_ptr, selected_prsm_ptrs);

    }
    long stop_t = clock();
    double time = (stop_t - start_t) /double (CLOCKS_PER_SEC);
    total_time += time;
    std::cout << std::flush <<  "PTM search is processing " << cnt 
        << " of " << spectra_num << " spectra.\r";
    
    if (log_file_name.length() != 0){
      if (logfile.is_open()) {
        logfile << 0.203 + (double) cnt / spectra_num * 0.17 << std::endl;
      }
    }
    
    //std::cout << std::flush << "Ptm searching is processing " << cnt 
    //    << " of " << spectra_num << " spectra." << " time " << time << " total time " << total_time << " seconds.\r ";

  }
  sp_reader.close();
  logfile.close();
  closeWriters();
  std::cout << std::endl;
}

inline void seleTopPrsms(const PrsmPtrVec &all_prsm_ptrs, 
                         PrsmPtrVec &sele_prsm_ptrs, int n_report) {
  int match_size = all_prsm_ptrs.size();
  if(all_prsm_ptrs.size()!=0){
    for(int r=0;r< n_report;r++){
      if(r >= match_size){
        break;
      }
      if(all_prsm_ptrs[r]->getMatchFragNum() > 0){
        sele_prsm_ptrs.push_back(all_prsm_ptrs[r]);
      }
    }
  }
  std::sort(sele_prsm_ptrs.begin(), sele_prsm_ptrs.end(), prsmProteoformIdUp);
}


void PtmProcessor::processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                                      SimplePrsmPtrVec simple_prsm_ptrs) {

  PtmSlowFilterPtr slow_filter_ptr = PtmSlowFilterPtr(
      new PtmSlowFilter(spectrum_set_ptr,simple_prsm_ptrs,comp_shift_ptr_,mng_ptr_));
  //LOG_DEBUG("init filter complete");
  for (int s = 1; s <= mng_ptr_->n_unknown_shift_; s++) {
    PrsmPtrVec complete_prsm_ptrs = slow_filter_ptr->getPrsms(
        s-1, SemiAlignTypeFactory::getCompletePtr());
    std::sort(complete_prsm_ptrs.begin(), complete_prsm_ptrs.end(), 
              prsmCompPreMatchFragDown);
    PrsmPtrVec sele_complete_prsm_ptrs;
    seleTopPrsms(complete_prsm_ptrs, sele_complete_prsm_ptrs, mng_ptr_->n_report_);
    complete_writer_ptrs_[s-1]->writeVector(sele_complete_prsm_ptrs);
    all_writer_ptr_->writeVector(sele_complete_prsm_ptrs);

    PrsmPtrVec prefix_prsm_ptrs = slow_filter_ptr->getPrsms(
        s-1, SemiAlignTypeFactory::getPrefixPtr());
    std::sort(prefix_prsm_ptrs.begin(), prefix_prsm_ptrs.end(), 
              prsmCompPreMatchFragDown);
    PrsmPtrVec sele_prefix_prsm_ptrs;
    seleTopPrsms(prefix_prsm_ptrs, sele_prefix_prsm_ptrs, mng_ptr_->n_report_);
    prefix_writer_ptrs_[s-1]->writeVector(sele_prefix_prsm_ptrs);
    all_writer_ptr_->writeVector(sele_prefix_prsm_ptrs);

    PrsmPtrVec suffix_prsm_ptrs = slow_filter_ptr->getPrsms(
        s-1, SemiAlignTypeFactory::getSuffixPtr());
    std::sort(suffix_prsm_ptrs.begin(), suffix_prsm_ptrs.end(), prsmMatchFragmentDown);
    PrsmPtrVec sele_suffix_prsm_ptrs;
    seleTopPrsms(suffix_prsm_ptrs, sele_suffix_prsm_ptrs, mng_ptr_->n_report_);
    suffix_writer_ptrs_[s-1]->writeVector(sele_suffix_prsm_ptrs);
    all_writer_ptr_->writeVector(sele_suffix_prsm_ptrs);

    PrsmPtrVec internal_prsm_ptrs = slow_filter_ptr->getPrsms(
        s-1, SemiAlignTypeFactory::getInternalPtr());
    std::sort(internal_prsm_ptrs.begin(), internal_prsm_ptrs.end(), prsmMatchFragmentDown);
    PrsmPtrVec sele_internal_prsm_ptrs;
    seleTopPrsms(internal_prsm_ptrs, sele_internal_prsm_ptrs, mng_ptr_->n_report_);
    internal_writer_ptrs_[s-1]->writeVector(sele_internal_prsm_ptrs);
    all_writer_ptr_->writeVector(sele_internal_prsm_ptrs);
  }
}


} /* namespace prot */
