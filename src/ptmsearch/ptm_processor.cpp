
#include "base/prot_mod.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {

PtmProcessor::PtmProcessor(PtmMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  initWriters();
  initData();
}

PtmProcessor::~PtmProcessor() {
  fai_destroy(fai_);
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
  comp_shift_ptr_ = CompShiftLowMemPtr(new CompShiftLowMem());
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  fai_ = fai_load(prsm_para_ptr->getSearchDbFileName().c_str());
}

// process ptm search
void PtmProcessor::process(){
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = basename(sp_file_name)+"."+mng_ptr_->input_file_ext_;
  SimplePrsmReader simple_prsm_reader(input_file_name);
  SimplePrsmPtr prsm_ptr = simple_prsm_reader.readOnePrsm();

  //init variables
  int spectrum_num = getSpNum (sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ResiduePtrVec residue_ptr_vec = prsm_para_ptr->getFixModResiduePtrVec();
  ProtModPtrVec prot_mod_ptr_vec = prsm_para_ptr->getAllowProtModPtrVec();

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num);
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      SimplePrsmPtrVec selected_prsm_ptrs;
      int spec_id = spec_set_ptr->getSpecId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        prsm_ptr->addProteoformPtr(fai_, residue_ptr_vec, prot_mod_ptr_vec);
        selected_prsm_ptrs.push_back(prsm_ptr);
        prsm_ptr = simple_prsm_reader.readOnePrsm();
      }
      if (selected_prsm_ptrs.size() > 0) {
        //LOG_DEBUG("start processing one spectrum.");
        processOneSpectrum(spec_set_ptr, selected_prsm_ptrs);
      }
    }
    std::cout << std::flush <<  "PTM search is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
    
    WebLog::percent_log(0.09 + (double) cnt / spectrum_num * 0.17);
  }
  LOG_DEBUG("Search completed");
  sp_reader.close();
  simple_prsm_reader.close();
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
