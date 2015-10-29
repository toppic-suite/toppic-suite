#include "base/prot_mod.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "base/web_logger.hpp"
#include "spec/deconv_ms.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/prsm.hpp"
#include "oneptmsearch/one_ptm_search_processor.hpp"
#include "oneptmsearch/one_ptm_slow_filter.hpp"

namespace prot {

OnePtmSearchProcessor::OnePtmSearchProcessor(OnePtmSearchMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  initWriters();
  initData();
}

OnePtmSearchProcessor::~OnePtmSearchProcessor() {
  fai_destroy(fai_);
}

void OnePtmSearchProcessor::initReaders(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = basename(sp_file_name)+"."+mng_ptr_->input_file_ext_;

  std::string file_name = input_file_name+"_"
      + SemiAlignTypeFactory::getCompletePtr()->getName();
  complete_reader_ptr_ = SimplePrsmReaderPtr(new SimplePrsmReader (file_name));
  file_name = input_file_name+"_"
      + SemiAlignTypeFactory::getPrefixPtr()->getName();
  prefix_reader_ptr_ = SimplePrsmReaderPtr(new SimplePrsmReader (file_name));
  file_name = input_file_name+"_"
      + SemiAlignTypeFactory::getSuffixPtr()->getName();
  suffix_reader_ptr_ = SimplePrsmReaderPtr(new SimplePrsmReader (file_name));
  file_name = input_file_name+"_"
      + SemiAlignTypeFactory::getInternalPtr()->getName();
  internal_reader_ptr_ = SimplePrsmReaderPtr(new SimplePrsmReader (file_name));
}


// initialize writers 
void OnePtmSearchProcessor::initWriters(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string output_file_name = basename(sp_file_name)+"."+mng_ptr_->output_file_ext_;

  std::string file_name = output_file_name+"_"
      + SemiAlignTypeFactory::getCompletePtr()->getName();
  complete_writer_ptr_ = PrsmWriterPtr (new PrsmWriter (file_name));
  file_name = output_file_name+"_"
      + SemiAlignTypeFactory::getPrefixPtr()->getName();
  prefix_writer_ptr_ = PrsmWriterPtr (new PrsmWriter (file_name));
  file_name = output_file_name+"_"
      + SemiAlignTypeFactory::getSuffixPtr()->getName();
  suffix_writer_ptr_ = PrsmWriterPtr (new PrsmWriter (file_name));
  file_name = output_file_name+"_"
      + SemiAlignTypeFactory::getInternalPtr()->getName();
  internal_writer_ptr_ = PrsmWriterPtr (new PrsmWriter (file_name));
}

// close readers
void OnePtmSearchProcessor::closeReaders() {
  complete_reader_ptr_->close();
  prefix_reader_ptr_->close();
  suffix_reader_ptr_->close();
  internal_reader_ptr_->close();
}

// close writers
void OnePtmSearchProcessor::closeWriters() {
  complete_writer_ptr_->close();
  prefix_writer_ptr_->close();
  suffix_writer_ptr_->close();
  internal_writer_ptr_->close();
}

// initialize data 
void OnePtmSearchProcessor::initData() {
  comp_shift_ptr_ = CompShiftLowMemPtr(new CompShiftLowMem());
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  fai_ = fai_load(prsm_para_ptr->getSearchDbFileName().c_str());
}

// process ptm search
void OnePtmSearchProcessor::process(){
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  initReaders();
  SimplePrsmPtr comp_prsm_ptr = complete_reader_ptr_->readOnePrsm();
  SimplePrsmPtr pref_prsm_ptr = prefix_reader_ptr_->readOnePrsm();
  SimplePrsmPtr suff_prsm_ptr = suffix_reader_ptr_->readOnePrsm();
  SimplePrsmPtr internal_prsm_ptr = internal_reader_ptr_->readOnePrsm();

  //init variables
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int spectrum_num = getSpNum (sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ResiduePtrVec residue_ptr_vec = prsm_para_ptr->getFixModResiduePtrVec();
  ProtModPtrVec prot_mod_ptr_vec = prsm_para_ptr->getAllowProtModPtrVec();

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num);
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  //LOG_DEBUG("Start search");
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      SimplePrsmPtrVec comp_selected_prsm_ptrs;
      int spec_id = spec_set_ptr->getSpecId();
      while (comp_prsm_ptr != nullptr && comp_prsm_ptr->getSpectrumId() == spec_id) {
        comp_prsm_ptr->addProteoformPtr(fai_, residue_ptr_vec, prot_mod_ptr_vec);
        comp_selected_prsm_ptrs.push_back(comp_prsm_ptr);
        comp_prsm_ptr = complete_reader_ptr_->readOnePrsm();
      }

      if (comp_selected_prsm_ptrs.size() > 0) {
        //LOG_DEBUG("start processing one spectrum.");
        SemiAlignTypePtr type_ptr=SemiAlignTypeFactory::getCompletePtr(); 
        processOneSpectrum(spec_set_ptr, comp_selected_prsm_ptrs, type_ptr);
      }
    }
    std::cout << std::flush <<  "PTM search is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
    //WebLog::percentLog(cnt, spectrum_num, WebLog::PtmTime());
  }
  LOG_DEBUG("Search completed");
  sp_reader.close();
  closeReaders();
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


void OnePtmSearchProcessor::processOneSpectrum(SpectrumSetPtr spectrum_set_ptr, 
                                               SimplePrsmPtrVec simple_prsm_ptrs,
                                               SemiAlignTypePtr align_ptr) {
  OnePtmSlowFilterPtr slow_filter_ptr = OnePtmSlowFilterPtr(
      new OnePtmSlowFilter(spectrum_set_ptr,simple_prsm_ptrs,comp_shift_ptr_, align_ptr, mng_ptr_));

  /*
  PrsmPtrVec complete_prsm_ptrs = slow_filter_ptr->getPrsms(SemiAlignTypeFactory::getCompletePtr());
  std::sort(complete_prsm_ptrs.begin(), complete_prsm_ptrs.end(), 
            prsmCompPreMatchFragDown);
  PrsmPtrVec sele_complete_prsm_ptrs;
  seleTopPrsms(complete_prsm_ptrs, sele_complete_prsm_ptrs, mng_ptr_->n_report_);
  complete_writer_ptr_->writeVector(sele_complete_prsm_ptrs);

  PrsmPtrVec prefix_prsm_ptrs = slow_filter_ptr->getPrsms(SemiAlignTypeFactory::getPrefixPtr());
  std::sort(prefix_prsm_ptrs.begin(), prefix_prsm_ptrs.end(), 
            prsmCompPreMatchFragDown);
  PrsmPtrVec sele_prefix_prsm_ptrs;
  seleTopPrsms(prefix_prsm_ptrs, sele_prefix_prsm_ptrs, mng_ptr_->n_report_);
  prefix_writer_ptr_->writeVector(sele_prefix_prsm_ptrs);

  PrsmPtrVec suffix_prsm_ptrs = slow_filter_ptr->getPrsms(SemiAlignTypeFactory::getSuffixPtr());
  std::sort(suffix_prsm_ptrs.begin(), suffix_prsm_ptrs.end(), prsmMatchFragmentDown);
  PrsmPtrVec sele_suffix_prsm_ptrs;
  seleTopPrsms(suffix_prsm_ptrs, sele_suffix_prsm_ptrs, mng_ptr_->n_report_);
  suffix_writer_ptr_->writeVector(sele_suffix_prsm_ptrs);

  PrsmPtrVec internal_prsm_ptrs = slow_filter_ptr->getPrsms(SemiAlignTypeFactory::getInternalPtr());
  std::sort(internal_prsm_ptrs.begin(), internal_prsm_ptrs.end(), prsmMatchFragmentDown);
  PrsmPtrVec sele_internal_prsm_ptrs;
  seleTopPrsms(internal_prsm_ptrs, sele_internal_prsm_ptrs, mng_ptr_->n_report_);
  internal_writer_ptr_->writeVector(sele_internal_prsm_ptrs);
  */
}


} /* namespace prot */
