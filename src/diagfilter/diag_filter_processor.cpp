#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_writer.hpp"
#include "diagfilter/diag_filter_processor.hpp"

namespace prot {

DiagFilterProcessor::DiagFilterProcessor(DiagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  ProteoformPtrVec proteoform_ptrs 
      = readFastaToProteoform(prsm_para_ptr->getSearchDbFileName(),
                              prsm_para_ptr->getFixModResiduePtrVec());
  LOG_DEBUG("start init filter.");
  filter_ptr_ = DiagFilterBlockPtr(new DiagFilterBlock(proteoform_ptrs, mng_ptr_));
  LOG_DEBUG("init filter is done.");
}

void DiagFilterProcessor::process(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int n_spectrum = countSpNum(sp_file_name);
  for(int i=0;i< filter_ptr_->getBlockSize();i++){
    std::cout << "Fast filtering block " << (i+1) << " out of " 
        << filter_ptr_->getBlockSize() << " starts" << std::endl; 
    processBlock(i, sp_file_name, n_spectrum);
  }
  combineBlock(sp_file_name, filter_ptr_->getBlockSize(), mng_ptr_->output_file_ext_, 
               mng_ptr_->ptm_fast_filter_result_num_);
}

void DiagFilterProcessor::processBlock(int block, const std::string &sp_file_name,
                                       int n_spectra){
  filter_ptr_->initBlock(block);
  MsAlignReader reader(sp_file_name);
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string output_file_name = basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+ std::to_string(block);
      
  std::string log_file_name = prsm_para_ptr->getLogFileName();

  std::ofstream logfile;

  if (log_file_name.length() != 0){
    logfile.open(log_file_name, std::ios::out | std::ios::app);
  }
      
  SimplePrsmWriter writer(output_file_name);
  DeconvMsPtr deconv_ms_ptr;
  int cnt = 0;
  while((deconv_ms_ptr = reader.getNextMs()) != nullptr){
    cnt++;
    SpectrumSetPtr spectrum_set_ptr = getSpectrumSet(deconv_ms_ptr,0,
                                                     prsm_para_ptr->getSpParaPtr());
    if(spectrum_set_ptr != nullptr){
      std::string scan = deconv_ms_ptr->getHeaderPtr()->getScansString();
      SimplePrsmPtrVec match_ptrs = filter_ptr_->getBestMathBatch(spectrum_set_ptr);
      writer.write(match_ptrs);
    }
    
    if (log_file_name.length() != 0){
      if (logfile.is_open()) {
        logfile << 0.053 + (double) cnt / n_spectra * 0.075 << std::endl;
      }
    }
    
    std::cout << std::flush << "Fast filtering block " << (block +1) 
        << " is processing " << cnt << " of " << n_spectra << " spectra.\r";
  }
  reader.close();
  writer.close();
  logfile.close();
  std::cout << std::endl << "Fast filtering block " << (block +1) 
      << " finished. " << std::endl;
}


} /* namespace prot */
