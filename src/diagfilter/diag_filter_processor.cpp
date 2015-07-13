#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "base/web_logger.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "diagfilter/diag_filter_processor.hpp"

namespace prot {

DiagFilterProcessor::DiagFilterProcessor(DiagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
}

void DiagFilterProcessor::process(){
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = readDbBlockIndex(db_file_name);

  for(size_t i=0; i< db_block_ptr_vec.size(); i++){
    std::cout << "Diagonal filtering block " << (i+1) << " out of " 
        << db_block_ptr_vec.size() << " started." << std::endl; 
    processBlock(db_block_ptr_vec[i], db_block_ptr_vec.size());
    std::cout << "Diagonal filtering block " << (i +1) 
        << " finished. " << std::endl;
  }

  std::cout << "Diagonal filtering: combining blocks started." << std::endl; 

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int block_num = db_block_ptr_vec.size();

  SimplePrsmStrCombinePtr combine_ptr(new SimplePrsmStrCombine(sp_file_name, mng_ptr_->output_file_ext_,
                                                               block_num, mng_ptr_->output_file_ext_, 
                                                               mng_ptr_->ptm_fast_filter_result_num_));
  combine_ptr->process();

  std::cout << "Diagonal filtering: combining blocks finished." << std::endl; 
}

void DiagFilterProcessor::processBlock(DbBlockPtr block_ptr, int total_block_num) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName() 
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(db_block_file_name, 
                              prsm_para_ptr->getFixModResiduePtrVec(),
                              block_ptr->getSeqIdx());
  DiagFilterPtr filter_ptr(new DiagFilter(raw_forms, mng_ptr_));

  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName());
  std::string output_file_name = basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+ std::to_string(block_ptr->getBlockIdx());
  
     
  SimplePrsmWriter writer(output_file_name);
  DeconvMsPtr deconv_ms_ptr;
  int spectrum_num = getSpNum (prsm_para_ptr->getSpectrumFileName());
  int cnt = 0;
  while((deconv_ms_ptr = reader.getNextMs()) != nullptr){
    cnt++;
    SpectrumSetPtr spectrum_set_ptr = getSpectrumSet(deconv_ms_ptr,0,
                                                     prsm_para_ptr->getSpParaPtr());
    if(spectrum_set_ptr != nullptr){
      PrmMsPtr ms_ptr = spectrum_set_ptr->getMsTwoPtr();
      SimplePrsmPtrVec match_ptrs = filter_ptr->getBestMatch(ms_ptr);
      writer.write(match_ptrs);
    }
    
    WebLog::percentLog(cnt, spectrum_num, block_ptr->getBlockIdx(), total_block_num, 
                       WebLog::DiagFilterTime());
    
    std::cout << std::flush << "Diagonal filtering is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
  }
  std::cout << std::endl;
  reader.close();
  writer.close();

}

} /* namespace prot */
