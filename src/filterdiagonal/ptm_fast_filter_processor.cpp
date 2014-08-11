#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_writer.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"

namespace prot {

PtmFastFilterProcessor::PtmFastFilterProcessor(PtmFastFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  ProteoformPtrVec proteoform_ptrs 
      = readFastaToProteoform(prsm_para_ptr->getSearchDbFileName(),
                              prsm_para_ptr->getFixModResiduePtrVec());
  LOG_DEBUG("start init filter.");
  filter_ptr_ = PtmFastFilterBlockPtr(new PtmFastFilterBlock(proteoform_ptrs, mng_ptr_));
  LOG_DEBUG("init filter is done.");
}

void PtmFastFilterProcessor::process(){
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int n_spectrum = countSpNum(sp_file_name);
  for(int i=0;i< filter_ptr_->getBlockSize();i++){
    std::cout << "Fast filtering block " << (i+1) << " out of " 
        << filter_ptr_->getBlockSize() << " starts" << std::endl; 
    processBlock(i, sp_file_name, n_spectrum);
  }
  combineBlock();
}

void PtmFastFilterProcessor::processBlock(int block, const std::string &sp_file_name,
                                          int n_spectra){
  filter_ptr_->initBlock(block);
  MsAlignReader reader(sp_file_name);
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string output_file_name = basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+ std::to_string(block);
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
    std::cout << std::flush << "Fast filtering block " << (block +1) 
        << " is processing " << cnt << " of " << n_spectra << " spectra.\r";
  }
  reader.close();
  writer.close();
  std::cout << std::endl << "Fast filtering block " << (block +1) 
      << " finished. " << std::endl;
}

void PtmFastFilterProcessor::combineBlock(){
  SimplePrsmPtrVec2D match_ptrs;

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();

  for(int i=0;i<filter_ptr_->getBlockSize();i++){
    std::string block_file_name = basename(sp_file_name) + 
        "." + mng_ptr_->output_file_ext_+"_"+std::to_string(i);
    match_ptrs.push_back(readSimplePrsms(block_file_name));
  }

  std::vector<int> pointers(filter_ptr_->getBlockSize());

  MsAlignReader reader(sp_file_name);
  std::string output_file_name = basename(sp_file_name) 
      + "." + mng_ptr_->output_file_ext_+"_COMBINED";
  SimplePrsmWriter writer(output_file_name);
  DeconvMsPtr deconv_ms_ptr;
  while((deconv_ms_ptr = reader.getNextMs()) != nullptr){
    SimplePrsmPtrVec selected_match_ptrs;
    for(size_t i = 0;i<match_ptrs.size();i++){
      for(size_t j = pointers[i]; j <match_ptrs[i].size(); j++){
        if(match_ptrs[i][j]->isSameSpectrum(deconv_ms_ptr->getHeaderPtr())){
          selected_match_ptrs.push_back(match_ptrs[i][j]);
        }
        if (match_ptrs[i][j]->isLargerSpectrumId(deconv_ms_ptr->getHeaderPtr())) {
          pointers[i] = j;
          break;
        }
      }
    }
    std::sort(selected_match_ptrs.begin(),selected_match_ptrs.end(),simplePrsmDown);
    SimplePrsmPtrVec result_match_ptrs;
    for(size_t i=0; i < selected_match_ptrs.size(); i++){
      if( i >= mng_ptr_-> ptm_fast_filter_result_num_){
        break;
      }
      result_match_ptrs.push_back(selected_match_ptrs[i]);
    }

    writer.write(result_match_ptrs);
  }
  reader.close();
  writer.close();
}

} /* namespace prot */
