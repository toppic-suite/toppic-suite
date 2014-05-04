/*
 * ptm_fast_filter_processor.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <sstream>

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
  ProteoformPtrVec proteoforms 
      = readFastaToProteoform(prsm_para_ptr->getSearchDbFileName(),
                              prsm_para_ptr->getFixModResiduePtrVec());
    filter_ptr_ = PtmFastFilterBlockPtr(new PtmFastFilterBlock(proteoforms,mng_ptr_));
}

void PtmFastFilterProcessor::process(){
    std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
    int n_spectrum = prot::countSpNum(sp_file_name.c_str());
    for(int i=0;i< filter_ptr_->getBlockSize();i++){
        processBlock(i,sp_file_name,n_spectrum);
    }
    combineBlock(sp_file_name);
}

void PtmFastFilterProcessor::processBlock(int block,std::string sp_file_name,
                                          int n_spectra){
    //system.out
    filter_ptr_->initBlock(block);
    MsAlignReader reader(sp_file_name);
    std::stringstream block_s;
    block_s<<block;
    PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
    std::string output_file_name = basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_+"_"+block_s.str();
    SimplePrSMWriter writer(output_file_name.c_str());
    DeconvMsPtr deconv_sp;
    int cnt = 0;
    while((deconv_sp = reader.getNextMs()) != nullptr){
      cnt++;
      SpectrumSetPtr spectrum_set = getSpectrumSet(deconv_sp,0,
                                                   prsm_para_ptr->getSpParaPtr());
      if(spectrum_set != nullptr){
        std::string scan = deconv_sp->getHeaderPtr()->getScansString();
        SimplePrSMPtrVec matches = filter_ptr_->getBestMathBatch(spectrum_set);
        writer.write(matches);
      }
    }
    reader.close();
    writer.close();
}

void PtmFastFilterProcessor::combineBlock(std::string sp_file_name){
    //system.out
    SimplePrSMPtrVec2D matches;

    //pravate readsimplePrsm
    for(int i=0;i<filter_ptr_->getBlockSize();i++){
        std::stringstream block_s;
        block_s<<i;
        std::string block_file_name = basename(mng_ptr_->prsm_para_ptr_->getSpectrumFileName()) 
        + "." + mng_ptr_->output_file_ext_+"_"+block_s.str();
        matches.push_back(prot::readSimplePrSM(block_file_name.c_str()));
    }

    MsAlignReader reader(sp_file_name.c_str());
    std::string output_file_name = basename(mng_ptr_->prsm_para_ptr_->getSpectrumFileName()) 
      + "." + mng_ptr_->output_file_ext_+"_COMBINED";
    SimplePrSMWriter writer(output_file_name.c_str());
    DeconvMsPtr deconv_sp;
    while((deconv_sp = reader.getNextMs()) != nullptr){
        SimplePrSMPtrVec selected_matches;
        for(unsigned int i=0;i<matches.size();i++){
            for(unsigned int j=0;j<matches[i].size();j++){
                if(matches[i][j]->isMatch(deconv_sp->getHeaderPtr())){
                    selected_matches.push_back(matches[i][j]);
                }
            }
        }
        writer.write(selected_matches);
    }
    reader.close();
    writer.close();
    //system.out
}

} /* namespace prot */
