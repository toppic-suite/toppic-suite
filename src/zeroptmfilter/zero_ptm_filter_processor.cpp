// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/file_util.hpp"
#include "base/web_logger.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"
#include "zeroptmfilter/mass_zero_ptm_filter.hpp"

namespace prot {

ZeroPtmFilterProcessor::ZeroPtmFilterProcessor(ZeroPtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
}

void ZeroPtmFilterProcessor::process(){
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  for(size_t i=0; i< db_block_ptr_vec.size(); i++){
    std::cout << "Zero PTM filtering - block " << (i+1) << " out of " 
        << db_block_ptr_vec.size() << " started." << std::endl; 
    processBlock(db_block_ptr_vec[i], db_block_ptr_vec.size());
    std::cout << "Zero PTM filtering - block " << (i +1) 
        << " finished. " << std::endl;
  }

  std::cout << "Zero PTM filtering - combining blocks started." << std::endl; 

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int block_num = db_block_ptr_vec.size();

  SimplePrsmStrCombine comp_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_COMPLETE",
                                    block_num, mng_ptr_->output_file_ext_ + "_COMPLETE", 
                                    mng_ptr_->comp_num_);
  comp_combine.process();

  SimplePrsmStrCombine pref_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_PREFIX",
                                    block_num, mng_ptr_->output_file_ext_ + "_PREFIX", 
                                    mng_ptr_->pref_suff_num_);
  pref_combine.process();

  SimplePrsmStrCombine suff_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_SUFFIX",
                                    block_num, mng_ptr_->output_file_ext_ + "_SUFFIX", 
                                    mng_ptr_->pref_suff_num_);
  suff_combine.process();

  SimplePrsmStrCombine internal_combine(sp_file_name, mng_ptr_->output_file_ext_ + "_INTERNAL",
                                        block_num, mng_ptr_->output_file_ext_ + "_INTERNAL", 
                                        mng_ptr_->inte_num_);
  internal_combine.process();

  std::cout << "Zero PTM filtering - combining blocks finished." << std::endl; 
}

void ZeroPtmFilterProcessor::processBlock(DbBlockPtr block_ptr, int total_block_num) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName() 
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms 
      = ProteoformFactory::readFastaToProteoformPtrVec(db_block_file_name, 
                                                       prsm_para_ptr->getFixModPtrVec());
  LOG_DEBUG("read fasta complete");
  MassZeroPtmFilterPtr filter_ptr(new MassZeroPtmFilter(raw_forms, mng_ptr_));
  LOG_DEBUG("filter inited");

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(),
                       group_spec_num,
                       prsm_para_ptr->getSpParaPtr()->getActivationPtr());
  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_;
  std::string block_str = std::to_string(block_ptr->getBlockIdx());

  SimplePrsmXmlWriter comp_writer(output_file_name + "_COMPLETE_" + block_str);
  SimplePrsmXmlWriter pref_writer(output_file_name + "_PREFIX_" + block_str);
  SimplePrsmXmlWriter suff_writer(output_file_name + "_SUFFIX_" + block_str);
  SimplePrsmXmlWriter internal_writer(output_file_name + "_INTERNAL_" + block_str);

  SpectrumSetPtr spec_set_ptr;
  int spectrum_num = MsAlignUtil::getSpNum (prsm_para_ptr->getSpectrumFileName());
  int cnt = 0;
  while((spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)) != nullptr){
    cnt+= group_spec_num;
    LOG_DEBUG("spec set ptr valid " << spec_set_ptr->isValid());
    if(spec_set_ptr->isValid()){
      ExtendMsPtrVec ms_ptr_vec = spec_set_ptr->getMsThreePtrVec();
      filter_ptr->computeBestMatch(ms_ptr_vec);
      comp_writer.write(filter_ptr->getCompMatchPtrs());
      pref_writer.write(filter_ptr->getPrefMatchPtrs());
      suff_writer.write(filter_ptr->getSuffMatchPtrs());
      internal_writer.write(filter_ptr->getInternalMatchPtrs());
    }
    WebLog::percentLog(cnt, spectrum_num, block_ptr->getBlockIdx(), total_block_num,  
                       WebLog::ZeroPtmFilterTime());
    std::cout << std::flush << "Zero PTM filtering - processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";
  }
  std::cout << std::endl;
  reader.close();
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
}

} /* namespace prot */
