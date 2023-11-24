//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <iomanip>
#include <iostream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/base/mod_util.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/factory/prm_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"

#include "prsm/simple_prsm_xml_writer_set.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/oneptm/one_ptm_filter.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"

namespace toppic {

namespace one_ptm_filter_processor {

inline void filterBlock(const ProteoformPtrVec & raw_forms,
                        int block_idx, 
                        OnePtmFilterMngPtr mng_ptr,
                        const std::vector<double> & mod_mass_list) {
  std::string block_str = str_util::toString(block_idx);

  OnePtmFilterPtr filter_ptr = std::make_shared<OnePtmFilter>(raw_forms, mng_ptr, block_str);
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr reader_ptr = std::make_shared<MsAlignReader>(sp_file_name, 
                                                                group_spec_num,
                                                                sp_para_ptr->getActivationPtr());
  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr->output_file_ext_ + "_" + block_str;
  SimplePrsmXmlWriterSet writers(output_file_name);

  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  std::vector<double> prec_error_vec = sp_para_ptr->getOneShiftSearchPrecErrorVec(); 
  while(deconv_ms_ptr_vec.size() != 0) { 
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      // allow one dalton error
      SpectrumSetPtrVec spec_set_vec 
        = spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError(deconv_ms_ptr_vec, 
                                                                   sp_para_ptr,
                                                                   prec_error_vec);
      for (size_t k = 0; k < spec_set_vec.size(); k++) {
        SpectrumSetPtr spec_set_ptr = spec_set_vec[k];
        if (spec_set_ptr->isValid()) {
          if (mng_ptr->var_num_ == 0) {
            PrmMsPtrVec prm_ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
            PrmMsPtrVec srm_ms_ptr_vec = spec_set_ptr->getSuffixMsTwoPtrVec();
            filter_ptr->computeBestMatch(prm_ms_ptr_vec, srm_ms_ptr_vec);
            writers.getCompleteWriterPtr()->write(filter_ptr->getCompMatchPtrs());
            writers.getPrefixWriterPtr()->write(filter_ptr->getPrefMatchPtrs());
            writers.getSuffixWriterPtr()->write(filter_ptr->getSuffMatchPtrs());
            writers.getInternalWriterPtr()->write(filter_ptr->getInternalMatchPtrs());
          } 
          else {
            std::vector<double> mod_mass(3);
            for (size_t i = 0; i < mod_mass_list.size(); i++) {
              for (size_t k1 = 0; k1 < mod_mass.size(); k1++) {
                std::fill(mod_mass.begin(), mod_mass.end(), 0.0);
                mod_mass[k1] += mod_mass_list[i];
                DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
                double prec_mono_mass = spec_set_ptr->getPrecMonoMass();
                double n_term_label_mass = spec_set_ptr->getNTermLabelMass();
                PrmMsPtrVec prm_ms_ptr_vec = prm_ms_factory::geneMsTwoPtrVec(deconv_ms_ptr_vec,
                                                                             sp_para_ptr,
                                                                             prec_mono_mass, 
                                                                             n_term_label_mass,
                                                                             mod_mass);
                PrmMsPtrVec srm_ms_ptr_vec 
                  = prm_ms_factory::geneSuffixMsTwoPtrVec(deconv_ms_ptr_vec, sp_para_ptr,
                                                          prec_mono_mass, n_term_label_mass, mod_mass);
                filter_ptr->computeBestMatch(prm_ms_ptr_vec, srm_ms_ptr_vec);
                writers.getCompleteWriterPtr()->write(filter_ptr->getCompMatchPtrs());
                writers.getPrefixWriterPtr()->write(filter_ptr->getPrefMatchPtrs());
                writers.getSuffixWriterPtr()->write(filter_ptr->getSuffMatchPtrs());
                writers.getInternalWriterPtr()->write(filter_ptr->getInternalMatchPtrs());
              }
            }
          }
        }
      }
    }

    mng_ptr->cnts_[block_idx] = mng_ptr->cnts_[block_idx] + group_spec_num;
    int cnt_sum = 0; 
    for (size_t i = 0; i < mng_ptr->cnts_.size(); i++) {
      cnt_sum = cnt_sum + mng_ptr->cnts_[i];
    }
    double perc = cnt_sum * 100.0 / mng_ptr->n_spec_block_;
    std::stringstream msg;
    msg << std::flush << "One unexpected shift filtering - processing " << std::setprecision(3) <<  perc << "%.     \r";
    mng_ptr->mutex_.lock();
    std::cout << msg.str();
    mng_ptr->mutex_.unlock();
    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
  writers.close();
}

std::function<void()> geneTask(int block_idx, 
                               const std::vector<double> &mod_mass_list, 
                               OnePtmFilterMngPtr mng_ptr) {
  return[block_idx, mod_mass_list, mng_ptr] () {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder()
      + "_" + str_util::toString(block_idx);
    ProteoformPtrVec raw_forms
        = proteoform_factory::readFastaToProteoformPtrVec(db_block_file_name,
                                                          prsm_para_ptr->getFixModPtrVec());
    filterBlock(raw_forms, block_idx, mng_ptr, mod_mass_list);
  };
}

void process(OnePtmFilterMngPtr mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::vector<double> mod_mass_list;
  if (mng_ptr->residueModFileName_ != "") {
    mod_mass_list = mod_util::getModMassVec(mod_util::readModTxt(mng_ptr->residueModFileName_)[2]);
  }
  int spec_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  mng_ptr->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);
  int block_num = db_block_ptr_vec.size();
  mng_ptr->cnts_.resize(block_num, 0);
  for (int i = 0; i < block_num; i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      boost::this_thread::sleep(boost::posix_time::milliseconds(100));
    }
    pool_ptr->Enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(), mod_mass_list, mng_ptr));
  }
  pool_ptr->ShutDown();
  std::cout << std::endl;

  std::cout << "One unexpected shift filtering - combining blocks started." << std::endl;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_pref = mng_ptr->output_file_ext_;
  SimplePrsmStrMerge::mergeBlockResults(sp_file_name, input_pref, block_num,  
                                        mng_ptr->comp_num_, mng_ptr->pref_suff_num_, 
                                        mng_ptr->inte_num_);
  std::cout << "One unexpected shift filtering - combining blocks finished." << std::endl;
}

}

} /* namespace toppic */
