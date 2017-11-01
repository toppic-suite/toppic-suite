//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include <string>
#include <vector>

#include "base/proteoform.hpp"
#include "base/proteoform_factory.hpp"
#include "base/file_util.hpp"
#include "base/mod_util.hpp"
#include "spec/msalign_util.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"
#include "oneptmfilter/mass_one_ptm_filter.hpp"

namespace prot {

void OnePtmFilterProcessor::process() {
  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residueModFileName_ != "") {
    mod_mass_list = ModUtil::getModMassVec(ModUtil::readModTxt(mng_ptr_->residueModFileName_)[2]);
  }

  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    std::cout << "One PTM filtering - block " << (i + 1) << " out of "
        << db_block_ptr_vec.size() << " started." << std::endl;
    processBlock(db_block_ptr_vec[i], db_block_ptr_vec.size(), mod_mass_list);
    std::cout << "One PTM filtering - block " << (i + 1) << " finished. " << std::endl;
  }

  std::cout << "One PTM filtering - combining blocks started." << std::endl;

  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  int block_num = db_block_ptr_vec.size();

  LOG_DEBUG("comp number " << mng_ptr_->comp_num_);
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

  std::cout << "One PTM filtering - combining blocks finished." << std::endl;
}

void OnePtmFilterProcessor::processBlock(DbBlockPtr block_ptr, int total_block_num,
                                         const std::vector<double> & mod_mass_list) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string db_block_file_name = prsm_para_ptr->getSearchDbFileName()
      + "_" + std::to_string(block_ptr->getBlockIdx());
  ProteoformPtrVec raw_forms
      = ProteoformFactory::readFastaToProteoformPtrVec(db_block_file_name,
                                                       prsm_para_ptr->getFixModPtrVec());
  MassOnePtmFilterPtr filter_ptr = std::make_shared<MassOnePtmFilter>(raw_forms, mng_ptr_);

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr =  mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  sp_para_ptr->prec_error_ = 0;
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(), group_spec_num,
                       sp_para_ptr->getActivationPtr(),
                       sp_para_ptr->getSkipList());
  std::string output_file_name = FileUtil::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_;
  std::string block_str = std::to_string(block_ptr->getBlockIdx());

  SimplePrsmXmlWriter comp_writer(output_file_name + "_COMPLETE_" + block_str);
  SimplePrsmXmlWriter pref_writer(output_file_name + "_PREFIX_" + block_str);
  SimplePrsmXmlWriter suff_writer(output_file_name + "_SUFFIX_" + block_str);
  SimplePrsmXmlWriter internal_writer(output_file_name + "_INTERNAL_" + block_str);

  SpectrumSetPtr spec_set_ptr;
  int spectrum_num = MsAlignUtil::getSpNum(prsm_para_ptr->getSpectrumFileName());
  int cnt = 0;
  while ((spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr)[0]) != nullptr) {
    cnt+= group_spec_num;
    if (spec_set_ptr->isValid()) {
      if (mng_ptr_->var_num_ == 0) {
        PrmMsPtrVec prm_ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
        PrmMsPtrVec srm_ms_ptr_vec = spec_set_ptr->getSuffixMsTwoPtrVec();
        filter_ptr->computeBestMatch(prm_ms_ptr_vec, srm_ms_ptr_vec);
        comp_writer.write(filter_ptr->getCompMatchPtrs());
        pref_writer.write(filter_ptr->getPrefMatchPtrs());
        suff_writer.write(filter_ptr->getSuffMatchPtrs());
        internal_writer.write(filter_ptr->getInternalMatchPtrs());
      } else {
        for (size_t i = 0; i < mod_mass_list.size(); i++) {
          for (size_t k1 = 0; k1 < sp_para_ptr->mod_mass_.size(); k1++) {
            std::fill(sp_para_ptr->mod_mass_.begin(), sp_para_ptr->mod_mass_.end(), 0.0);
            sp_para_ptr->mod_mass_[k1] += mod_mass_list[i];
            PrmMsPtrVec prm_ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec(sp_para_ptr);
            PrmMsPtrVec srm_ms_ptr_vec = spec_set_ptr->getSuffixMsTwoPtrVec(sp_para_ptr);
            filter_ptr->computeBestMatch(prm_ms_ptr_vec, srm_ms_ptr_vec);
            comp_writer.write(filter_ptr->getCompMatchPtrs());
            pref_writer.write(filter_ptr->getPrefMatchPtrs());
            suff_writer.write(filter_ptr->getSuffMatchPtrs());
            internal_writer.write(filter_ptr->getInternalMatchPtrs());
          }
        }
      }
    }
    std::cout << std::flush << "One PTM filtering - processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
  }
  std::cout << std::endl;
  reader.close();
  comp_writer.close();
  pref_writer.close();
  suff_writer.close();
  internal_writer.close();
  LOG_DEBUG("block end");
}

} /* namespace prot */
