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

#include <algorithm>
#include <cmath>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/base/prot_mod_util.hpp"
#include "seq/proteoform.hpp"
#include "seq/proteoform_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "search/varptmsearch/var_ptm_slow_match.hpp"
#include "search/varptmsearch/var_ptm_search_processor.hpp"

namespace toppic {

ProtModPtr findMatchedProtModPtr(ProteoformPtr db_proteo_ptr, 
                                 ProtModPtrVec &mod_ptr_vec, 
                                 double n_term_shift) {
  ResSeqPtr db_res_seq_ptr = db_proteo_ptr->getResSeqPtr();
  ResiduePtrVec residues = db_res_seq_ptr->getResidues();
  for (size_t i = 0; i < mod_ptr_vec.size(); i++) {
    double mod_shift = mod_ptr_vec[i]->getProtShift();
    // allow small errors introduced in double addition
    if (abs(mod_shift - n_term_shift) < 0.0001) {
      bool valid_mod = prot_mod_util::allowMod(mod_ptr_vec[i], residues); 
      if (valid_mod) {
        return mod_ptr_vec[i];
      }
    }
  }
  LOG_ERROR("Failed to find protein modification!" << n_term_shift);
  exit(EXIT_FAILURE);
}

int findMatchedPos(ProteoformPtr proteo_ptr, double n_term_trunc) {
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  for (size_t i = 0; i < prms.size(); i++) {
    // allow small errors introduced in double addition
    if (abs(prms[i] - n_term_trunc) < 0.0001) {
      return i;
    }
  }
  LOG_ERROR("Failed to find start posisiton!" << n_term_trunc);
  exit(EXIT_FAILURE);
}


std::function<void()> geneTask(SpectrumSetPtr spec_set_ptr,
                               SimplePrsmPtrVec simple_prsm_ptr_vec,
                               FastaIndexReaderPtr reader_ptr,
                               VarPtmSearchMngPtr mng_ptr,
                               ProteoformTypePtr type_ptr, 
                               SimpleThreadPoolPtr pool_ptr, 
                               PrsmXmlWriterPtrVec writer_ptr_vec) {
  return [spec_set_ptr, simple_prsm_ptr_vec, reader_ptr, mng_ptr, type_ptr,  pool_ptr, writer_ptr_vec]() {
    ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
    ProtModPtrVec prot_mod_ptr_vec = mng_ptr->prsm_para_ptr_->getProtModPtrVec();
    ProteoformPtrVec proteoform_ptr_vec;
    for (size_t i = 0; i < simple_prsm_ptr_vec.size(); i++) {
      if (std::abs(spec_set_ptr->getPrecMonoMass() - simple_prsm_ptr_vec[i]->getPrecMass()) 
          > std::pow(10, -4)) {
        // When precursor error is allowed, if the adjusted precursor of the
        // spectrum set does not match the adjusted precursor mass in the
        // filtering result, the spectrum proteoform match is ignored. 
        // A small error is allowed for errors introduced in writing real numbers
        // to xml files. 
        continue;
      }
      SimplePrsmPtr simple_prsm_ptr = simple_prsm_ptr_vec[i];
      std::string seq_name = simple_prsm_ptr->getSeqName();
      std::string seq_desc = simple_prsm_ptr->getSeqDesc();
      ProteoformPtr db_proteo_ptr
        = proteoform_factory::readFastaToProteoformPtr(reader_ptr, seq_name, seq_desc, fix_mod_list);

      double n_term_shift = simple_prsm_ptr->getNTermShifts()[0];
      double c_term_shift = simple_prsm_ptr->getCTermShifts()[0];
      LOG_DEBUG("type " << type_ptr->getName() << " proteoform mass " << db_proteo_ptr->getMass()
                << " n_term_shift " << n_term_shift << " c term shift " << c_term_shift);

      ProteoformPtr proteo_with_prot_mod_ptr = db_proteo_ptr;
      // start position is relative to the database sequence
      int start_pos = 0;
      if (type_ptr == ProteoformType::COMPLETE || type_ptr == ProteoformType::PREFIX) {
        if (n_term_shift != 0) {
          // Get the proteoform with N-terminal modificaiton
          ProtModPtr prot_mod_ptr = findMatchedProtModPtr(db_proteo_ptr, prot_mod_ptr_vec, n_term_shift);
          LOG_DEBUG("Db seq: " << db_proteo_ptr->getProteoformMatchSeq());
          proteo_with_prot_mod_ptr = proteoform_factory::geneProtModProteoform(db_proteo_ptr, prot_mod_ptr);
          if (proteo_with_prot_mod_ptr == nullptr) {
            LOG_ERROR("Null proteoform is reported!");
            exit(EXIT_FAILURE);
          }
          start_pos = proteo_with_prot_mod_ptr->getStartPos();
          LOG_DEBUG("Start pos with N terminal mod: " << start_pos);
        }
      }
      else {
        double n_term_trunc = - n_term_shift;
        start_pos = findMatchedPos(db_proteo_ptr, n_term_trunc);
      }
      int end_pos = db_proteo_ptr->getEndPos();
      if (type_ptr == ProteoformType::PREFIX || type_ptr == ProteoformType::INTERNAL) {
        double proteo_minus_water_mass = db_proteo_ptr->getMinusWaterMass();
        double n_term_trunc = proteo_minus_water_mass + c_term_shift;
        // the reported position is the position of matched breakpoint
        // the sequence end position is breakpoint position - 1
        end_pos = findMatchedPos(db_proteo_ptr, n_term_trunc) - 1; 
      }
      LOG_DEBUG("Start pos " << start_pos << " end position " << end_pos); 

      // get the start and end positions relative to the proteoform with
      // N-terminal modification
      int local_start_pos = start_pos - proteo_with_prot_mod_ptr->getStartPos();
      int local_end_pos = end_pos - proteo_with_prot_mod_ptr->getStartPos();
      FastaSeqPtr fasta_seq_ptr = db_proteo_ptr->getFastaSeqPtr();
      ProteoformPtr sub_proteo_ptr = proteoform_factory::geneSubProteoform(proteo_with_prot_mod_ptr,
                                                                           fasta_seq_ptr,
                                                                           local_start_pos, 
                                                                           local_end_pos);
      LOG_DEBUG("sub seq: " << sub_proteo_ptr->getProteoformMatchSeq());
      LOG_DEBUG("protein mass " << sub_proteo_ptr->getMass() 
                << "proteoform mass diff " << (simple_prsm_ptr->getPrecMass() - sub_proteo_ptr->getMass()));
      proteoform_ptr_vec.push_back(sub_proteo_ptr);
    }

    PrsmPtrVec prsms;
    for (size_t i = 0; i < proteoform_ptr_vec.size(); i++) {
      VarPtmSlowMatch slow_match(proteoform_ptr_vec[i], spec_set_ptr, mng_ptr);
      if (slow_match.isSuccessInit()) {
        PrsmPtr tmp = slow_match.compute();
        LOG_DEBUG("Slow match computation finish!");
        if (tmp != nullptr) {
          prsms.push_back(tmp);
        }
      }
    }
    std::sort(prsms.begin(), prsms.end(),
              Prsm::cmpMatchFragDecMatchPeakDecProtInc);
    if (prsms.size() > 0) {
      prsms.erase(prsms.begin() + mng_ptr->n_report_, prsms.end());
    }

    boost::thread::id thread_id = boost::this_thread::get_id();
    int writer_id = pool_ptr->getId(thread_id);
    writer_ptr_vec[writer_id]->writeVector(prsms);
  };
}

void VarPtmSearchProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  SimplePrsmReader comp_prsm_reader(input_file_name + "_" + ProteoformType::COMPLETE->getName());
  SimplePrsmReader pref_prsm_reader(input_file_name + "_" + ProteoformType::PREFIX->getName());
  SimplePrsmReader suff_prsm_reader(input_file_name + "_" + ProteoformType::SUFFIX->getName());
  SimplePrsmReader internal_prsm_reader(input_file_name + "_" + ProteoformType::INTERNAL->getName());
  SimplePrsmPtr comp_prsm_ptr = comp_prsm_reader.readOnePrsm();
  SimplePrsmPtr pref_prsm_ptr = pref_prsm_reader.readOnePrsm();
  SimplePrsmPtr suff_prsm_ptr = suff_prsm_reader.readOnePrsm();
  SimplePrsmPtr internal_prsm_ptr = internal_prsm_reader.readOnePrsm();

  std::string output_file_name = file_util::basename(prsm_para_ptr->getSpectrumFileName())
      + "." + mng_ptr_->output_file_ext_;

  PrsmXmlWriterPtrVec comp_writer_ptr_vec;
  PrsmXmlWriterPtrVec pref_writer_ptr_vec;
  PrsmXmlWriterPtrVec suff_writer_ptr_vec;
  PrsmXmlWriterPtrVec internal_writer_ptr_vec;
  for (int i = 0; i < mng_ptr_->thread_num_; i++) { 
    std::string idx = str_util::toString(i);
    PrsmXmlWriterPtr comp_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name 
                                                                   + "_" + ProteoformType::COMPLETE->getName() 
                                                                   + "_" + idx);
    comp_writer_ptr_vec.push_back(comp_writer_ptr);
    PrsmXmlWriterPtr pref_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name 
                                                                       + "_" + ProteoformType::PREFIX->getName()
                                                                       + "_" + idx);
    pref_writer_ptr_vec.push_back(pref_writer_ptr);

    PrsmXmlWriterPtr suff_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name 
                                                                       + "_" + ProteoformType::SUFFIX->getName()
                                                                       + "_" + idx);
    suff_writer_ptr_vec.push_back(suff_writer_ptr);

    PrsmXmlWriterPtr internal_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name 
                                                                           + "_" + ProteoformType::INTERNAL->getName()
                                                                           + "_" + idx);
    internal_writer_ptr_vec.push_back(internal_writer_ptr);
  }


  SimpleThreadPoolPtr pool_ptr = std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);

  // init variables
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();

  FastaIndexReaderPtr reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);
  int spectrum_num = msalign_util::getSpNum(sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr msalign_reader_ptr = std::make_shared<MsAlignReader>(sp_file_name, 
                                                                        group_spec_num,
                                                                        sp_para_ptr->getActivationPtr());
  int cnt = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = msalign_reader_ptr->getNextMsPtrVec(); 
  std::vector<double> prec_error_vec = sp_para_ptr->getVarPtmSearchPrecErrorVec();
  //LOG_ERROR("start ");
  while (deconv_ms_ptr_vec.size() > 0) {
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      std::vector<SpectrumSetPtr> spec_set_vec 
        = spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError(deconv_ms_ptr_vec, 
                                                                   sp_para_ptr,
                                                                   prec_error_vec);
      if (spec_set_vec.size() == 0) {
        LOG_ERROR("Spectrum set size is 0!");
      }

      if (spec_set_vec[0]->isValid()) {
        int spec_id = spec_set_vec[0]->getSpectrumId();
        // complete
        SimplePrsmPtrVec comp_selected_prsm_ptrs;
        while (comp_prsm_ptr != nullptr && comp_prsm_ptr->getSpectrumId() == spec_id) {
          comp_selected_prsm_ptrs.push_back(comp_prsm_ptr);
          comp_prsm_ptr = comp_prsm_reader.readOnePrsm();
        }
        //LOG_ERROR("complete list size " << comp_selected_prsm_ptrs.size());
        if (comp_selected_prsm_ptrs.size() > 0) {
          // LOG_DEBUG("start processing one spectrum.");
          for (size_t k = 0; k < spec_set_vec.size(); k++) {
            while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
              boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(spec_set_vec[k], 
                                       comp_selected_prsm_ptrs,
                                       reader_ptr, 
                                       mng_ptr_,
                                       ProteoformType::COMPLETE,
                                       pool_ptr,
                                       comp_writer_ptr_vec));
          }
        }

        // prefix
        SimplePrsmPtrVec pref_selected_prsm_ptrs;
        while (pref_prsm_ptr != nullptr && pref_prsm_ptr->getSpectrumId() == spec_id) {
          pref_selected_prsm_ptrs.push_back(pref_prsm_ptr);
          pref_prsm_ptr = pref_prsm_reader.readOnePrsm();
        }
        LOG_DEBUG("prefix list size " << pref_selected_prsm_ptrs.size());
        if (pref_selected_prsm_ptrs.size() > 0) {
          // LOG_DEBUG("start processing one spectrum.");
          for (size_t k = 0; k < spec_set_vec.size(); k++) {
            while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
              boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(spec_set_vec[k], 
                                       pref_selected_prsm_ptrs,
                                       reader_ptr, 
                                       mng_ptr_,
                                       ProteoformType::PREFIX,
                                       pool_ptr,
                                       pref_writer_ptr_vec));
          }
        }

        // suffix
        SimplePrsmPtrVec suff_selected_prsm_ptrs;
        while (suff_prsm_ptr != nullptr && suff_prsm_ptr->getSpectrumId() == spec_id) {
          suff_selected_prsm_ptrs.push_back(suff_prsm_ptr);
          suff_prsm_ptr = suff_prsm_reader.readOnePrsm();
        }
        LOG_DEBUG("suffix list size " << suff_selected_prsm_ptrs.size());
        if (suff_selected_prsm_ptrs.size() > 0) {
          // LOG_DEBUG("start processing one spectrum.");
          for (size_t k = 0; k < spec_set_vec.size(); k++) {
            while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
              boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(spec_set_vec[k], 
                                       suff_selected_prsm_ptrs,
                                       reader_ptr, 
                                       mng_ptr_,
                                       ProteoformType::SUFFIX,
                                       pool_ptr,
                                       suff_writer_ptr_vec));
          }
        }

        // internal
        SimplePrsmPtrVec internal_selected_prsm_ptrs;
        while (internal_prsm_ptr != nullptr && internal_prsm_ptr->getSpectrumId() == spec_id) {
          internal_selected_prsm_ptrs.push_back(internal_prsm_ptr);
          internal_prsm_ptr = internal_prsm_reader.readOnePrsm();
        }
        LOG_DEBUG("internal list size " << internal_selected_prsm_ptrs.size());
        if (internal_selected_prsm_ptrs.size() > 0) {
          // LOG_DEBUG("start processing one spectrum.");
          for (size_t k = 0; k < spec_set_vec.size(); k++) {
            while (pool_ptr->getQueueSize() >= mng_ptr_->thread_num_ * 2) {
              boost::this_thread::sleep(boost::posix_time::milliseconds(100));
            }
            pool_ptr->Enqueue(geneTask(spec_set_vec[k], 
                                       internal_selected_prsm_ptrs,
                                       reader_ptr, 
                                       mng_ptr_,
                                       ProteoformType::INTERNAL,
                                       pool_ptr,
                                       internal_writer_ptr_vec));
          }
        }
      }
    }
    cnt = cnt + group_spec_num; 
    std::cout << std::flush <<  "Variable PTM search - processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
    deconv_ms_ptr_vec = msalign_reader_ptr->getNextMsPtrVec(); 
  }
  pool_ptr->ShutDown();
  int remainder = spectrum_num - cnt;
  if (prsm_para_ptr->getGroupSpecNum() > remainder && remainder > 0){
      //if there are spectrum remaining because they were not combined due to not having enough pairs
      //fix the message as the processing is completed.
      //this code avoids error when no combined spectra is used but a scan is remaining unprocessed
      //because then it will not satisfy the first condition
      std::cout << std::flush <<  "Variable PTM search - processing " << spectrum_num
        << " of " << spectrum_num << " spectra.\r";
  } 
  comp_prsm_reader.close();
  pref_prsm_reader.close();
  suff_prsm_reader.close();
  internal_prsm_reader.close();
  for (int i = 0; i < mng_ptr_->thread_num_; i++) { 
    comp_writer_ptr_vec[i]->close();
    pref_writer_ptr_vec[i]->close();
    suff_writer_ptr_vec[i]->close();
    internal_writer_ptr_vec[i]->close();
  }
  std::cout << std::endl;

  bool norm = false;
  bool remove_dup = true;
  PrsmStrMergePtr merge_ptr
    = std::make_shared<PrsmStrMerge>(sp_file_name, 
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::COMPLETE->getName(),
                                     mng_ptr_->thread_num_,
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::COMPLETE->getName(),
                                     mng_ptr_->n_report_, norm, remove_dup);
  merge_ptr->process();
  merge_ptr = nullptr;

  merge_ptr
    = std::make_shared<PrsmStrMerge>(sp_file_name, 
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::PREFIX->getName(),
                                     mng_ptr_->thread_num_,
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::PREFIX->getName(),
                                     mng_ptr_->n_report_, norm, remove_dup);
  merge_ptr->process();
  merge_ptr = nullptr; 

  merge_ptr
    = std::make_shared<PrsmStrMerge>(sp_file_name, 
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::SUFFIX->getName(),
                                     mng_ptr_->thread_num_,
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::SUFFIX->getName(),
                                     mng_ptr_->n_report_, norm, remove_dup);
  merge_ptr->process();
  merge_ptr = nullptr; 

  merge_ptr
    = std::make_shared<PrsmStrMerge>(sp_file_name, 
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::INTERNAL->getName(),
                                     mng_ptr_->thread_num_,
                                     mng_ptr_->output_file_ext_ + "_" +
                                     ProteoformType::INTERNAL->getName(),
                                     mng_ptr_->n_report_, norm, remove_dup);
  merge_ptr->process();
  merge_ptr = nullptr; 

  //remove temporary files
  file_util::cleanTempFiles(sp_file_name, 
                            mng_ptr_->output_file_ext_ + "_" 
                            + ProteoformType::COMPLETE->getName() + "_");
  file_util::cleanTempFiles(sp_file_name, 
                            mng_ptr_->output_file_ext_ + "_" +
                            ProteoformType::PREFIX->getName() + "_");
  file_util::cleanTempFiles(sp_file_name, 
                            mng_ptr_->output_file_ext_ + "_" +
                            ProteoformType::SUFFIX->getName() + "_");
  file_util::cleanTempFiles(sp_file_name, 
                            mng_ptr_->output_file_ext_ + "_" +
                            ProteoformType::INTERNAL->getName() + "_");
}

}  // namespace toppic
