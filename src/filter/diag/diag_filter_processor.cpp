// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "filter/diag/diag_filter_processor.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "common/base/mod_util.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "filter/diag/diag_filter.hpp"
#include "ms/factory/prm_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_xml_writer_util.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace {

// The text record of a candidate as SimplePrsmReader reads it back from a
// file written by SimplePrsmXmlWriter: the serialized element, split into
// trimmed non-empty lines.
SimplePrsmStrPtr toPrsmStr(const SimplePrsmPtr& prsm_ptr) {
  XmlDOMDocument doc("simple_prsm_list");
  XmlDOMElement root = doc.getDocumentElement();
  XmlDOMElement element = prsm_ptr->toXml(&doc, root);
  std::string str = xml_dom_util::writeToString(element);
  std::vector<std::string> lines;
  std::string line;
  std::istringstream stream(str);
  while (std::getline(stream, line)) {
    str_util::trim(line);
    if (line != "") {
      lines.push_back(line);
    }
  }
  return std::make_shared<SimplePrsmStr>(lines);
}

}  // namespace

inline void filterBlock(const ProteoformPtrVec& raw_forms, int block_idx,
                        const DiagFilterMngPtr& mng_ptr,
                        const std::vector<double>& mod_mass_list,
                        DiagFilterProcessor* processor) {
  std::string block_str = std::to_string(block_idx);

  DiagFilterPtr filter_ptr =
      std::make_shared<DiagFilter>(raw_forms, mng_ptr, block_str);

  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  int group_spec_num = mng_ptr->prsm_para_ptr_->getGroupSpecNum();

  MsAlignReaderPtr reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr(),
      sp_para_ptr->getEnvCnnCutoff());

  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  std::vector<double> prec_error_vec =
      sp_para_ptr->getMultiShiftSearchPrecErrorVec();
  while (deconv_ms_ptr_vec.size() != 0) {
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      // allow one dalton error
      SpectrumSetPtrVec spec_set_vec =
          spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError(
              deconv_ms_ptr_vec, sp_para_ptr, prec_error_vec);
      for (size_t k = 0; k < spec_set_vec.size(); k++) {
        SpectrumSetPtr spec_set_ptr = spec_set_vec[k];
        if (spec_set_ptr->isValid()) {
          int spec_id = spec_set_ptr->getSpectrumId();
          if (mng_ptr->var_num_ == 0) {
            PrmMsPtrVec ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
            SimplePrsmPtrVec match_ptrs = filter_ptr->getBestMatch(ms_ptr_vec);
            processor->addResults(spec_id, match_ptrs);
          } else {
            DeconvMsPtrVec deconv_ms_ptr_vec =
                spec_set_ptr->getDeconvMsPtrVec();
            double prec_mono_mass = spec_set_ptr->getPrecMonoMass();
            double n_term_label_mass = spec_set_ptr->getNTermLabelMass();
            std::vector<double> mod_mass(3);
            for (size_t i = 0; i < mod_mass_list.size(); i++) {
              for (size_t k1 = 0; k1 < mod_mass.size(); k1++) {
                std::fill(mod_mass.begin(), mod_mass.end(), 0.0);
                mod_mass[k1] += mod_mass_list[i];
                PrmMsPtrVec ms_ptr_vec = prm_ms_factory::geneMsTwoPtrVec(
                    deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass,
                    n_term_label_mass, mod_mass);
                SimplePrsmPtrVec match_ptrs =
                    filter_ptr->getBestMatch(ms_ptr_vec);
                processor->addResults(spec_id, match_ptrs);
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
    msg << std::flush << "Multiple unexpected shifts filtering - processing "
        << std::setprecision(3) << perc << "%.     \r";
    mng_ptr->mutex_.lock();
    std::cout << msg.str();
    mng_ptr->mutex_.unlock();

    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
}

std::function<void()> geneTask(int block_idx,
                               const std::vector<double>& mod_mass_list,
                               const DiagFilterMngPtr& mng_ptr,
                               DiagFilterProcessor* processor) {
  return [block_idx, mod_mass_list, mng_ptr, processor]() {
    PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
    std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
    std::string db_block_file_name =
        prsm_para_ptr->getSearchDbFileNameWithFolder() + "_" +
        std::to_string(block_idx);
    ProteoformPtrVec raw_forms =
        proteoform_factory::readFastaToProteoformPtrVec(
            db_block_file_name, prsm_para_ptr->getFixModPtrVec());
    filterBlock(raw_forms, block_idx, mng_ptr, mod_mass_list, processor);
  };
}

// The candidates of a spectrum are combined in memory instead of in one file
// per block: the running list of a spectrum is cut back to the top_num best
// (score decreasing, sequence name increasing, the order of the final file)
// whenever it grows past twice that number, which cannot drop a record of
// the final top list.
void DiagFilterProcessor::addResults(int spec_id,
                                     const SimplePrsmPtrVec& match_ptrs) {
  if (match_ptrs.size() == 0) {
    return;
  }
  SimplePrsmStrPtrVec str_ptrs;
  for (size_t i = 0; i < match_ptrs.size(); i++) {
    str_ptrs.push_back(toPrsmStr(match_ptrs[i]));
  }
  size_t top_num = mng_ptr_->filter_result_num_;
  std::lock_guard<std::mutex> lock(
      spec_mutexes_[spec_id % spec_mutexes_.size()]);
  if (static_cast<size_t>(spec_id) >= spec_results_.size()) {
    std::lock_guard<std::mutex> resize_lock(resize_mutex_);
    spec_results_.resize(spec_id + 1);
  }
  SimplePrsmStrPtrVec& results = spec_results_[spec_id];
  results.insert(results.end(), str_ptrs.begin(), str_ptrs.end());
  if (results.size() > 2 * top_num) {
    std::sort(results.begin(), results.end(),
              SimplePrsmStr::cmpScoreDecSeqInc);
    results.resize(top_num);
  }
}

void DiagFilterProcessor::writeResults(const std::string& output_file_name) {
  SimplePrsmXmlWriter writer(output_file_name);
  size_t top_num = mng_ptr_->filter_result_num_;
  for (size_t spec_id = 0; spec_id < spec_results_.size(); spec_id++) {
    SimplePrsmStrPtrVec& results = spec_results_[spec_id];
    std::sort(results.begin(), results.end(),
              SimplePrsmStr::cmpScoreDecSeqInc);
    for (size_t i = 0; i < results.size() && i < top_num; i++) {
      writer.write(results[i]);
    }
    results.clear();
  }
  writer.close();
}

void DiagFilterProcessor::process() {
  std::string db_file_name =
      mng_ptr_->prsm_para_ptr_->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);

  std::vector<double> mod_mass_list;
  if (mng_ptr_->residue_mod_file_name_ != "") {
    mod_mass_list = mod_util::getModMassVec(
        mod_util::readModTxt(mng_ptr_->residue_mod_file_name_)[2]);
  }

  int spec_num =
      msalign_util::getSpNum(mng_ptr_->prsm_para_ptr_->getSpectrumFileName());
  mng_ptr_->n_spec_block_ = spec_num * db_block_ptr_vec.size();
  spec_results_.resize(spec_num);
  std::vector<std::mutex> mutexes(64);
  spec_mutexes_.swap(mutexes);

  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);
  int block_num = db_block_ptr_vec.size();
  mng_ptr_->cnts_.resize(block_num, 0);

  for (size_t i = 0; i < db_block_ptr_vec.size(); i++) {
    while (pool_ptr->getQueueSize() > 0 || pool_ptr->getIdleThreadNum() == 0) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    pool_ptr->enqueue(geneTask(db_block_ptr_vec[i]->getBlockIdx(),
                               mod_mass_list, mng_ptr_, this));
  }
  pool_ptr->shutDown();
  std::cout << std::endl;

  std::cout
      << "Multiple unexpected shifts filtering - combining blocks started."
      << std::endl;
  std::string sp_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  writeResults(file_util::basename(sp_file_name) + "." +
               mng_ptr_->output_file_ext_);
  std::cout
      << "Multiple unexpected shifts filtering - combining blocks finished."
      << std::endl;
}

}  // namespace toppic
