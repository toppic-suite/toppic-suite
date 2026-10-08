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

#include "stat/tdgf/evalue_processor.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <mutex>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/prsm_xml_writer_util.hpp"

namespace toppic {

namespace {
// CountTestNum is a deterministic pure function of these parameters (a full
// target+decoy database scan, about 9 s on the human database with a decoy),
// and one toppic -s 2 run used to build the identical object three times:
// the first-pass E-value computation, the multiple-shift E-value computation
// and the cluster protein correction. It is immutable after construction and
// already shared read-only across the worker threads of one stage, so one
// instance serves all three consumers; reuse changes no arithmetic, order or
// tie-break (toppic_test experiment 66). The key holds every input the
// constructor reads: its three numeric parameters and, from PrsmPara, the
// search database file and the fixed/N-terminal modification lists (each
// pipeline step builds its own PrsmPara object from the same arguments, so
// the key must compare values, not the PrsmPara pointer).
std::mutex count_test_num_cache_mutex;
std::map<std::string, CountTestNumPtr> count_test_num_cache;

std::string countTestNumKey(const TdgfMngPtr& mng_ptr) {
  const PrsmParaPtr& para_ptr = mng_ptr->prsm_para_ptr_;
  std::string key = std::to_string(mng_ptr->convert_ratio_) + "|" +
                    std::to_string(mng_ptr->max_ptm_mass_) + "|" +
                    std::to_string(mng_ptr->max_prec_mass_) + "|" +
                    para_ptr->getSearchDbFileNameWithFolder() + "|";
  for (const ModPtr& mod_ptr : para_ptr->getFixModPtrVec()) {
    key += mod_ptr->getOriResiduePtr()->toString() + ">" +
           mod_ptr->getModResiduePtr()->toString() + ",";
  }
  key += "|";
  for (const ProtModPtr& prot_mod_ptr : para_ptr->getProtModPtrVec()) {
    key += prot_mod_ptr->getName() + ",";
  }
  return key;
}
}  // namespace

void EValueProcessor::init() {
  const std::string key = countTestNumKey(mng_ptr_);
  {
    std::lock_guard<std::mutex> lock(count_test_num_cache_mutex);
    auto it = count_test_num_cache.find(key);
    if (it != count_test_num_cache.end()) {
      test_num_ptr_ = it->second;
    }
  }
  if (test_num_ptr_ == nullptr) {
    test_num_ptr_ = std::make_shared<CountTestNum>(
        mng_ptr_->convert_ratio_, mng_ptr_->max_ptm_mass_,
        mng_ptr_->max_prec_mass_, mng_ptr_->prsm_para_ptr_);
    std::lock_guard<std::mutex> lock(count_test_num_cache_mutex);
    count_test_num_cache[key] = test_num_ptr_;
  }

  ResFreqPtrVec residue_freqs = test_num_ptr_->getResFreqPtrVec();
  comp_pvalue_ptr_ = std::make_shared<CompPValueArray>(test_num_ptr_, mng_ptr_);
}

std::function<void()> geneTask(const SpectrumSetPtr& spec_set_ptr,
                               const PrsmPtrVec& sele_prsm_ptrs,
                               bool is_separate, const TdgfMngPtr& mng_ptr,
                               const CountTestNumPtr& test_num_ptr,
                               const SimpleThreadPoolPtr& pool_ptr,
                               const PrsmXmlWriterPtrVec& writer_ptr_vec) {
  return [spec_set_ptr, sele_prsm_ptrs, is_separate, mng_ptr, test_num_ptr,
          pool_ptr, writer_ptr_vec]() {
    PrsmPtrVec prsm_vec;  // copy sele_prsm_ptrs
    for (size_t i = 0; i < sele_prsm_ptrs.size(); i++) {
      prsm_vec.push_back(std::make_shared<Prsm>(*sele_prsm_ptrs[i].get()));
    }
    // one CompPValueArray per worker thread: its CompProbValue allocates
    // a 138 MB page table and two 5.5 MB score arrays, which were allocated
    // and freed for every spectrum
    static thread_local CompPValueArrayPtr comp_pvalue_ptr;
    if (comp_pvalue_ptr == nullptr) {
      comp_pvalue_ptr = std::make_shared<CompPValueArray>(test_num_ptr,
                                                          mng_ptr);
    }
    // A PrSM with at most comp_evalue_min_match_theo_frag_num_ matched
    // fragments gets the maximum E-value below whatever the generating
    // function says, so a spectrum with no other PrSM needs no generating
    // function (toppic_test experiment 150).
    bool evalue_needed = false;
    for (size_t i = 0; i < prsm_vec.size(); i++) {
      if (std::round(prsm_vec[i]->getMatchTheoFragNum()) >
          std::round(mng_ptr->comp_evalue_min_match_theo_frag_num_)) {
        evalue_needed = true;
        break;
      }
    }
    if (evalue_needed) {
      comp_pvalue_ptr->process(spec_set_ptr, prsm_vec, is_separate);
    }

    for (unsigned i = 0; i < prsm_vec.size(); i++) {
      if (std::round(prsm_vec[i]->getMatchTheoFragNum()) <=
          std::round(mng_ptr->comp_evalue_min_match_theo_frag_num_)) {
        prsm_vec[i]->setExpectedValuePtr(ExpectedValue::getMaxEvaluePtr());
      } else {
        if (prsm_vec[i]->getEValue() == 0.0) {
          LOG_WARN("Invalid e value!");
          prsm_vec[i]->setExpectedValuePtr(ExpectedValue::getMaxEvaluePtr());
        }
      }
    }

    std::thread::id thread_id = std::this_thread::get_id();
    int writer_id = pool_ptr->getId(thread_id);
    for (size_t i = 0; i < prsm_vec.size(); i++) {
      writer_ptr_vec[writer_id]->write(prsm_vec[i]);
    }
  };
}

// Compute E-value. Separate: compute E-values for PrSMs separately or not
void EValueProcessor::process(bool is_separate) {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string output_file_name =
      file_util::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;
  PrsmXmlWriter writer(output_file_name);

  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  FastaIndexReaderPtr seq_reader =
      std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
  std::string input_file_name =
      file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  // init variables
  int spectrum_num = msalign_util::getSpNum(sp_file_name);
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr(),
      sp_para_ptr->getEnvCnnCutoff(),
        sp_para_ptr->getFragDedupPpo());

  PrsmXmlWriterPtrVec writer_ptr_vec = prsm_xml_writer_util::geneWriterPtrVec(
      output_file_name, mng_ptr_->thread_num_);

  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(mng_ptr_->thread_num_);

  int cnt = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() != 0) {
    cnt += group_spec_num;
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
    if (header_ptr->containsPrec()) {
      // the observed precursor mass: the spectrum set subtracts the N-terminal
      // label (toppic -N) itself, for its PRM and extended spectra
      double prec_mono_mass = header_ptr->getFirstPrecMonoMass();
      SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
          deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        PrsmPtrVec selected_prsm_ptrs;
        while (prsm_ptr != nullptr &&
               prsm_ptr->getSpectrumId() == spec_set_ptr->getSpectrumId()) {
          selected_prsm_ptrs.push_back(prsm_ptr);
          prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
        }
        if (checkPrsms(selected_prsm_ptrs)) {
          while (pool_ptr->getQueueSize() >=
                 static_cast<std::size_t>(mng_ptr_->thread_num_) + 2) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
          }
          pool_ptr->enqueue(geneTask(spec_set_ptr, selected_prsm_ptrs,
                                     is_separate, mng_ptr_, test_num_ptr_,
                                     pool_ptr, writer_ptr_vec));
        }
      }
    }
    std::cout << std::flush << "E-value computation - processing " << cnt
              << " of " << spectrum_num << " spectra.\r";
    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
  int remainder = spectrum_num - cnt;
  if (prsm_para_ptr->getGroupSpecNum() > remainder && remainder > 0) {
    // if there are spectrum remaining because they were not combined due to not
    // having enough pairs fix the message as the processing is completed. this
    // code avoids error when no combined spectra is used but a scan is
    // remaining unprocessed because then it will not satisfy the first
    // condition
    std::cout << std::flush << "E-value computation - processing "
              << spectrum_num << " of " << spectrum_num << " spectra.\r";
  }
  pool_ptr->shutDown();
  std::cout << std::endl;
  prsm_reader.close();
  prsm_xml_writer_util::closeWriterPtrVec(writer_ptr_vec);
  writer.close();

  int prsm_top_num = std::numeric_limits<int>::max();
  bool norm = false;
  bool remove_dup = false;
  std::vector<std::string> input_exts;
  for (int t = 0; t < mng_ptr_->thread_num_; t++) {
    input_exts.push_back(mng_ptr_->output_file_ext_ + "_" + std::to_string(t));
  }
  PrsmStrMergePtr merge_ptr = std::make_shared<PrsmStrMerge>(
      sp_file_name, input_exts, mng_ptr_->output_file_ext_, prsm_top_num, norm,
      remove_dup);
  merge_ptr->process();
  merge_ptr = nullptr;

  // remove temporary files
  file_util::cleanTempFiles(sp_file_name, mng_ptr_->output_file_ext_ + "_");
}

bool EValueProcessor::checkPrsms(const PrsmPtrVec& prsm_ptrs) {
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    ExpectedValuePtr extreme_value_ptr = prsm_ptrs[i]->getExpectedValuePtr();
    if (extreme_value_ptr != nullptr) {
      double evalue = extreme_value_ptr->getEValue();
      double frag_num = prsm_ptrs[i]->getMatchTheoFragNum();
      if (evalue <= mng_ptr_->computation_evalue_cutoff &&
          frag_num >= mng_ptr_->computation_frag_num_cutoff) {
        return false;
      }
    }
  }
  return true;
}

// if matched peak number is too small or E-value is 0, replace it
// with a max evalue.
void EValueProcessor::fixupEvalues(PrsmPtrVec& sele_prsm_ptrs) {
  for (unsigned i = 0; i < sele_prsm_ptrs.size(); i++) {
    if (std::round(sele_prsm_ptrs[i]->getMatchTheoFragNum()) <=
        std::round(mng_ptr_->comp_evalue_min_match_theo_frag_num_)) {
      // LOG_DEBUG("Set max e value ");
      sele_prsm_ptrs[i]->setExpectedValuePtr(ExpectedValue::getMaxEvaluePtr());
    } else {
      if (sele_prsm_ptrs[i]->getEValue() == 0.0) {
        LOG_WARN("Invalid e value!");
        sele_prsm_ptrs[i]->setExpectedValuePtr(
            ExpectedValue::getMaxEvaluePtr());
      }
    }
  }
}

void EValueProcessor::compEvalues(const SpectrumSetPtr& spec_set_ptr,
                                  PrsmPtrVec& sele_prsm_ptrs,
                                  bool is_separate) {
  comp_pvalue_ptr_->process(spec_set_ptr, sele_prsm_ptrs, is_separate);
  fixupEvalues(sele_prsm_ptrs);
}

void EValueProcessor::compEvaluesThreadLocal(const SpectrumSetPtr& spec_set_ptr,
                                             PrsmPtrVec& sele_prsm_ptrs,
                                             bool is_separate) {
  // one CompPValueArray per worker thread, exactly as geneTask above: its
  // CompProbValue tables are large and not shareable between threads
  static thread_local CompPValueArrayPtr comp_pvalue_ptr;
  static thread_local const CountTestNum* comp_pvalue_built_for = nullptr;
  if (comp_pvalue_ptr == nullptr ||
      comp_pvalue_built_for != test_num_ptr_.get()) {
    comp_pvalue_ptr = std::make_shared<CompPValueArray>(test_num_ptr_,
                                                        mng_ptr_);
    comp_pvalue_built_for = test_num_ptr_.get();
  }
  comp_pvalue_ptr->process(spec_set_ptr, sele_prsm_ptrs, is_separate);
  fixupEvalues(sele_prsm_ptrs);
}

}  // namespace toppic
