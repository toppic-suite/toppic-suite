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

#include "filter/tag/tag_filter_processor.hpp"

#include <iostream>
#include <string>
#include <vector>

#include "common/base/proteoform_type.hpp"
#include "common/thread/ordered_task_runner.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "filter/tag/tag_filter.hpp"
#include "filter/tag/tag_index.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace toppic {

namespace tag_filter_processor {

namespace {

struct SpecResult {
  std::vector<SimplePrsmPtrVec> prsms;
  TagFilter::Stats stats;
  int set_num = 0;
};

}  // namespace

void process(const TagFilterMngPtr& mng_ptr) {
  PrsmParaPtr prsm_para_ptr = mng_ptr->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string stage_name = TagFilterMng::getStageName(mng_ptr->stage_);
  std::string prefix = "Tag filtering (" + stage_name + ") - ";

  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  TagIndexPtr index_ptr = std::make_shared<TagIndex>(
      db_file_name, prsm_para_ptr->getFixModPtrVec());
  std::cout << prefix << "database indexed: " << index_ptr->getProtNum()
            << " proteins, " << index_ptr->getResidueNum() << " residues."
            << std::endl;
  TagFilter filter(index_ptr, mng_ptr);

  // the output files
  std::string output_base =
      file_util::basename(sp_file_name) + "." + mng_ptr->output_file_ext_;
  std::vector<SimplePrsmXmlWriterPtr> writer_ptrs;
  if (mng_ptr->writeTypeFiles()) {
    const std::vector<ProteoformTypePtr> type_ptrs = {
        ProteoformType::COMPLETE, ProteoformType::PREFIX,
        ProteoformType::SUFFIX, ProteoformType::INTERNAL};
    for (size_t t = 0; t < type_ptrs.size(); t++) {
      writer_ptrs.push_back(std::make_shared<SimplePrsmXmlWriter>(
          output_base + "_" + type_ptrs[t]->getName()));
    }
  } else {
    writer_ptrs.push_back(std::make_shared<SimplePrsmXmlWriter>(output_base));
  }

  int spectrum_num = msalign_util::getSpNum(sp_file_name);
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr(),
      sp_para_ptr->getEnvCnnCutoff(), sp_para_ptr->getFragDedupPpo());
  const std::vector<double>& prec_error_vec = mng_ptr->getPrecErrorVec();

  // the spectra are scored in parallel and their candidates written on
  // this thread in spectrum order
  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(mng_ptr->thread_num_);
  OrderedTaskRunner<SpecResult> runner(pool_ptr, mng_ptr->thread_num_ * 4);
  int cnt = 0;
  int scored_spec_num = 0;
  int skipped_spec_num = 0;
  int cand_spec_num = 0;
  size_t cand_num = 0;
  size_t tag_num = 0;
  size_t occurrence_num = 0;
  size_t set_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() > 0) {
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      int spec_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getSpecId();
      if (!mng_ptr->searchSpectrum(spec_id)) {
        skipped_spec_num++;
      } else {
        SpectrumSetPtrVec spec_set_vec =
            spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError(
                deconv_ms_ptr_vec, sp_para_ptr, prec_error_vec);
        SpectrumSetPtrVec valid_sets;
        for (size_t k = 0; k < spec_set_vec.size(); k++) {
          if (spec_set_vec[k]->isValid()) {
            valid_sets.push_back(spec_set_vec[k]);
          }
        }
        if (!valid_sets.empty()) {
          scored_spec_num++;
          runner.submit(
              [&filter, valid_sets]() {
                SpecResult result;
                result.prsms.resize(4);
                for (size_t k = 0; k < valid_sets.size(); k++) {
                  TagFilter::Stats stats;
                  std::vector<SimplePrsmPtrVec> prsms =
                      filter.process(valid_sets[k], &stats);
                  for (size_t t = 0; t < prsms.size(); t++) {
                    result.prsms[t].insert(result.prsms[t].end(),
                                           prsms[t].begin(), prsms[t].end());
                  }
                  result.stats.tag_num += stats.tag_num;
                  result.stats.occurrence_num += stats.occurrence_num;
                  result.set_num++;
                }
                return result;
              },
              [&](SpecResult& result) {
                size_t spec_cand_num = 0;
                for (size_t t = 0; t < writer_ptrs.size(); t++) {
                  writer_ptrs[t]->write(result.prsms[t]);
                  spec_cand_num += result.prsms[t].size();
                }
                if (spec_cand_num > 0) {
                  cand_spec_num++;
                }
                cand_num += spec_cand_num;
                tag_num += result.stats.tag_num;
                occurrence_num += result.stats.occurrence_num;
                set_num += result.set_num;
              });
        }
      }
    }
    cnt += group_spec_num;
    std::cout << std::flush << prefix << "processing " << cnt << " of "
              << spectrum_num << " spectra.\r";
    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
  runner.finish();
  pool_ptr->shutDown();
  std::cout << std::endl;
  for (size_t t = 0; t < writer_ptrs.size(); t++) {
    writer_ptrs[t]->close();
  }
  std::cout << prefix << scored_spec_num << " spectra scored";
  if (skipped_spec_num > 0) {
    std::cout << ", " << skipped_spec_num << " skipped";
  }
  std::cout << " (" << set_num << " spectrum sets): " << tag_num
            << " tags, " << occurrence_num << " occurrences, " << cand_num
            << " candidates for " << cand_spec_num << " spectra." << std::endl;
}

}  // namespace tag_filter_processor

}  // namespace toppic
