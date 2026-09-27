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

#include "stat/tdgf/cluster_prot_corrector.hpp"

#include <chrono>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "common/base/mod.hpp"
#include "common/thread/simple_thread_pool.hpp"
#include "common/util/file_util.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "para/sp_para.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "seq/fasta_index_reader.hpp"
#include "seq/proteoform.hpp"
#include "stat/tdgf/evalue_processor.hpp"

namespace toppic {

namespace cluster_prot_corrector {

int process(const PrsmParaPtr& prsm_para_ptr, const TdgfMngPtr& tdgf_mng_ptr,
            const std::string& input_file_ext,
            const std::string& output_file_ext) {
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string base_name = file_util::basename(sp_file_name);
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(
      prsm_para_ptr->getSearchDbFileNameWithFolder());
  ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
  PrsmPtrVec prsm_ptrs = prsm_reader_util::readAllPrsms(
      base_name + "." + input_file_ext, seq_reader, fix_mod_list);

  // the best PrSM of every proteoform cluster
  std::map<int, PrsmPtr> best_ptrs;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int cluster_id = prsm_ptrs[i]->getProteoformPtr()->getProteoClusterId();
    auto it = best_ptrs.find(cluster_id);
    if (it == best_ptrs.end() ||
        Prsm::cmpEValueIncProtInc(prsm_ptrs[i], it->second)) {
      best_ptrs[cluster_id] = prsm_ptrs[i];
    }
  }
  // the PrSMs to correct, by spectrum id, with the proteoform they get
  std::map<int, std::vector<std::pair<PrsmPtr, ProteoformPtr>>> corrections;
  int corrected_num = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    ProteoformPtr form_ptr = prsm_ptrs[i]->getProteoformPtr();
    PrsmPtr best_ptr = best_ptrs[form_ptr->getProteoClusterId()];
    if (form_ptr->getProtId() != best_ptr->getProteoformPtr()->getProtId()) {
      corrections[prsm_ptrs[i]->getSpectrumId()].emplace_back(
          prsm_ptrs[i], best_ptr->getProteoformPtr());
      corrected_num++;
    }
  }

  if (corrected_num > 0) {
    EValueProcessorPtr evalue_processor =
        std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    evalue_processor->init();
    SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
    double ppo = sp_para_ptr->getPrecTolerancePtr()->getPpo();
    // the spectra as the search and the E-value computation saw them: masses
    // below the EnvCNN cutoff removed, close masses combined with
    // --frag-mass-dedup
    MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(
        sp_file_name, prsm_para_ptr->getGroupSpecNum(),
        sp_para_ptr->getActivationPtr(), sp_para_ptr->getEnvCnnCutoff(),
        sp_para_ptr->getFragDedupPpo());
    // The E-value dynamic programs ran one after another on this thread
    // while the pool sat idle (8.5 s of the stage's 18 s on the 480
    // benchmark, toppic_test experiment 66); they now run as one task per
    // spectrum, like the E-value stage's own pool. The tasks are
    // independent: a PrSM belongs to one spectrum, so the sets they write
    // E-values into are disjoint, and the proteoform reads are the same
    // sharing the E-value stage's tasks already do. The output file is
    // written after the pool drains, from prsm_ptrs in its original order,
    // so the result is byte-identical.
    SimpleThreadPoolPtr pool_ptr =
        std::make_shared<SimpleThreadPool>(tdgf_mng_ptr->thread_num_);
    DeconvMsPtrVec deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
    while (deconv_ms_ptr_vec.size() != 0) {
      MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
      auto it = corrections.find(header_ptr->getSpecId());
      if (it != corrections.end() && header_ptr->containsPrec()) {
        double prec_mono_mass = header_ptr->getFirstPrecMonoMass() -
                                sp_para_ptr->getNTermLabelMass();
        SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
            deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
        if (spec_set_ptr->isValid()) {
          PrsmPtrVec sele_prsm_ptrs;
          for (size_t i = 0; i < it->second.size(); i++) {
            PrsmPtr prsm_ptr = it->second[i].first;
            ProteoformPtr form_ptr = it->second[i].second;
            prsm_ptr->setDeconvMsPtrVec(spec_set_ptr->getDeconvMsPtrVec());
            // the refined spectra are built for the proteoform's mass, as in
            // the search (the refined spectra do not exist yet, so only the
            // member is set here)
            prsm_ptr->setAdjustedPrecMass(form_ptr->getMass());
            // sets the proteoform, builds the refined spectra and recounts
            // the matched fragment numbers
            prsm_ptr->setProteoformPtr(form_ptr, sp_para_ptr);
            sele_prsm_ptrs.push_back(prsm_ptr);
          }
          while (pool_ptr->getQueueSize() >=
                 static_cast<std::size_t>(tdgf_mng_ptr->thread_num_) + 2) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
          }
          pool_ptr->enqueue([evalue_processor, spec_set_ptr, sele_prsm_ptrs,
                             ppo]() mutable {
            evalue_processor->compEvaluesThreadLocal(spec_set_ptr,
                                                     sele_prsm_ptrs, ppo,
                                                     false);
          });
        }
      }
      deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
    }
    pool_ptr->shutDown();
  }

  PrsmXmlWriter writer(base_name + "." + output_file_ext);
  writer.writeVector(prsm_ptrs);
  writer.close();
  return corrected_num;
}

}  // namespace cluster_prot_corrector

}  // namespace toppic
