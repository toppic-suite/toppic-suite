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

#include "base/acid_base.hpp"
#include "base/mod_util.hpp"
#include "base/file_util.hpp"
#include "base/residue_util.hpp"
#include "base/neutral_loss.hpp"
#include "base/base_algo.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/theo_peak.hpp"
#include "spec/theo_peak_util.hpp"
#include "spec/extend_ms_factory.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_xml_writer.hpp"

#include "mcmc/mcmc_dpr_processor.hpp"

namespace prot {

ProteoformPtr DprProcessor::randomTrans(ProteoformPtr prot_form) {
  return prot_form;
}

void DprProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  sp_para_ptr->prec_error_ = 0;
  ppo_ = sp_para_ptr->getPeakTolerancePtr()->getPpo();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->residue_mod_file_;
  residue_vec_ = residue_util::convertStrToResiduePtrVec("ARNDCEQGHILKMFPSTWYVUO",
                                                         prsm_para_ptr->getFixModPtrVec());
  ModPtrVec var_mod_ptr_vec = mod_util::readModTxt(var_mod_file_name)[2];
  for (size_t j = 0; j < var_mod_ptr_vec.size(); j++) {
    residue_vec_.push_back(var_mod_ptr_vec[j]->getModResiduePtr());
  }

  residue_dist_ = std::uniform_int_distribution<int32_t>(0, residue_vec_.size() - 1);

  LOG_DEBUG("residue_set_ size: " << residue_vec_.size());

  PrsmReaderPtr prsm_reader
      = std::make_shared<PrsmReader>(file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_);

  FastaIndexReaderPtr fasta_reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());

  // no multi-spec support now
  MsAlignReaderPtr sp_reader_ptr = std::make_shared<MsAlignReader>(sp_file_name,
                                                                   1,  // prsm_para_ptr->getGroupSpecNum()
                                                                   sp_para_ptr->getActivationPtr(),
                                                                   sp_para_ptr->getSkipList());

  SpectrumSetPtr spec_set_ptr = sp_reader_ptr->getNextSpectrumSet(sp_para_ptr)[0];

  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      prec_mass_ = spec_set_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr()->getPrecMonoMass();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        pos_dist_ = std::uniform_int_distribution<int32_t>(0, prsm_ptr->getProteoformPtr()->getLen() - 1);
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        TheoPeakPtrVec theo_peak_ptrs
            = theo_peak_util::geneProteoformTheoPeak(prsm_ptr->getProteoformPtr(),
                                                     deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getActivationPtr(),
                                                     sp_para_ptr->getMinMass());

        std::vector<double> theo_masses = theo_peak_util::getTheoMassVec(theo_peak_ptrs);
        std::sort(theo_masses.begin(), theo_masses.end());
        ExtendMsPtrVec refine_ms_ptr_vec = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, 
                                                                              sp_para_ptr, 
                                                                              prsm_ptr->getAdjustedPrecMass());
        std::vector<double> ms_masses = ExtendMs::getExtendMassVec(refine_ms_ptr_vec[0]);
        std::sort(ms_masses.begin(), ms_masses.end());
        int scr = base_algo::compNumMatchedTheoMasses(ms_masses, theo_masses, ppo_);
        //scr = prot::peak_ion_pair_util::compMatchFragNum(pairs);
        LOG_DEBUG("matching score: " << scr);
        // mu
        std::vector<double> mu(theo_masses.size(), 1.0);
        std::vector<double> p(theo_masses.size(), 1.0);
        std::vector<int> n(theo_masses.size(), 0);

        prsm_ptr = prsm_reader->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());
      }
    }
    spec_set_ptr = sp_reader_ptr->getNextSpectrumSet(sp_para_ptr)[0];
  }
  sp_reader_ptr->close();
}

void DprProcessor::simulateDPR(ProteoformPtr p, const DeconvMsPtrVec & deconv_ms_ptr_vec,
                               const std::vector<double> & ms_masses,
                               long omega, int N, const std::vector<double> & mu) {
  TheoPeakPtrVec theo_peak_ptrs
      = theo_peak_util::geneProteoformTheoPeak(p, deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getActivationPtr(),
                                               mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getMinMass());

  std::vector<double> theo_masses = theo_peak_util::getTheoMassVec(theo_peak_ptrs);

  int score = base_algo::compNumMatchedTheoMasses(ms_masses, theo_masses, ppo_);

  ProteoformPtr p2;
  while (this->z_ < N) {
    p2 = randomTrans(p);
    theo_peak_ptrs
        = theo_peak_util::geneProteoformTheoPeak(p2, deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getActivationPtr(),
                                                 mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getMinMass());

    theo_masses = theo_peak_util::getTheoMassVec(theo_peak_ptrs);

    ExtendMsPtrVec refine_ms_ptr_vec = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, 
                                                                          mng_ptr_->prsm_para_ptr_->getSpParaPtr(),
                                                                          p2->getMass());

    std::vector<double> spec_masses = ExtendMs::getExtendMassVec(refine_ms_ptr_vec[0]);

    int score2 = base_algo::compNumMatchedTheoMasses(spec_masses, theo_masses, ppo_);

    if (mu[score2] < omega) return;

    if (mu[score2] > mu[score]) {
      long Y = std::round(mu[score2] / mu[score]); 
      std::uniform_int_distribution<long> omega_dist(static_cast<long>(mu[score]), static_cast<long>(mu[score2]));
      for (long i = 1; i < Y; i++) {
        long omega2 = omega_dist(mt_);
        simulateDPR(p2, deconv_ms_ptr_vec, ms_masses, omega2, N, mu);
      }
    }
    this->z_++;
    score_vec_[z_] = score2;
  }
}

}  // namespace prot
