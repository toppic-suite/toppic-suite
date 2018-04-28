//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include <random>
#include <string>
#include <algorithm>
#include <vector>
#include <map>

#include <boost/algorithm/string.hpp>

#include "base/ptm_util.hpp"
#include "base/mod_util.hpp"
#include "base/file_util.hpp"
#include "base/residue_util.hpp"
#include "base/neutral_loss.hpp"
#include "base/base_algo.hpp"
#include "base/extreme_value.hpp"
#include "base/thread_pool.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/theo_peak.hpp"
#include "spec/theo_peak_util.hpp"
#include "spec/extend_ms_factory.hpp"

#include "prsm/prsm_reader.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_xml_writer.hpp"

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/count_test_num.hpp"

#include "mcmc/mcmc_dpr_processor.hpp"
#include "mcmc/comp_pvalue_mcmc.hpp"

namespace prot {

void DprProcessor::init() {
  std::string var_mod_file_name = mng_ptr_->residue_mod_file_;

  ptm_vec_ = ptm_util::readPtmTxt(var_mod_file_name);

  ModPtrVec var_mod_ptr_vec = mod_util::readModTxt(var_mod_file_name)[2];

  for (size_t i = 0; i < var_mod_ptr_vec.size(); i++) {
    ptm_residue_map_[var_mod_ptr_vec[i]->getModResiduePtr()->getPtmPtr()].push_back(var_mod_ptr_vec[i]->getOriResiduePtr());
  }

  read_mass_table();

  TdgfMngPtr tdgf_mng_ptr
      = std::make_shared<TdgfMng>(mng_ptr_->prsm_para_ptr_, 0, 0.0, false, false, 1, "", "");

  test_num_ptr_ = std::make_shared<CountTestNum>(tdgf_mng_ptr);

  ptm_mass_vec2d_ = compPtmComb();
}

std::vector<std::vector<double> > DprProcessor::compPtmComb() {
  std::vector<std::vector<int> > mass_ptm_vec2d_int(mng_ptr_->max_known_mods_ + 1);
  std::vector<std::vector<double> > ptm_mass_vec2d(mng_ptr_->max_known_mods_ + 1);
  int inte_tole = mng_ptr_->getIntTolerance();

  for (size_t i = 0; i < ptm_vec_.size(); i++) {
    int mass = std::ceil(ptm_vec_[i]->getMonoMass() * mng_ptr_->convert_ratio_);
    bool found = false;
    for (auto it = mass_ptm_vec2d_int[1].begin(); it != mass_ptm_vec2d_int[1].end(); it++) {
      if (std::abs(*it - mass) <= inte_tole) {
        found = true;
        break;
      }
    }

    if (!found) {
      mass_ptm_vec2d_int[1].push_back(mass);
    }
  }

  for (int k = 2; k <= mng_ptr_->max_known_mods_; k++) {
    std::vector<int> cur_vec_mass;
    for (auto it = mass_ptm_vec2d_int[k - 1].begin(); it != mass_ptm_vec2d_int[k - 1].end(); it++) {
      cur_vec_mass.push_back(*it);
    }

    for (size_t i = 0; i < ptm_vec_.size(); i++) {
      int mass = std::ceil(ptm_vec_[i]->getMonoMass() * mng_ptr_->convert_ratio_);
      for (size_t cnt = 0; cnt < cur_vec_mass.size(); cnt++) {
        int new_mass = mass + cur_vec_mass[cnt];
        bool found = false;
        auto it2 = mass_ptm_vec2d_int[k].begin();
        while (it2 != mass_ptm_vec2d_int[k].end() && !found) {
          if (std::abs(*it2 - new_mass) <= inte_tole) {
            found = true;
          }
          it2++;
        }
        if (!found) {
          mass_ptm_vec2d_int[k].push_back(new_mass);
        }
      }
    }
    std::sort(mass_ptm_vec2d_int[k].begin(), mass_ptm_vec2d_int[k].end());
  }

  for (size_t i = 0; i < mass_ptm_vec2d_int.size(); i++) {
    for (size_t k = 0; k < mass_ptm_vec2d_int[i].size(); k++) {
      ptm_mass_vec2d[i].push_back(mass_ptm_vec2d_int[i][k] / mng_ptr_->convert_ratio_);
    }
  }

  return ptm_mass_vec2d;
}

void DprProcessor::read_mass_table() {
  std::string mass_table_file
      = mng_ptr_->prsm_para_ptr_->getResourceDir() + "/base_data/mass_table.txt";

  std::ifstream infile(mass_table_file);
  std::string line;
  while (std::getline(infile, line)) {
    std::vector<std::string> strs;
    boost::split(strs, line, boost::is_any_of("\t"));
    int m = std::stoi(strs[0]);
    for (size_t i = 1; i < strs.size(); i++) {
      mass_table_[m].push_back(strs[i]);
    }
  }
  infile.close();
}

void DprProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;

  sp_para_ptr_ = prsm_para_ptr->getSpParaPtr();

  sp_para_ptr_->prec_error_ = 0;

  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_;
  std::string output_file_name = file_util::basename(sp_file_name) + "." + mng_ptr_->output_file_ext_;

  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);

  PrsmXmlWriterPtr prsm_writer
      = std::make_shared<PrsmXmlWriter>(output_file_name + "_" + std::to_string(mng_ptr_->thread_num_));

  pool_ptr_ = std::make_shared<ThreadPool<PrsmXmlWriter> >(mng_ptr_->thread_num_ , output_file_name);

  FastaIndexReaderPtr fasta_reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());

  // no multi-spec support now
  MsAlignReaderPtr sp_reader_ptr = std::make_shared<MsAlignReader>(sp_file_name,
                                                                   1,  // prsm_para_ptr->getGroupSpecNum()
                                                                   sp_para_ptr_->getActivationPtr(),
                                                                   sp_para_ptr_->getSkipList());

  int spectrum_num = msalign_util::getSpNum(sp_file_name);

  int cnt = 0;

  SpectrumSetPtr spec_set_ptr = sp_reader_ptr->getNextSpectrumSet(sp_para_ptr_)[0];

  while (spec_set_ptr != nullptr) {
    cnt += prsm_para_ptr->getGroupSpecNum();
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();

        ExtendMsPtrVec refine_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,
                                                   sp_para_ptr_,
                                                   prsm_ptr->getAdjustedPrecMass());

        processOnePrsm(prsm_ptr, spec_set_ptr, prsm_writer);

        prsm_ptr = prsm_reader->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());
      }
    }
    spec_set_ptr = sp_reader_ptr->getNextSpectrumSet(sp_para_ptr_)[0];

    std::cout << std::flush << "E-value computation - processing " << cnt << " of "
        << spectrum_num << " spectra.\r";
  }
  pool_ptr_->ShutDown();
  std::cout << std::endl;
  sp_reader_ptr->close();
  prsm_reader->close();
  prsm_writer->close();

  PrsmXmlWriterPtr all_writer_ptr = std::make_shared<PrsmXmlWriter>(output_file_name);
  PrsmPtrVec prsm_vec;

  for (int i = 0; i <= mng_ptr_->thread_num_; i++) {
    PrsmReaderPtr all_reader_ptr = std::make_shared<PrsmReader>(output_file_name + "_" + std::to_string(i));
    PrsmPtr p = all_reader_ptr->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());

    while (p != nullptr) {
      prsm_vec.push_back(p);
      p = all_reader_ptr->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());
    }

    all_reader_ptr->close();
  }

  std::sort(prsm_vec.begin(), prsm_vec.end(), Prsm::cmpSpectrumIdIncEvalueInc);

  all_writer_ptr->writeVector(prsm_vec);

  all_writer_ptr->close();
}

std::function<void()> geneTask(SpectrumSetPtr spec_set_ptr,
                               PrsmPtr prsm_ptr,
                               MCMCMngPtr mng_ptr,
                               const std::map<PtmPtr, std::vector<ResiduePtr> > & ptm_residue_map,
                               const std::map<int, std::vector<std::string> > & mass_table,
                               CountTestNumPtr test_num_ptr,
                               const std::vector<std::vector<double> > & ptm_mass_vec2d,
                               std::shared_ptr<ThreadPool<PrsmXmlWriter> > pool_ptr) {
  return [spec_set_ptr, prsm_ptr, mng_ptr, ptm_residue_map, mass_table, test_num_ptr, ptm_mass_vec2d, pool_ptr]() {
    CompPValueMCMCPtr comp_mcmc_ptr
        = std::make_shared<prot::CompPValueMCMC>(mng_ptr, ptm_residue_map, mass_table);
    DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
    ExtendMsPtrVec refine_ms_ptr_vec
        = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,
                                               mng_ptr->prsm_para_ptr_->getSpParaPtr(),
                                               prsm_ptr->getAdjustedPrecMass());

    double ppo = mng_ptr->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr()->getPpo();

    double tolerance = refine_ms_ptr_vec[0]->getMsHeaderPtr()->getErrorTolerance(ppo);

    std::vector<double> ms_masses = extend_ms::getExtendMassVec(refine_ms_ptr_vec[0]);

    std::vector<int> ms_mass_int(ms_masses.size());

    ActivationPtr act = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getActivationPtr();

    for (size_t k = 0; k < ms_masses.size(); k++) {
      ms_mass_int[k] = static_cast<int>(ms_masses[k] * mng_ptr->convert_ratio_) >> 5;
    }

    std::sort(ms_mass_int.begin(), ms_mass_int.end());

    ms_mass_int.erase(std::unique(ms_mass_int.begin(), ms_mass_int.end()), ms_mass_int.end());

    double one_prob = comp_mcmc_ptr->compOneProbMCMC(prsm_ptr, act, ms_mass_int);

    double cand_num;

    AlignTypePtr type_ptr = prsm_ptr->getProteoformPtr()->getAlignType();

    cand_num = test_num_ptr->compCandNum(type_ptr, 0, prsm_ptr->getAdjustedPrecMass() - mass_constant::getWaterMass(),
                                         tolerance);

    std::vector<double> mass_ptm_vec = ptm_mass_vec2d[prsm_ptr->getProteoformPtr()->getVariablePtmNum()];

    for (size_t k = 0; k < mass_ptm_vec.size(); k++) {
      cand_num += test_num_ptr->compCandNum(type_ptr, 0,
                                            prsm_ptr->getAdjustedPrecMass() - mass_constant::getWaterMass() - mass_ptm_vec[k],
                                            tolerance);
    }

    if (cand_num == 0) {cand_num = 1;}

    LOG_DEBUG("cand_num " << cand_num);
    ExtremeValuePtr evalue = std::make_shared<ExtremeValue>(one_prob, cand_num, 1);
    prsm_ptr->setExtremeValuePtr(evalue);

    boost::thread::id thread_id = boost::this_thread::get_id();
    PrsmXmlWriterPtr writer_ptr = pool_ptr->getWriter(thread_id);

    writer_ptr->write(prsm_ptr);
  };
}

void DprProcessor::processOnePrsm(PrsmPtr prsm_ptr, SpectrumSetPtr spec_set_ptr,
                                  PrsmXmlWriterPtr prsm_writer) {
  if (prsm_ptr->getMatchFragNum() < 4) {
    prsm_ptr->setExtremeValuePtr(ExtremeValue::getMaxEvaluePtr());
    prsm_writer->write(prsm_ptr);
    return;
  }

  int peak_num = 0;

  DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();

  for (size_t k = 0; k < deconv_ms_ptr_vec.size(); k++) {
    peak_num += deconv_ms_ptr_vec[k]->size();
  }

  if (peak_num > 500) {
    if (prsm_ptr->getProteoformPtr()->getVariablePtmNum() > 0 && prsm_ptr->getMatchFragNum() < peak_num / 10.0) {
      prsm_ptr->setExtremeValuePtr(ExtremeValue::getMaxEvaluePtr());
      prsm_writer->write(prsm_ptr);
      return;
    }
  }

  while (pool_ptr_->getQueueSize() >= mng_ptr_->thread_num_ + 2) {
    boost::this_thread::sleep(boost::posix_time::milliseconds(100));
  }
  pool_ptr_->Enqueue(geneTask(spec_set_ptr, prsm_ptr, mng_ptr_,
                              ptm_residue_map_, mass_table_,
                              test_num_ptr_, ptm_mass_vec2d_, pool_ptr_));
}

}  // namespace prot
