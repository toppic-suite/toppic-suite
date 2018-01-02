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

#include "base/ptm_util.hpp"
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

double getResidueVecMass(const ResiduePtrVec &residues) {
  double m = 0.0;
  for (size_t i = 0; i < residues.size(); i++) {
    m += residues[i]->getMass();
  }
  return m;
}

ResiduePtrVec DprProcessor::randomTrans(ResiduePtrVec &residues) {
  std::uniform_int_distribution<size_t> dis(0, residues.size() - 1);

  size_t p = dis(mt_);

  ResiduePtr res = residues[p];

  size_t idx = residue_idx_map_[res->getAminoAcidPtr()->getOneLetter()]; 

  double cur_mass = getResidueVecMass(residues);

  if (cur_mass < pep_mass_ - mng_ptr_->mass_limit_) {
    if (idx == residue_vec_.size() - 1) return residues;
    std::uniform_int_distribution<size_t> change_dis(idx + 1, residue_vec_.size() - 1);
    residues[idx] = residue_vec_[change_dis(mt_)];
  } else if (cur_mass > pep_mass_ + mng_ptr_->mass_limit_) {
    if (idx == 1) return residues;
    std::uniform_int_distribution<size_t> change_dis(1, idx - 1);
    residues[idx] = residue_vec_[change_dis(mt_)];
  } else {
    std::uniform_int_distribution<size_t> change_dis(1, residue_vec_.size() - 1);
    residues[idx] = residue_vec_[change_dis(mt_)];
  }
  return residues;
}

void DprProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;

  sp_para_ptr_ = prsm_para_ptr->getSpParaPtr();
  sp_para_ptr_->prec_error_ = 0;

  ppo_ = sp_para_ptr_->getPeakTolerancePtr()->getPpo();

  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string var_mod_file_name = mng_ptr_->residue_mod_file_;

  residue_vec_ = ResidueBase::getBaseNonePtmResiduePtrVec();
  std::sort(residue_vec_.begin(), residue_vec_.end(),
            [](ResiduePtr a, ResiduePtr b) {
            return a->getMass() < b->getMass();
            });

  for (size_t i = 0; i < residue_vec_.size(); i++) {
    LOG_DEBUG(residue_vec_[i]->getAminoAcidPtr()->getOneLetter() << " " << residue_vec_[i]->getMass());
    residue_idx_map_[residue_vec_[i]->getAminoAcidPtr()->getOneLetter()] = i;
  }

  ptm_vec_ = ptm_util::readPtmTxt(var_mod_file_name);
  LOG_DEBUG("ptm_vec_ size: " << ptm_vec_.size());

  ModPtrVec var_mod_ptr_vec = mod_util::readModTxt(var_mod_file_name)[2];

  for (size_t i = 0; i < var_mod_ptr_vec.size(); i++) {
    ptm_residue_map_[var_mod_ptr_vec[i]->getModResiduePtr()->getPtmPtr()].push_back(var_mod_ptr_vec[i]->getOriResiduePtr()); 
  }

  PrsmReaderPtr prsm_reader
      = std::make_shared<PrsmReader>(file_util::basename(sp_file_name) + "." + mng_ptr_->input_file_ext_);

  FastaIndexReaderPtr fasta_reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());

  // no multi-spec support now
  MsAlignReaderPtr sp_reader_ptr = std::make_shared<MsAlignReader>(sp_file_name,
                                                                   1,  // prsm_para_ptr->getGroupSpecNum()
                                                                   sp_para_ptr_->getActivationPtr(),
                                                                   sp_para_ptr_->getSkipList());

  SpectrumSetPtr spec_set_ptr = sp_reader_ptr->getNextSpectrumSet(sp_para_ptr_)[0];

  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      prec_mass_ = spec_set_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr()->getPrecMonoMass();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        ExtendMsPtrVec refine_ms_ptr_vec = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, 
                                                                              sp_para_ptr_, 
                                                                              prec_mass_);
        std::vector<double> ms_masses = ExtendMs::getExtendMassVec(refine_ms_ptr_vec[0]);
        std::sort(ms_masses.begin(), ms_masses.end());

        ActivationPtr act = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getActivationPtr();

        processOnePrsm(prsm_ptr, act, ms_masses);

        prsm_ptr = prsm_reader->readOnePrsm(fasta_reader_ptr, prsm_para_ptr->getFixModPtrVec());
      }
    }
    spec_set_ptr = sp_reader_ptr->getNextSpectrumSet(sp_para_ptr_)[0];
  }
  sp_reader_ptr->close();
}

std::vector<double> getNTheoMassVec(const ResiduePtrVec &residues, IonTypePtr n_ion_type_ptr, double min_mass) {
  std::vector<double> theo_masses;
  ResSeqPtr res_seq = std::make_shared<ResidueSeq>(residues);
  BpSpecPtr bp_spec = std::make_shared<BpSpec>(res_seq);
  BreakPointPtrVec bps = bp_spec->getBreakPointPtrVec();
  double max_mass = res_seq->getSeqMass() - min_mass;
  for (size_t i = 0; i < bps.size(); i++) {
    double n_mass = bps[i]->getNTermMass(n_ion_type_ptr);
    if (min_mass <= n_mass && n_mass <= max_mass) {
      theo_masses.push_back(n_mass); 
    }
  }
  return theo_masses;
}

std::vector<double> getCTheoMassVec(const ResiduePtrVec &residues, IonTypePtr c_ion_type_ptr, double min_mass) {
  std::vector<double> theo_masses;
  ResSeqPtr res_seq = std::make_shared<ResidueSeq>(residues);
  BpSpecPtr bp_spec = std::make_shared<BpSpec>(res_seq);
  BreakPointPtrVec bps = bp_spec->getBreakPointPtrVec();
  double max_mass = res_seq->getSeqMass() - min_mass;
  for (size_t i = 0; i < bps.size(); i++) {
    double n_mass = bps[i]->getCTermMass(c_ion_type_ptr);
    if (min_mass <= n_mass && n_mass <= max_mass) {
      theo_masses.push_back(n_mass); 
    }
  }
  return theo_masses;
}

void DprProcessor::processOnePrsm(PrsmPtr prsm_ptr, ActivationPtr act, const std::vector<double> & ms_masses) {
  ProteoformPtr prot_form = prsm_ptr->getProteoformPtr();
  pep_mass_ = getResidueVecMass(prot_form->getResSeqPtr()->getResidues());
  ChangePtrVec change_vec = prot_form->getChangePtrVec();
  PtmPtrVec ptm_vec(change_vec.size());
  for (size_t i = 0; i < change_vec.size(); i++) {
    ptm_vec[i] = change_vec[i]->getModPtr()->getModResiduePtr()->getPtmPtr();
  }

  std::vector<double> theo_masses = getNTheoMassVec(prot_form->getResSeqPtr()->getResidues(),
                                                    act->getNIonTypePtr(),
                                                    sp_para_ptr_->getMinMass());

  prot::ResiduePtrVec residues = prot_form->getResSeqPtr()->getResidues();

  int scr = getMaxScore(residues, ms_masses, act, ptm_vec);
  LOG_DEBUG("matching score: " << scr);
  std::cin.get();

  score_vec_.resize(mng_ptr_->N_);
  std::fill(score_vec_.begin(), score_vec_.end(), 0);

  std::vector<double> mu(mng_ptr_->n_, 1.0);
  std::vector<double> p(mng_ptr_->n_, 1.0);
  std::vector<int> n(mng_ptr_->n_, 0);

  for (int k = 0; k < mng_ptr_->k_; k++) {
    std::fill(n.begin(), n.end(), 0);
    std::fill(score_vec_.begin(), score_vec_.end(), 0);
    score_vec_[0] = scr;
    residues = prot_form->getResSeqPtr()->getResidues();

    if (k != 0) {
      for (size_t i = 0; i < mu.size(); i++) {
        mu[i] = 1 / (p[i] + 1e-5);
      }

    }
    this->z_ = 0;
    std::cout << "mu min " << *std::min_element(mu.begin(), mu.end()) << std::endl;

    simulateDPR(residues, ms_masses, act, ptm_vec,
                *std::min_element(mu.begin(), mu.end()), mu);

    for (size_t i = 0; i < score_vec_.size(); i++) {
      n[score_vec_[i]]++; 
    }

    for (size_t i = 0; i < n.size(); i++) {
      p[i] = n[i] * 1.0 / mu[i]; 
    }

    double sum = std::accumulate(p.begin(), p.end(), 0.0);
    for (size_t i = 0; i < p.size(); i++) {
      p[i] = p[i] / sum;
    }

    std::cout << "p-value: " << std::accumulate(p.begin() + scr, p.end(), 0.0) << std::endl;
  }
}

int compScore(const std::vector<double> & ms_masses,
              std::vector<double> n_theo_masses, std::vector<double> c_theo_masses,
              const std::vector<size_t> & change_pos, const PtmPtrVec & ptm_vec, double ppo) {
  std::vector<double> change_masses(n_theo_masses.size(), 0.0);
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    double m = ptm_vec[i]->getMonoMass();
    std::for_each(change_masses.begin() + change_pos[i],
                  change_masses.end(), [m](double& d) { d += m;});
  }

  for (size_t i = 0; i < n_theo_masses.size(); i++) {
    n_theo_masses[i] += change_masses[i]; 
  }

  int scr = base_algo::compNumMatchedTheoMasses(ms_masses, n_theo_masses, ppo);

  change_masses.resize(c_theo_masses.size());
  std::fill(change_masses.begin(), change_masses.end(), 0.0);
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    double m = ptm_vec[i]->getMonoMass();
    std::for_each(change_masses.begin() + change_masses.size() - change_pos[i],
                  change_masses.end(), [m](double& d) { d += m;});
  }

  for (size_t i = 0; i < c_theo_masses.size(); i++) {
    c_theo_masses[i] += change_masses[i];
  }

  scr += base_algo::compNumMatchedTheoMasses(ms_masses, c_theo_masses, ppo);

  return scr;
}

int DprProcessor::getMaxScore(const ResiduePtrVec &residues, const std::vector<double> & ms_masses,
                              ActivationPtr act, const PtmPtrVec & ptm_vec) {
  std::vector<double> n_theo_masses = getNTheoMassVec(residues, act->getNIonTypePtr(),
                                                      sp_para_ptr_->getMinMass());
  std::sort(n_theo_masses.begin(), n_theo_masses.end());

  std::vector<double> c_theo_masses = getCTheoMassVec(residues, act->getCIonTypePtr(),
                                                      sp_para_ptr_->getMinMass());
  std::sort(c_theo_masses.begin(), c_theo_masses.end());

  std::vector<std::vector<size_t> > possible_change_pos(ptm_vec.size());
  std::vector<size_t> change_pos(ptm_vec.size());
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    ResiduePtrVec possible_res = ptm_residue_map_[ptm_vec[i]];
    for (size_t k = 0; k < residues.size(); k++) {
      if (std::find(possible_res.begin(), possible_res.end(), residues[k]) != possible_res.end()) {
        possible_change_pos[i].push_back(k);
      } 
    }
  }

  // random init positions
  for (size_t i = 0; i < possible_change_pos.size(); i++) {
    std::uniform_int_distribution<size_t> dis(0, possible_change_pos[i].size() - 1);
    change_pos[i] = dis(mt_);
  }

  int max_scr = compScore(ms_masses, n_theo_masses, c_theo_masses,
                          change_pos, ptm_vec, ppo_);

  // search the max score using greedy method
  for (size_t i = 0; i < possible_change_pos.size(); i++) {
    size_t max_idx = change_pos[i];
    for (size_t k = 0; k < possible_change_pos[i].size(); k++) {
      change_pos[i] = possible_change_pos[i][k];
      int scr = compScore(ms_masses, n_theo_masses, c_theo_masses,
                          change_pos, ptm_vec, ppo_);
      if (scr > max_scr) {
        max_idx = possible_change_pos[i][k];
        max_scr = scr;
      }
    } 
    change_pos[i] = max_idx;
  }

  return max_scr;
}

void DprProcessor::simulateDPR(ResiduePtrVec &residues, const std::vector<double> & ms_masses,
                               ActivationPtr act, const PtmPtrVec & ptm_vec,
                               long omega, const std::vector<double> & mu) {
  int score = getMaxScore(residues, ms_masses, act, ptm_vec);
  while (this->z_ < mng_ptr_->N_) {
    while (true) {
      residues = randomTrans(residues);
      if (std::abs(getResidueVecMass(residues) - pep_mass_) <= prec_mass_ * ppo_) {
        break;
      }
    }

    int score2 = getMaxScore(residues, ms_masses, act, ptm_vec);

    if (mu[score2] < omega) return;

    if (mu[score2] > mu[score]) {
      long Y = std::round(mu[score2] / mu[score]); 
      std::uniform_int_distribution<long> omega_dist(static_cast<long>(mu[score]), static_cast<long>(mu[score2]));
      for (long i = 1; i < Y; i++) {
        long omega2 = omega_dist(mt_);
        simulateDPR(residues, ms_masses, act, ptm_vec, omega2, mu);
      } 
    }
    this->z_++;
    this->score_vec_[this->z_] = score2;
  }
}

}  // namespace prot
