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

#include <string>
#include <algorithm>
#include <vector>
#include <functional>

#include "base/ptm.hpp"
#include "base/ptm_util.hpp"
#include "base/mod_util.hpp"
#include "base/residue_util.hpp"
#include "base/local_anno.hpp"
#include "base/mass_shift.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_base.hpp"
#include "base/proteoform_factory.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_util.hpp"

#include "prsm/prsm.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/peak_ion_pair_util.hpp"

#include "local_processor.hpp"
#include "local_util.hpp"

namespace prot {

void LocalProcessor::init() {
  ptm_vec_ = ptm_util::readPtmTxt(mng_ptr_->residueModFileName_);

  for (size_t i = 0; i < ptm_vec_.size(); i++) {
    for (size_t j = 0; j < ptm_vec_.size(); j++) {
      ptm_pair_vec_.push_back(std::make_pair(ptm_vec_[i], ptm_vec_[j]));
    }
  }

  std::vector<ModPtrVec> mod_ptr_vec2d = mod_util::readModTxt(mng_ptr_->residueModFileName_);
  mod_list_N_ = mod_ptr_vec2d[0];
  mod_list_C_ = mod_ptr_vec2d[1];
  mod_list_any_ = mod_ptr_vec2d[2];
}

void LocalProcessor::process() {
  std::string spec_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();

  std::string input_file_name = file_util::basename(spec_file_name) + "." + mng_ptr_->input_file_ext_;

  std::string output_file_name = file_util::basename(spec_file_name) + "." + mng_ptr_->output_file_ext_;

  PrsmXmlWriterPtr prsm_writer = std::make_shared<PrsmXmlWriter>(output_file_name);

  std::string db_file_name = mng_ptr_->prsm_para_ptr_->getSearchDbFileName();

  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);

  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);

  ModPtrVec fix_mod_list = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_list);

  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();

  MsAlignReader sp_reader(spec_file_name, group_spec_num,
                          mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getActivationPtr(),
                          mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getSkipList());

  SpectrumSetPtr spec_set_ptr;

  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();

  int spectrum_num = msalign_util::getSpNum(mng_ptr_->prsm_para_ptr_->getSpectrumFileName());

  int cnt = 0;

  while ((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0])!= nullptr) {
    cnt += group_spec_num;
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);

        if (prsm_ptr->getProteoformPtr()->getMassShiftNum(MassShiftType::UNEXPECTED) > 0) {
          prsm_ptr = processOnePrsm(prsm_ptr);
        }

        prsm_writer->write(prsm_ptr);
        prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_list);
      }
    }
    std::cout << std::flush << "PTM characterization is processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
  }

  sp_reader.close();
  prsm_reader->close();
  prsm_writer->close();
  std::cout << std::endl;
}

PrsmPtr LocalProcessor::processOnePrsm(PrsmPtr prsm) {
  int mass_shift_num = prsm->getProteoformPtr()->getMassShiftNum(MassShiftType::UNEXPECTED);
  if (mass_shift_num == 1) {
    double mass = prsm->getProteoformPtr()->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0]->getMassShift();

    // if the mass shift is between [-1, 1], we don't characterize it
    if (std::abs(mass) <= 1 + prsm->getAdjustedPrecMass() * ppo_) return prsm;

    return processOnePtm(prsm);
  } else if (mass_shift_num == 2) {
    return processTwoPtm(prsm);
  }
  return prsm;
}

PrsmPtr LocalProcessor::processOnePtm(PrsmPtr prsm) {
  int ori_num_match_ion = local_util::compMatchFragNum(prsm->getProteoformPtr(),
                                                       prsm->getRefineMsPtrVec(),
                                                       mng_ptr_->min_mass_);

  // we will get a nullptr if the mass shift can't be explained by a known variable ptm
  ProteoformPtr one_known_proteoform = processOneKnownPtm(prsm);

  if (one_known_proteoform != nullptr) {
    int new_num_match_ion = local_util::compMatchFragNum(one_known_proteoform,
                                                         prsm->getRefineMsPtrVec(),
                                                         mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - DESC_MATCH_LIMIT) {
      prsm->setProteoformPtr(one_known_proteoform, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  ProteoformPtr two_known_prsm = processTwoKnownPtm(prsm);

  if (two_known_prsm != nullptr) {
    double new_num_match_ion = local_util::compMatchFragNum(two_known_prsm,
                                                            prsm->getRefineMsPtrVec(),
                                                            mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - DESC_MATCH_LIMIT) {
      prsm->setProteoformPtr(two_known_prsm, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  return prsm;
}

// we will get a nullptr if the mass shift can't be explained by a variable ptm
ProteoformPtr LocalProcessor::processOneKnownPtm(PrsmPtr prsm) {
  ProteoformPtr ori_prot_form = prsm->getProteoformPtr();

  MassShiftPtrVec ori_mass_shift_vec = ori_prot_form->getMassShiftPtrVec();

  MassShiftPtrVec unexpected_shift_vec
      = ori_prot_form->getMassShiftPtrVec(MassShiftType::UNEXPECTED);

  MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec,
                                                                   MassShiftType::UNEXPECTED);

  double shift_mass = local_util::compMassShift(unexpected_shift_vec);

  MassShiftPtr unexpected_shift = local_util::geneMassShift(unexpected_shift_vec[0], shift_mass,
                                                            MassShiftType::UNEXPECTED);

  expected_shift_vec.push_back(unexpected_shift);

  std::sort(expected_shift_vec.begin(), expected_shift_vec.end(), MassShift::cmpPosInc);

  ProteoformPtr one_shift_proteoform
      = std::make_shared<Proteoform>(ori_prot_form->getFastaSeqPtr(),
                                     ori_prot_form->getProtModPtr(),
                                     ori_prot_form->getStartPos(),
                                     ori_prot_form->getEndPos(),
                                     ori_prot_form->getResSeqPtr(),
                                     expected_shift_vec);

  double err = prsm->getAdjustedPrecMass() * ppo_;

  PtmPtrVec ptm_vec = local_util::getPtmPtrVecByMass(shift_mass, err, ptm_vec_);

  // try to adjust the N/C terminus
  if (ptm_vec.size() == 0) {
    if (one_shift_proteoform->getProtModPtr()->getType() == ProtModBase::getType_NME_ACETYLATION()
        || one_shift_proteoform->getProtModPtr()->getType() == ProtModBase::getType_M_ACETYLATION()) {
      expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec,
                                                       MassShiftType::UNEXPECTED);
      for (size_t k = 0; k < expected_shift_vec.size(); k++) {
        if (expected_shift_vec[k]->getTypePtr() == MassShiftType::PROTEIN_VARIABLE) {
          shift_mass += expected_shift_vec[k]->getMassShift();
          expected_shift_vec.erase(expected_shift_vec.begin() + k);
          break;
        }
      }

      unexpected_shift = local_util::geneMassShift(unexpected_shift_vec[0], shift_mass,
                                                   MassShiftType::UNEXPECTED);

      expected_shift_vec.push_back(unexpected_shift);

      std::sort(expected_shift_vec.begin(), expected_shift_vec.end(), MassShift::cmpPosInc);

      one_shift_proteoform
          = std::make_shared<Proteoform>(ori_prot_form->getFastaSeqPtr(),
                                         prot::ProtModBase::getProtModPtr_NONE(),
                                         ori_prot_form->getStartPos(),
                                         ori_prot_form->getEndPos(),
                                         ori_prot_form->getResSeqPtr(),
                                         expected_shift_vec);
    }

    one_shift_proteoform = onePtmTermAdjust(one_shift_proteoform, prsm->getRefineMsPtrVec(), shift_mass, err);

    if (std::abs(shift_mass) < mass_constant::getIsotopeMass() + err) {
      one_shift_proteoform->setProteoClusterId(ori_prot_form->getProteoClusterId());
      one_shift_proteoform->setProtId(ori_prot_form->getProtId());
      return one_shift_proteoform;
    }

    ptm_vec = local_util::getPtmPtrVecByMass(shift_mass, err, ptm_vec_);
  }

  // the mass shift can't be explained after adjusting
  if (ptm_vec.size() == 0) return nullptr;

  double raw_scr;
  std::vector<double> scr_vec;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  compOnePtmScr(one_shift_proteoform, extend_ms_ptr_vec, scr_vec, raw_scr, ptm_vec);
  // mass shift can be explained by a variable ptm, but no modifiable site
  if (scr_vec.size() == 0) return nullptr;

  int bgn, end;
  double conf;
  local_util::scrFilter(scr_vec, bgn, end, conf, threshold_);

  if (bgn == -1) return nullptr;

  // it is known ptm, raw_scr * theta_; otherwise raw_scr * (1 - theta_)
  LocalAnnoPtr anno = std::make_shared<LocalAnno>(bgn, end, conf, scr_vec,
                                                  raw_scr * theta_, ptm_vec[0]);

  ChangePtr change = std::make_shared<Change>(anno->getLeftBpPos(),
                                              anno->getRightBpPos() + 1,
                                              MassShiftType::UNEXPECTED, shift_mass,
                                              std::make_shared<Mod>(ResidueBase::getEmptyResiduePtr(),
                                                                    ResidueBase::getEmptyResiduePtr()));

  change->setLocalAnno(anno);

  MassShiftPtr mass_shift = std::make_shared<MassShift>(change->getLeftBpPos(),
                                                        change->getRightBpPos(),
                                                        MassShiftType::UNEXPECTED);

  mass_shift->setChangePtr(change);

  MassShiftPtrVec mass_shift_vec = one_shift_proteoform->getMassShiftPtrVec();

  for (size_t k = 0; k < mass_shift_vec.size(); k++) {
    if (mass_shift_vec[k]->getTypePtr() == MassShiftType::UNEXPECTED) {
      mass_shift_vec[k] = mass_shift;
      break;
    }
  }

  one_shift_proteoform = std::make_shared<Proteoform>(one_shift_proteoform->getFastaSeqPtr(),
                                                      one_shift_proteoform->getProtModPtr(),
                                                      one_shift_proteoform->getStartPos(),
                                                      one_shift_proteoform->getEndPos(),
                                                      one_shift_proteoform->getResSeqPtr(),
                                                      mass_shift_vec);

  one_shift_proteoform->setProteoClusterId(ori_prot_form->getProteoClusterId());
  one_shift_proteoform->setProtId(ori_prot_form->getProtId());
  return one_shift_proteoform;
}

ProteoformPtr LocalProcessor::onePtmTermAdjust(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                               double & shift_mass, double err) {
  int n_trunc_min, n_trunc_max, c_trunc_min, c_trunc_max;
  getNtermTruncRange(proteoform, extend_ms_ptr_vec, n_trunc_min, n_trunc_max);
  getCtermTruncRange(proteoform, extend_ms_ptr_vec, c_trunc_min, c_trunc_max);

  MassShiftPtr mass_shift_ptr = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0];
  MassShiftPtrVec fix_mass_shift_vec
      = local_util::massShiftFilter(proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED),
                                    MassShiftType::UNEXPECTED);

  double ori_mass = mass_shift_ptr->getMassShift();
  shift_mass = ori_mass;
  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  std::string n_seq, c_seq;
  std::vector<bool> ptm_known_vec;
  std::vector<int> c_vec, n_vec;
  std::vector<double> raw_scr_vec;

  for (int i = n_trunc_min; i <= n_trunc_max; i++) {
    for (int j = c_trunc_min; j <= c_trunc_max; j++) {
      if (i < 0) {
        n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + i, -i);
        shift_mass = ori_mass - residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      } else {
        n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, i);
        shift_mass = ori_mass + residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }

      if (j >= 0) {
        c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, j);
        shift_mass = shift_mass - residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      } else {
        c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + j, -j);
        shift_mass = shift_mass + residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }

      if ((shift_mass >= mng_ptr_->max_ptm_mass_ || shift_mass <= mng_ptr_->min_ptm_mass_) && (i != 0 || j != 0))
        continue;

      if (std::abs(shift_mass) < mass_constant::getIsotopeMass() + err) {
        mass_shift_ptr = local_util::geneMassShift(mass_shift_ptr, shift_mass,
                                                   MassShiftType::UNEXPECTED);
        fix_mass_shift_vec.push_back(mass_shift_ptr);
        std::sort(fix_mass_shift_vec.begin(), fix_mass_shift_vec.end(), MassShift::cmpPosInc);
        return proteoform_factory::geneProteoform(proteoform,
                                                  ori_start + i,
                                                  ori_end + j,
                                                  fix_mass_shift_vec,
                                                  mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }

      n_vec.push_back(i);
      c_vec.push_back(j);
      PtmPtrVec ptm_vec_tmp = local_util::getPtmPtrVecByMass(shift_mass, err, ptm_vec_);
      ptm_known_vec.push_back(ptm_vec_tmp.size() > 0);

      mass_shift_ptr = local_util::geneMassShift(mass_shift_ptr, shift_mass,
                                                 MassShiftType::UNEXPECTED);

      MassShiftPtrVec new_mass_shift_vec;
      new_mass_shift_vec.push_back(mass_shift_ptr);
      new_mass_shift_vec.insert(new_mass_shift_vec.end(),
                                fix_mass_shift_vec.begin(), fix_mass_shift_vec.end());
      std::sort(new_mass_shift_vec.begin(), new_mass_shift_vec.end(), MassShift::cmpPosInc);

      proteoform = proteoform_factory::geneProteoform(proteoform,
                                                      ori_start + i,
                                                      ori_end + j,
                                                      new_mass_shift_vec,
                                                      mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

      double raw_scr;
      std::vector<double> scr_vec;
      compOnePtmScr(proteoform, extend_ms_ptr_vec, scr_vec, raw_scr, ptm_vec_tmp);
      raw_scr_vec.push_back(raw_scr);
    }
  }

  for (size_t i = 0; i < raw_scr_vec.size(); i++) {
    if (ptm_known_vec[i])
      raw_scr_vec[i] = raw_scr_vec[i] * mng_ptr_->theta_;
    else
      raw_scr_vec[i] = raw_scr_vec[i] * (1 - mng_ptr_->theta_);
  }

  int idx = std::distance(raw_scr_vec.begin(), std::max_element(raw_scr_vec.begin(), raw_scr_vec.end()));

  if (n_vec[idx] < 0) {
    n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + n_vec[idx], -n_vec[idx]);
    shift_mass = ori_mass - residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, n_vec[idx]);
    shift_mass = ori_mass + residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  }

  if (c_vec[idx] >= 0) {
    c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, c_vec[idx]);
    shift_mass = shift_mass - residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + c_vec[idx], -c_vec[idx]);
    shift_mass = shift_mass + residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  }

  mass_shift_ptr = local_util::geneMassShift(mass_shift_ptr, shift_mass,
                                             MassShiftType::UNEXPECTED);

  MassShiftPtrVec new_mass_shift_vec;
  new_mass_shift_vec.push_back(mass_shift_ptr);
  new_mass_shift_vec.insert(new_mass_shift_vec.end(),
                            fix_mass_shift_vec.begin(), fix_mass_shift_vec.end());
  std::sort(new_mass_shift_vec.begin(), new_mass_shift_vec.end(), MassShift::cmpPosInc);

  return proteoform_factory::geneProteoform(proteoform,
                                            ori_start + n_vec[idx],
                                            ori_end + c_vec[idx],
                                            new_mass_shift_vec,
                                            mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
}

ProteoformPtr LocalProcessor::twoPtmTermAdjust(ProteoformPtr proteoform, int num_match,
                                               const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                                               double & mass1, double & mass2) {
  int n_trunc_min, n_trunc_max, c_trunc_min, c_trunc_max;
  getNtermTruncRange(proteoform, extend_ms_ptr_vec, n_trunc_min, n_trunc_max);
  getCtermTruncRange(proteoform, extend_ms_ptr_vec, c_trunc_min, c_trunc_max);

  MassShiftPtrVec ori_mass_shift_vec = proteoform->getMassShiftPtrVec();

  MassShiftPtrVec expected_shift_vec
      = local_util::massShiftFilter(ori_mass_shift_vec, MassShiftType::UNEXPECTED);

  MassShiftPtr mass_shift1 = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0];
  MassShiftPtr mass_shift2 = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[1];

  double err = prec_mass * ppo_;

  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();

  double ori_mass1 = mass_shift1->getMassShift();
  double ori_mass2 = mass_shift2->getMassShift();
  mass1 = ori_mass1;
  mass2 = ori_mass2;

  std::vector<bool> ptm1_known_vec, ptm2_known_vec;
  std::vector<int> c_vec, n_vec;
  std::string n_seq, c_seq;
  std::vector<double> raw_scr_vec;

  for (int i = n_trunc_min; i <= n_trunc_max; i++) {
    for (int j = c_trunc_min; j <= c_trunc_max; j++) {
      if (i < 0) {
        n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + i, -i);
        mass1 = ori_mass1 - residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      } else {
        n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, i);
        mass1 = ori_mass1 + residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }
      if (j >= 0) {
        c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, j);
        mass2 = ori_mass2 - residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      } else {
        c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + j, -j);
        mass2 = ori_mass2 + residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }

      if ((mass1 > mng_ptr_->max_ptm_mass_ || mass1 < mng_ptr_->min_ptm_mass_
           || mass2 > mng_ptr_->max_ptm_mass_ || mass2 < mng_ptr_->min_ptm_mass_) && (i != 0 || j != 0))
        continue;

      n_vec.push_back(i);
      c_vec.push_back(j);

      PtmPairVec ptm_pair_vec = local_util::getPtmPairVecByMass(mass1, mass2, err, ptm_pair_vec_);
      if (ptm_pair_vec.size() == 0) {
        ptm_pair_vec.push_back(std::make_pair(nullptr, nullptr));
      }

      ptm1_known_vec.push_back(local_util::getPtmPtrVecByMass(mass1, err, ptm_vec_).size() > 0);
      ptm2_known_vec.push_back(local_util::getPtmPtrVecByMass(mass2, err, ptm_vec_).size() > 0);

      mass_shift1 = local_util::geneMassShift(mass_shift1, mass1, MassShiftType::UNEXPECTED);
      mass_shift2 = local_util::geneMassShift(mass_shift2, mass2, MassShiftType::UNEXPECTED);

      MassShiftPtrVec new_mass_shift_vec;
      new_mass_shift_vec.push_back(mass_shift1);
      new_mass_shift_vec.push_back(mass_shift2);
      new_mass_shift_vec.insert(new_mass_shift_vec.end(),
                                expected_shift_vec.begin(), expected_shift_vec.end());

      std::sort(new_mass_shift_vec.begin(), new_mass_shift_vec.end(), MassShift::cmpPosInc);

      proteoform = proteoform_factory::geneProteoform(proteoform,
                                                      ori_start + i,
                                                      ori_end + j,
                                                      new_mass_shift_vec,
                                                      mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      double raw_scr;
      compTwoPtmScr(proteoform, num_match, extend_ms_ptr_vec, prec_mass, raw_scr, ptm_pair_vec);
      raw_scr_vec.push_back(raw_scr);
    }
  }

  for (size_t i = 0; i < raw_scr_vec.size(); i++) {
    if (ptm1_known_vec[i] && ptm2_known_vec[i]) {
      raw_scr_vec[i] = raw_scr_vec[i] * mng_ptr_->theta_ * mng_ptr_->theta_;
    } else if (ptm1_known_vec[i] || ptm2_known_vec[i]) {
      raw_scr_vec[i] = raw_scr_vec[i] * mng_ptr_->theta_ * (1 - mng_ptr_->theta_);
    } else {
      raw_scr_vec[i] = raw_scr_vec[i] * (1 - mng_ptr_->theta_) * (1 - mng_ptr_->theta_);
    }
  }

  int idx = std::distance(raw_scr_vec.begin(), std::max_element(raw_scr_vec.begin(), raw_scr_vec.end()));

  if (n_vec[idx] < 0) {
    n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + n_vec[idx], -n_vec[idx]);
    mass1 = ori_mass1 - residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    n_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, n_vec[idx]);
    mass1 = ori_mass1 + residue_util::compResiduePtrVecMass(n_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  }

  if (c_vec[idx] >= 0) {
    c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, c_vec[idx]);
    mass2 = ori_mass2 - residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  } else {
    c_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + c_vec[idx], -c_vec[idx]);
    mass2 = ori_mass2 + residue_util::compResiduePtrVecMass(c_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  }

  mass_shift1 = local_util::geneMassShift(mass_shift1, mass1, MassShiftType::UNEXPECTED);
  mass_shift2 = local_util::geneMassShift(mass_shift2, mass2, MassShiftType::UNEXPECTED);

  MassShiftPtrVec new_mass_shift_vec;
  new_mass_shift_vec.push_back(mass_shift1);
  new_mass_shift_vec.push_back(mass_shift2);
  new_mass_shift_vec.insert(new_mass_shift_vec.end(),
                            expected_shift_vec.begin(), expected_shift_vec.end());

  std::sort(new_mass_shift_vec.begin(), new_mass_shift_vec.end(), MassShift::cmpPosInc);

  return proteoform_factory::geneProteoform(proteoform,
                                            ori_start + n_vec[idx],
                                            ori_end + c_vec[idx],
                                            new_mass_shift_vec,
                                            mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
}

void LocalProcessor::getNtermTruncRange(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                        int & min, int & max) {
  int left_sup, tmp;
  MassShiftPtr mass_shift_ptr = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0];
  local_util::compSupPeakNum(proteoform, extend_ms_ptr_vec, mass_shift_ptr, mng_ptr_->min_mass_, left_sup, tmp);

  if (left_sup > LEFT_SUP_LIMIT) {
    min = max = 0;
    return;
  }

  double ori_mass = mass_shift_ptr->getMassShift();
  double mass = ori_mass;
  int ori_start = proteoform->getStartPos();
  int ori_end = proteoform->getEndPos();
  max = min = 0;

  while (ori_start + min > 0) {
    min--;
    std::string t_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start + min, -min);
    mass = ori_mass - residue_util::compResiduePtrVecMass(t_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    if (mass <= mng_ptr_->min_ptm_mass_ || ori_start + min <= 0) {
      min++;
      break;
    }
  }

  mass = ori_mass;
  while (ori_start + max < ori_end) {
    max++;
    std::string t_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_start, max);
    mass = ori_mass + residue_util::compResiduePtrVecMass(t_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    if (mass >= mng_ptr_->max_ptm_mass_) {
      max--;
      break;
    }
  }
}

void LocalProcessor::getCtermTruncRange(ProteoformPtr proteoform, const ExtendMsPtrVec & extend_ms_ptr_vec,
                                        int & min, int & max) {
  MassShiftPtrVec mass_shift_vec = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED);
  MassShiftPtr mass_shift_ptr = mass_shift_vec[mass_shift_vec.size() - 1];

  int right_sup, tmp;
  local_util::compSupPeakNum(proteoform, extend_ms_ptr_vec, mass_shift_ptr, mng_ptr_->min_mass_, tmp, right_sup);

  if (right_sup > RIGHT_SUP_LIMIT) {
    min = max = 0;
    return;
  }

  double ori_mass = mass_shift_ptr->getMassShift();
  double mass = ori_mass;
  int ori_end = proteoform->getEndPos();
  int raw_seq_len = static_cast<int>(proteoform->getFastaSeqPtr()->getRawSeq().length());

  max = min = 0;

  while (ori_end + max < raw_seq_len - 1) {
    max++;
    std::string t_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1, max);
    mass = ori_mass - residue_util::compResiduePtrVecMass(t_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    if (mass < mng_ptr_->min_ptm_mass_ || ori_end + max >= raw_seq_len - 1) {
      max--;
      break;
    }
  }

  mass = ori_mass;
  while (ori_mass < mng_ptr_->max_ptm_mass_) {
    min--;
    std::string t_seq = proteoform->getFastaSeqPtr()->getRawSeq().substr(ori_end + 1 + min, -min);
    mass = ori_mass + residue_util::compResiduePtrVecMass(t_seq, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
    if (mass >= mng_ptr_->max_ptm_mass_) {
      min++;
      break;
    }
  }
}

// similar to processOneKnownPtm, we might get a nullptr from this function
ProteoformPtr LocalProcessor::processTwoKnownPtm(PrsmPtr prsm) {
  ProteoformPtr ori_prot_form = prsm->getProteoformPtr();

  MassShiftPtrVec ori_mass_shift_vec = ori_prot_form->getMassShiftPtrVec();

  MassShiftPtrVec unexpected_shift_vec
      = ori_prot_form->getMassShiftPtrVec(MassShiftType::UNEXPECTED);

  MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec,
                                                                   MassShiftType::UNEXPECTED);

  MassShiftPtr unexpected_shift1 = local_util::geneMassShift(unexpected_shift_vec[0],
                                                             unexpected_shift_vec[0]->getMassShift(),
                                                             MassShiftType::UNEXPECTED);

  MassShiftPtr unexpected_shift2;

  if (unexpected_shift_vec.size() == 2) {
    unexpected_shift2 = local_util::geneMassShift(unexpected_shift_vec[1],
                                                  unexpected_shift_vec[1]->getMassShift(),
                                                  MassShiftType::UNEXPECTED);
  } else {
    unexpected_shift2 = local_util::geneMassShift(unexpected_shift_vec[0],
                                                  0.0, MassShiftType::UNEXPECTED);
  }

  double shift_mass1 = unexpected_shift1->getMassShift();

  double shift_mass2 = unexpected_shift2->getMassShift();

  expected_shift_vec.push_back(unexpected_shift1);

  expected_shift_vec.push_back(unexpected_shift2);

  std::sort(expected_shift_vec.begin(), expected_shift_vec.end(), MassShift::cmpPosInc);

  ProteoformPtr two_shift_proteoform
      = std::make_shared<Proteoform>(ori_prot_form->getFastaSeqPtr(),
                                     ori_prot_form->getProtModPtr(),
                                     ori_prot_form->getStartPos(),
                                     ori_prot_form->getEndPos(),
                                     ori_prot_form->getResSeqPtr(),
                                     expected_shift_vec);

  double err = prsm->getAdjustedPrecMass() * ppo_;

  PtmPairVec ptm_pair_vec = local_util::getPtmPairVecByMass(shift_mass1, shift_mass2,
                                                            err, ptm_pair_vec_);

  if (ptm_pair_vec.size() == 0) {
    if (two_shift_proteoform->getProtModPtr()->getType() == ProtModBase::getType_NME_ACETYLATION()
        || two_shift_proteoform->getProtModPtr()->getType() == ProtModBase::getType_M_ACETYLATION()) {
      MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec,
                                                                       MassShiftType::UNEXPECTED);
      for (size_t k = 0; k < expected_shift_vec.size(); k++) {
        if (expected_shift_vec[k]->getTypePtr() == MassShiftType::PROTEIN_VARIABLE) {
          shift_mass1 += expected_shift_vec[k]->getMassShift();
          expected_shift_vec.erase(expected_shift_vec.begin() + k);
          break;
        }
      }

      unexpected_shift1 = local_util::geneMassShift(unexpected_shift1, shift_mass1,
                                                    MassShiftType::UNEXPECTED);

      expected_shift_vec.push_back(unexpected_shift1);

      expected_shift_vec.push_back(unexpected_shift2);

      std::sort(expected_shift_vec.begin(), expected_shift_vec.end(), MassShift::cmpPosInc);

      two_shift_proteoform = std::make_shared<Proteoform>(ori_prot_form->getFastaSeqPtr(),
                                                          ori_prot_form->getProtModPtr(),
                                                          ori_prot_form->getStartPos(),
                                                          ori_prot_form->getEndPos(),
                                                          ori_prot_form->getResSeqPtr(),
                                                          expected_shift_vec);
    }

    twoPtmTermAdjust(two_shift_proteoform, prsm->getMatchPeakNum(), prsm->getRefineMsPtrVec(),
                     prsm->getAdjustedPrecMass(), shift_mass1, shift_mass2);

    ptm_pair_vec = local_util::getPtmPairVecByMass(shift_mass1, shift_mass2, err, ptm_pair_vec_);
  }

  if (ptm_pair_vec.size() == 0) return nullptr;

  double raw_scr;
  ExtendMsPtrVec extend_ms_ptr_vec = prsm->getRefineMsPtrVec();

  compTwoPtmScr(two_shift_proteoform, prsm->getMatchPeakNum(), extend_ms_ptr_vec,
                prsm->getAdjustedPrecMass(),
                raw_scr, ptm_pair_vec);

  // can be explained by a variable ptm, but no modifiable site
  if (raw_scr == 0) return nullptr;

  two_shift_proteoform->setProteoClusterId(ori_prot_form->getProteoClusterId());
  two_shift_proteoform->setProtId(ori_prot_form->getProtId());

  PtmPtr ptm1 = ptm_pair_vec[0].first;
  PtmPtr ptm2 = ptm_pair_vec[0].second;
  raw_scr = raw_scr * theta_ * theta_ * (1 - beta_);

  std::vector<double> empty_scr_vec;
  LocalAnnoPtr anno1 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm1);
  two_shift_proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0]->getChangePtr(0)->setLocalAnno(anno1);

  LocalAnnoPtr anno2 = std::make_shared<LocalAnno>(0, 0, 0, empty_scr_vec, raw_scr, ptm2);
  two_shift_proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[1]->getChangePtr(0)->setLocalAnno(anno2);

  return compSplitPoint(two_shift_proteoform, prsm->getMatchPeakNum(),
                        extend_ms_ptr_vec, prsm->getAdjustedPrecMass());
}

PrsmPtr LocalProcessor::processTwoPtm(PrsmPtr prsm) {
  int ori_num_match_ion = local_util::compMatchFragNum(prsm->getProteoformPtr(),
                                                       prsm->getRefineMsPtrVec(),
                                                       mng_ptr_->min_mass_);

  ProteoformPtr two_known_prsm = processTwoKnownPtm(prsm);

  if (two_known_prsm != nullptr) {
    int new_num_match_ion = local_util::compMatchFragNum(two_known_prsm,
                                                         prsm->getRefineMsPtrVec(),
                                                         mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - DESC_MATCH_LIMIT) {
      prsm->setProteoformPtr(two_known_prsm, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  ProteoformPtr one_known_prsm = processOneKnownPtm(prsm);

  if (one_known_prsm != nullptr) {
    double new_num_match_ion = local_util::compMatchFragNum(one_known_prsm,
                                                            prsm->getRefineMsPtrVec(),
                                                            mng_ptr_->min_mass_);
    if (new_num_match_ion > ori_num_match_ion - DESC_MATCH_LIMIT) {
      prsm->setProteoformPtr(one_known_prsm, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
      return prsm;
    }
  }

  return prsm;
}

bool LocalProcessor::modifiable(ProteoformPtr proteoform_ptr, int i, PtmPtr ptm_ptr) {
  if (ptm_ptr == nullptr) return true;

  MassShiftPtrVec fixed_shift_vec = proteoform_ptr->getMassShiftPtrVec(MassShiftType::FIXED);

  for (size_t k = 0; k < fixed_shift_vec.size(); k++) {
    if (fixed_shift_vec[k]->getLeftBpPos() <= i && i < fixed_shift_vec[k]->getRightBpPos()) {
      return false;
    }
  }

  int start = proteoform_ptr->getStartPos();
  int end = proteoform_ptr->getEndPos();

  ResiduePtr residue_ptr = proteoform_ptr->getResSeqPtr()->getResiduePtr(i);

  ModPtrVec mod_list;

  if (i == 0 && start == 0) {
    mod_list = mod_list_N_;
  } else if (i + start == end) {
    mod_list = mod_list_C_;
  } else {
    mod_list = mod_list_any_;
  }

  for (size_t j = 0; j < mod_list.size(); j++) {
    if (mod_list[j]->getOriResiduePtr()->isSame(residue_ptr) &&
        mod_list[j]->getModResiduePtr()->getPtmPtr()->isSame(ptm_ptr))
      return true;
  }

  return false;
}

void LocalProcessor::compOnePtmScr(ProteoformPtr proteoform,
                                   const ExtendMsPtrVec & extend_ms_ptr_vec,
                                   std::vector<double> &scr_vec, double & raw_scr,
                                   PtmPtrVec & ptm_vec) {
  raw_scr = 0.0;

  MassShiftPtr mass_shift = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0];

  if (ptm_vec.size() == 0) ptm_vec.push_back(nullptr);

  std::vector<double> temp;
  std::vector<std::vector<double> > scr_vec2d;
  int n = proteoform->getLen();
  int count = 0;

  for (size_t i = 0; i < ptm_vec.size(); i++) {
    scr_vec.clear();
    for (int j = 0; j < n; j++) {
      if (modifiable(proteoform, j, ptm_vec[i])) {
        count++;
        mass_shift->setLeftBpPos(j);
        mass_shift->setRightBpPos(j + 1);
        int match = static_cast<int>(local_util::compMatchFragNum(proteoform,
                                                                  extend_ms_ptr_vec,
                                                                  mng_ptr_->min_mass_));

        scr_vec.push_back(std::pow(p1_, n - match) * std::pow(p2_, match));
      } else {
        scr_vec.push_back(0.0);
      }
    }
    scr_vec2d.push_back(scr_vec);
    temp.push_back(std::accumulate(scr_vec.begin(), scr_vec.end(), 0.0));
  }

  scr_vec.clear();
  int idx = std::distance(temp.begin(), std::max_element(temp.begin(), temp.end()));

  if (temp[idx] == 0) return;

  raw_scr = std::accumulate(scr_vec2d[idx].begin(), scr_vec2d[idx].end(), 0.0) / count;
  scr_vec = scr_vec2d[idx];
  local_util::normalize(scr_vec);

  PtmPtr p = ptm_vec[idx];
  ptm_vec.clear();
  ptm_vec.push_back(p);
}

void LocalProcessor::compTwoPtmScr(ProteoformPtr proteoform, int num_match,
                                   const ExtendMsPtrVec & extend_ms_ptr_vec, double prec_mass,
                                   double & raw_scr, PtmPairVec & ptm_pair_vec) {
  MassShiftPtrVec ori_mass_shift_vec = proteoform->getMassShiftPtrVec();
  double shift_mass1 = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0]->getMassShift();
  double shift_mass2 = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[1]->getMassShift();
  ProteoformPtr no_shift_proteoform
      = proteoform_factory::geneProteoform(proteoform,
                                           proteoform->getStartPos(),
                                           proteoform->getEndPos(),
                                           local_util::massShiftFilter(ori_mass_shift_vec, MassShiftType::UNEXPECTED),
                                           mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::vector<double> scr_vec;

  for (size_t k = 0; k < ptm_pair_vec.size(); k++) {
    local_util::ptmMassAdjust(shift_mass1, shift_mass2, ptm_pair_vec[k].first, ptm_pair_vec[k].second);
    scr_vec.push_back(dpTwoPtmScr(no_shift_proteoform, num_match, extend_ms_ptr_vec, prec_mass,
                                  shift_mass1, shift_mass2, ptm_pair_vec[k].first, ptm_pair_vec[k].second));
  }

  int idx = std::distance(scr_vec.begin(), std::max_element(scr_vec.begin(), scr_vec.end()));
  raw_scr = scr_vec[idx];
  PtmPtr p1 = ptm_pair_vec[idx].first;
  PtmPtr p2 = ptm_pair_vec[idx].second;
  ptm_pair_vec.clear();
  ptm_pair_vec.push_back(std::make_pair(p1, p2));
}

double LocalProcessor::dpTwoPtmScr(ProteoformPtr proteoform, int h,
                                   const ExtendMsPtrVec & extend_ms_ptr_vec,
                                   double prec_mass, double mass1, double mass2, PtmPtr ptm1, PtmPtr ptm2) {
  int g = proteoform->getLen();
  double scr = 0.0;
  int count = 0;

  for (size_t k = 0; k < extend_ms_ptr_vec.size(); k++) {
    std::vector<std::vector<double>> b_table(3);
    b_table[0] = local_util::geneNTheoMass(proteoform, extend_ms_ptr_vec[k], mng_ptr_->min_mass_);
    local_util::fillTableB(b_table, mass1, mass2);

    std::vector<std::vector<int>> s_table(3);
    local_util::fillTableS(b_table, s_table, extend_ms_ptr_vec[k], prec_mass);

    // fill D(f,g,h)
    int d_table[3][g + 1][h + 1];
    memset(d_table, 0, sizeof(int) * 3 * (g + 1) * (h + 1));
    d_table[0][0][0] = 1;

    for (int i = 1; i <= g; i++) {
      for (int j = 0; j <= h; j++) {
        if (j >= s_table[0][i - 1]) {
          d_table[0][i][j] = d_table[0][i - 1][j - s_table[0][i - 1]];
        } else {
          d_table[0][i][j] = 0;
        }
      }
    }

    for (int i = 1; i <= g; i++) {
      for (int j = 0; j <= h; j++) {
        if (modifiable(proteoform, i - 1, ptm1) && j >= s_table[1][i - 1]) {
          d_table[1][i][j] = d_table[0][i - 1][j - s_table[1][i - 1]] + d_table[1][i - 1][j - s_table[1][i - 1]];
        } else if (j >= s_table[1][i - 1]) {
          d_table[1][i][j] = d_table[1][i - 1][j - s_table[1][i - 1]];
        } else {
          d_table[1][i][j] = 0;
        }
      }
    }

    for (int i = 1; i <= g; i++) {
      for (int j = 0; j <= h; j++) {
        if (modifiable(proteoform, i - 1, ptm2) && j >= s_table[2][i - 1]) {
          d_table[2][i][j] = d_table[1][i - 1][j - s_table[2][i - 1]] + d_table[2][i - 1][j - s_table[2][i - 1]];
        } else if (j >= s_table[2][i - 1]) {
          d_table[2][i][j] = d_table[2][i - 1][j - s_table[2][i - 1]];
        } else {
          d_table[2][i][j] = 0;
        }
      }
    }

    for (int i = 0; i <= h; i++) {
      count += d_table[2][g][i];
      scr += d_table[2][g][i] * std::pow(p1_, g - i) * std::pow(p2_, i);
    }
  }

  return scr / count;
}

ProteoformPtr LocalProcessor::compSplitPoint(ProteoformPtr proteoform, int h,
                                             const ExtendMsPtrVec & extend_ms_ptr_vec,
                                             double prec_mass) {
  MassShiftPtrVec ori_mass_shift_vec = proteoform->getMassShiftPtrVec();

  MassShiftPtr shift_ptr1 = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[0];
  MassShiftPtr shift_ptr2 = proteoform->getMassShiftPtrVec(MassShiftType::UNEXPECTED)[1];

  double mass1 = shift_ptr1->getMassShift();
  double mass2 = shift_ptr2->getMassShift();

  PtmPtr ptm1 = shift_ptr1->getChangePtr(0)->getLocalAnno()->getPtmPtr();
  PtmPtr ptm2 = shift_ptr2->getChangePtr(0)->getLocalAnno()->getPtmPtr();

  local_util::ptmMassAdjust(mass1, mass2, ptm1, ptm2);

  shift_ptr1 = local_util::geneMassShift(shift_ptr1, mass1, MassShiftType::UNEXPECTED);

  shift_ptr2 = local_util::geneMassShift(shift_ptr2, mass2, MassShiftType::UNEXPECTED);

  MassShiftPtrVec expected_shift_vec = local_util::massShiftFilter(ori_mass_shift_vec, MassShiftType::UNEXPECTED);

  expected_shift_vec.push_back(shift_ptr1);

  expected_shift_vec.push_back(shift_ptr2);

  int prot_cluster_id = proteoform->getProteoClusterId();

  int prot_id = proteoform->getProtId();

  proteoform = proteoform_factory::geneProteoform(proteoform,
                                                  proteoform->getStartPos(),
                                                  proteoform->getEndPos(),
                                                  expected_shift_vec,
                                                  mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  proteoform->setProteoClusterId(prot_cluster_id);

  proteoform->setProtId(prot_id);

  int g = proteoform->getLen();

  ProteoformPtr no_shift_proteoform
      = proteoform_factory::geneProteoform(proteoform,
                                           proteoform->getStartPos(),
                                           proteoform->getEndPos(),
                                           local_util::massShiftFilter(ori_mass_shift_vec, MassShiftType::UNEXPECTED),
                                           mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  std::vector<double> split_scr_vec;
  for (int k = 1; k < g; k++) {
    double scr = 0.0;
    for (size_t i = 0; i < extend_ms_ptr_vec.size(); i++) {
      std::vector<std::vector<double>> b_table(3);
      b_table[0] = local_util::geneNTheoMass(no_shift_proteoform, extend_ms_ptr_vec[i], mng_ptr_->min_mass_);
      local_util::fillTableB(b_table, mass1, mass2);

      std::vector<std::vector<int>> s_table(3);
      local_util::fillTableS(b_table, s_table, extend_ms_ptr_vec[i], prec_mass);

      int d_table[3][g + 1][h + 1];

      memset(d_table, 0, sizeof(int) * 3 * (g + 1) * (h + 1));
      d_table[0][0][0] = 1;

      for (int i = 1; i <= g; i++) {
        for (int j = 0; j <= h; j++) {
          if (j >= s_table[0][i - 1]) {
            d_table[0][i][j] = d_table[0][i - 1][j - s_table[0][i - 1]];
          } else {
            d_table[0][i][j] = 0;
          }
        }
      }

      for (int i = 1; i <= g; i++) {
        for (int j = 0; j <= h; j++) {
          if (modifiable(proteoform, i - 1, ptm1) && j >= s_table[1][i - 1] && i <= k) {
            d_table[1][i][j] = d_table[0][i - 1][j - s_table[1][i - 1]] + d_table[1][i - 1][j - s_table[1][i - 1]];
          } else if (j >= s_table[1][i - 1]) {
            d_table[1][i][j] = d_table[1][i - 1][j - s_table[1][i - 1]];
          } else {
            d_table[1][i][j] = 0;
          }
        }
      }

      for (int i = 1; i <= g; i++) {
        for (int j = 0; j <= h; j++) {
          if (modifiable(proteoform, i - 1, ptm2) && j >= s_table[2][i - 1] && i > k) {
            d_table[2][i][j] = d_table[1][i - 1][j - s_table[2][i - 1]] + d_table[2][i - 1][j - s_table[2][i - 1]];
          } else if (j >= s_table[2][i - 1]) {
            d_table[2][i][j] = d_table[2][i - 1][j - s_table[2][i - 1]];
          } else {
            d_table[2][i][j] = 0;
          }
        }
      }

      for (int i = 0; i <= h; i++) {
        scr += d_table[2][g][i] * std::pow(p1_, i) * std::pow(p2_, i);
      }
    }
    split_scr_vec.push_back(scr);
  }

  int split_point = std::distance(split_scr_vec.begin(), std::max_element(split_scr_vec.begin(), split_scr_vec.end()));

  double split_max = *std::max_element(split_scr_vec.begin(), split_scr_vec.end());

  double split_scr = 0.0;

  int split_end = split_scr_vec.size();

  for (; split_end > 0; split_end--) {
    if (split_scr_vec[split_end] == split_max)
      break;
  }

  for (int i = split_point; i <= split_end; i++) {
    split_scr += split_scr_vec[i];
  }

  split_point = (split_point + split_end) / 2;

  split_scr = split_scr / std::accumulate(split_scr_vec.begin(), split_scr_vec.end(), 0.0);

  if (split_scr <= mng_ptr_->threshold_) {
    return nullptr;
  }

  std::vector<double> ptm_scr;

  for (int i = 0; i <= split_point; i++) {
    shift_ptr1->setLeftBpPos(i);
    shift_ptr1->setRightBpPos(i + 1);
    if (modifiable(proteoform, i, ptm1)) {
      int match = static_cast<int>(local_util::compMatchFragNum(proteoform,
                                                                extend_ms_ptr_vec,
                                                                mng_ptr_->min_mass_));
      ptm_scr.push_back(std::pow(p1_, split_point - match) * std::pow(p2_, match));
    } else {
      ptm_scr.push_back(0.0);
    }
  }

  local_util::normalize(ptm_scr);
  int bgn, end;
  double conf;
  std::transform(ptm_scr.begin(), ptm_scr.end(), ptm_scr.begin(),
                 std::bind1st(std::multiplies<double>(), split_scr));
  local_util::scrFilter(ptm_scr, bgn, end, conf, mng_ptr_->threshold_);

  if (bgn == -1) {
    return nullptr;
  } else {
    LocalAnnoPtr anno1 = std::make_shared<LocalAnno>(bgn, end, conf, ptm_scr, 0, ptm1);
    shift_ptr1->getChangePtr(0)->setLocalAnno(anno1);
    shift_ptr1->setLeftBpPos(anno1->getLeftBpPos());
    shift_ptr1->setRightBpPos(anno1->getRightBpPos() + 1);
  }

  ptm_scr.clear();
  int len = proteoform->getLen();
  for (int i = split_point + 1; i < len; i++) {
    shift_ptr2->setLeftBpPos(i);
    shift_ptr2->setRightBpPos(i + 1);
    if (modifiable(proteoform, i, ptm2)) {
      int match = static_cast<int>(local_util::compMatchFragNum(proteoform,
                                                                extend_ms_ptr_vec,
                                                                mng_ptr_->min_mass_));
      ptm_scr.push_back(std::pow(p1_, len - match) * std::pow(p2_, match));
    } else {
      ptm_scr.push_back(0.0);
    }
  }

  local_util::normalize(ptm_scr);

  std::transform(ptm_scr.begin(), ptm_scr.end(), ptm_scr.begin(),
                 std::bind1st(std::multiplies<double>(), split_scr));

  local_util::scrFilter(ptm_scr, bgn, end, conf, mng_ptr_->threshold_);

  if (bgn == -1) {
    return nullptr;
  } else {
    LocalAnnoPtr anno2 = std::make_shared<LocalAnno>(split_point + bgn + 1, split_point + end + 1, conf, ptm_scr, 0, ptm2);
    shift_ptr2->getChangePtr(0)->setLocalAnno(anno2);
    shift_ptr2->setLeftBpPos(anno2->getLeftBpPos());
    shift_ptr2->setRightBpPos(anno2->getRightBpPos() + 1);
  }

  return proteoform;
}

}  // namespace prot

