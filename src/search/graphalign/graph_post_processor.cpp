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

#include "search/graphalign/graph_post_processor.hpp"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/base/mod_base.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/base/prot_mod_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "search/graphalign/graph_align_mng.hpp"
#include "seq/mass_shift.hpp"
#include "seq/proteoform_factory.hpp"
#include "seq/residue_seq.hpp"

namespace toppic {

std::vector<double> mass_split(double mass,
                               const std::vector<PtmPtr>& ptm_vec) {
  std::vector<double> mass_vec(ptm_vec.size());
  double ptm_mass = 0.0;

  for (size_t k = 0; k < ptm_vec.size(); k++) {
    ptm_mass += ptm_vec[k]->getMonoMass();
  }

  double err = (mass - ptm_mass) / ptm_vec.size();

  for (size_t k = 0; k < ptm_vec.size(); k++) {
    mass_vec[k] = ptm_vec[k]->getMonoMass() + err;
  }

  return mass_vec;
}

void GraphPostProcessor::geneMassPtmMap() {
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;

  // The side-chain variable PTMs of the file (Position any) are the names of
  // the variable shifts. The acetylation of a protein N-terminal form and a
  // proteoform N-terminal PTM (Position N-term) are not (chooseNTermForm):
  // an N-term line used to name side-chain shifts too, and with the same
  // name as a side-chain line made the file unreadable here.
  ModPtrVec2D mod_ptr_vec2d = mod_util::readModTxt(var_mod_file_name);
  std::vector<PtmPtr> ptm_vec;
  for (size_t i = 0; i < mod_ptr_vec2d[2].size(); i++) {
    PtmPtr ptm_ptr = mod_ptr_vec2d[2][i]->getModResiduePtr()->getPtmPtr();
    if (std::find(ptm_vec.begin(), ptm_vec.end(), ptm_ptr) == ptm_vec.end()) {
      ptm_vec.push_back(ptm_ptr);
    }
  }
  n_term_mod_ptrs_ = mod_ptr_vec2d[3];

  std::sort(ptm_vec.begin(), ptm_vec.end(),
            [](const PtmPtr& a, const PtmPtr& b) {
              return a->getName() < b->getName();
            });

  std::vector<PtmPtr>::iterator it = std::unique(
      ptm_vec.begin(), ptm_vec.end(), [](const PtmPtr& a, const PtmPtr& b) {
        return a->getName() == b->getName();
      });
  ptm_vec.erase(it, ptm_vec.end());

  int inte_tole = mng_ptr_->getIntTolerance();

  for (size_t i = 0; i < ptm_vec.size(); i++) {
    int mass = std::ceil(ptm_vec[i]->getMonoMass() * mng_ptr_->convert_ratio_);
    bool found = false;
    for (auto it = mass_ptm_map_.begin(); it != mass_ptm_map_.end(); it++) {
      if (std::abs(it->first - mass) <= inte_tole) {
        found = true;
        break;
      }
    }
    if (!found) {
      mass_ptm_map_[mass].push_back(ptm_vec[i]);
    }
  }

  for (int k = 2; k <= mng_ptr_->max_known_mods_; k++) {
    std::vector<int> cur_map_mass;
    for (auto it = mass_ptm_map_.begin(); it != mass_ptm_map_.end(); it++) {
      cur_map_mass.push_back(it->first);
    }

    for (size_t i = 0; i < ptm_vec.size(); i++) {
      int mass =
          std::ceil(ptm_vec[i]->getMonoMass() * mng_ptr_->convert_ratio_);
      for (size_t cnt = 0; cnt < cur_map_mass.size(); cnt++) {
        int new_mass = mass + cur_map_mass[cnt];
        bool found = false;
        auto it2 = mass_ptm_map_.begin();
        while (it2 != mass_ptm_map_.end() && !found) {
          if (std::abs(it2->first - new_mass) <= inte_tole) {
            found = true;
          }
          it2++;
        }
        if (!found) {
          mass_ptm_map_[new_mass].insert(
              mass_ptm_map_[new_mass].end(),
              mass_ptm_map_[cur_map_mass[cnt]].begin(),
              mass_ptm_map_[cur_map_mass[cnt]].end());
          mass_ptm_map_[new_mass].push_back(ptm_vec[i]);
        }
      }
    }
  }
}

// Each shift takes the combination of mass_ptm_map_ nearest to its mass.
// When these exceed max_known_mods_ (--var-ptm) PTMs in all, the combinations
// are chosen together instead: the smallest sum of mass differences with at
// most max_known_mods_ PTMs, by a dynamic programme over the shifts and the
// number of PTMs used (for each shift and PTM number, the nearest combination
// of that number). Below the default 0.1 Da alignment tolerance the map keeps
// more multi-PTM combinations, and the nearest ones gave proteoforms up to 9
// variable PTMs of 5, which overran the MCMC E-value tables (indexed by the
// variable PTM number; toppic_test experiment 114). Choosing each shift's
// combination in turn within the PTMs left kept the limit but labelled a
// 70.06 Da shift with one Phospho; the joint choice explains it as Acetyl and
// two Methyl. While the nearest combinations keep the limit, as every PrSM
// did at 0.1 and 0.5 Da, they are the choice.
std::vector<PtmPtrVec> GraphPostProcessor::choosePtmCombs(
    const MassShiftPtrVec& shift_vec) {
  int max_num = mng_ptr_->max_known_mods_;
  size_t shift_num = shift_vec.size();
  std::vector<PtmPtrVec> nearest(shift_num);
  // best_err[k][c] and best_comb[k][c]: the nearest combination of c PTMs to
  // shift k (-1: none)
  std::vector<std::vector<int>> best_err(shift_num,
                                         std::vector<int>(max_num + 1, -1));
  std::vector<std::vector<const PtmPtrVec*>> best_comb(
      shift_num, std::vector<const PtmPtrVec*>(max_num + 1, nullptr));
  int total_num = 0;
  for (size_t k = 0; k < shift_num; k++) {
    int mass =
        std::ceil(shift_vec[k]->getMassShift() * mng_ptr_->convert_ratio_);
    int err = -1;
    for (auto it = mass_ptm_map_.begin(); it != mass_ptm_map_.end(); it++) {
      int e = std::abs(mass - it->first);
      if (err < 0 || e < err) {
        err = e;
        nearest[k] = it->second;
      }
      int c = static_cast<int>(it->second.size());
      if (c <= max_num && (best_err[k][c] < 0 || e < best_err[k][c])) {
        best_err[k][c] = e;
        best_comb[k][c] = &it->second;
      }
    }
    total_num += static_cast<int>(nearest[k].size());
  }
  if (total_num <= max_num) {
    return nearest;
  }
  // cost[k][u]: the smallest sum of mass differences of shifts 0..k-1 with u
  // PTMs (-1: not reachable); num[k][u]: the PTM number of shift k-1 there
  std::vector<std::vector<long long>> cost(
      shift_num + 1, std::vector<long long>(max_num + 1, -1));
  std::vector<std::vector<int>> num(shift_num + 1,
                                    std::vector<int>(max_num + 1, 0));
  cost[0][0] = 0;
  for (size_t k = 0; k < shift_num; k++) {
    for (int u = 0; u <= max_num; u++) {
      if (cost[k][u] < 0) {
        continue;
      }
      for (int c = 1; u + c <= max_num; c++) {
        if (best_err[k][c] < 0) {
          continue;
        }
        long long v = cost[k][u] + best_err[k][c];
        if (cost[k + 1][u + c] < 0 || v < cost[k + 1][u + c]) {
          cost[k + 1][u + c] = v;
          num[k + 1][u + c] = c;
        }
      }
    }
  }
  int best_u = -1;
  for (int u = 0; u <= max_num; u++) {
    if (cost[shift_num][u] >= 0 &&
        (best_u < 0 || cost[shift_num][u] < cost[shift_num][best_u])) {
      best_u = u;
    }
  }
  if (best_u < 0) {
    // more variable shifts than max_known_mods_, which the alignment does
    // not produce
    LOG_WARN("Variable shifts exceed the maximum number of variable PTMs.");
    return nearest;
  }
  std::vector<PtmPtrVec> combs(shift_num);
  int u = best_u;
  for (size_t k = shift_num; k > 0; k--) {
    int c = num[k][u];
    combs[k - 1] = *best_comb[k - 1][c];
    u -= c;
  }
  return combs;
}

int GraphPostProcessor::nearestPtmComb(int mass, int max_ptm_num,
                                       PtmPtrVec& ptm_vec) {
  int err = -1;
  ptm_vec.clear();
  for (auto it = mass_ptm_map_.begin(); it != mass_ptm_map_.end(); it++) {
    if (static_cast<int>(it->second.size()) > max_ptm_num) {
      continue;
    }
    if (err < 0 || std::abs(mass - it->first) < err) {
      err = std::abs(mass - it->first);
      ptm_vec = it->second;
    }
  }
  return err;
}

namespace {

double ptmMassSum(const PtmPtrVec& ptm_vec) {
  double mass = 0;
  for (size_t i = 0; i < ptm_vec.size(); i++) {
    mass += ptm_vec[i]->getMonoMass();
  }
  return mass;
}

// A mass shift of a proteoform that starts offset residues later: its
// positions moved by offset, a shift that began on a removed residue
// beginning on the new first residue.
MassShiftPtr moveShift(const MassShiftPtr& shift_ptr, int offset) {
  int left = std::max(0, shift_ptr->getLeftBpPos() - offset);
  int right = std::max(left, shift_ptr->getRightBpPos() - offset);
  MassShiftPtr new_shift_ptr =
      std::make_shared<MassShift>(left, right, shift_ptr->getMassShift());
  AlterPtrVec alter_ptrs = shift_ptr->getAlterPtrVec();
  AlterPtrVec new_alter_ptrs;
  for (size_t i = 0; i < alter_ptrs.size(); i++) {
    int alter_left = std::max(0, alter_ptrs[i]->getLeftBpPos() - offset);
    int alter_right =
        std::max(alter_left, alter_ptrs[i]->getRightBpPos() - offset);
    new_alter_ptrs.push_back(std::make_shared<Alter>(
        alter_left, alter_right, alter_ptrs[i]->getTypePtr(),
        alter_ptrs[i]->getMass(), alter_ptrs[i]->getModPtr()));
  }
  new_shift_ptr->setAlterPtrVec(new_alter_ptrs);
  return new_shift_ptr;
}

}  // namespace

// The protein N-terminal form of a PrSM, among those -n allows and allowMod
// accepts for its protein. The graph alignment's form edges are at residue 0
// (ProteoAnno::anno) and take a variable PTM layer, so an alignment that
// used one starts at the protein N-terminus and has a variable shift that
// begins at its first residue: its mass is the form's acetylation, alone or
// with variable PTMs of the residues after it. Start 0: NONE, M_ACETYLATION
// (the shift less the acetyl) or NME_ACETYLATION (the acetyl in place of the
// M: the shift less the acetyl plus the M; the proteoform starts at residue
// 1). Start 1: NME when allowed, else NONE (a truncation), or
// NME_ACETYLATION (the shift less the acetyl). Later starts: NONE, with or
// without a truncation N-terminal modification of the file's proteoform
// N-terminal PTMs (the shift less its mass), whose consistent pairs the
// graph alignment starts with (GraphAlign::getConsistentPairs).
//
// The shift named as it is (the nearest PTM combination, as choosePtmCombs
// names it while the proteoform keeps --var-ptm) is the form without an
// acetylation. An acetylated form leaves
// the rest of the shift: nothing, or the nearest PTM combination. It is
// taken when it explains the shift at least as well, with a mass error (in
// Da) within the alignment tolerance or no larger than the shift's own error
// up to twice the tolerance, and with fewer variable PTMs, or as many and a
// smaller error. The errors are compared in Da: the integer masses round the
// two explanations of one shift differently. The shift's own error allows
// for the diagonals' mass estimates, which miss a PTM's mass by more than
// the tolerance at times (0.10-0.15 Da for an acetylation).
GraphNTermForm GraphPostProcessor::chooseNTermForm(
    const ProteoformPtr& form_ptr, const MassShiftPtr& first_shift_ptr) {
  GraphNTermForm plain;
  plain.prot_mod_ptr = ProtModBase::getProtModPtr_NONE();
  int start = form_ptr->getStartPos();
  std::vector<GraphNTermForm> mod_forms;
  if (start >= 2) {
    // A truncated proteoform: its truncation N-terminal modifications (the
    // proteoform N-terminal PTMs of the file whose amino acid is its first
    // residue's), as Proteoform::setTruncNTermModPtr allows them.
    if (first_shift_ptr != nullptr && form_ptr->getResSeqPtr()->getLen() > 0) {
      AminoAcidPtr acid_ptr =
          form_ptr->getResSeqPtr()->getResiduePtr(0)->getAminoAcidPtr();
      for (size_t i = 0; i < n_term_mod_ptrs_.size(); i++) {
        const ModPtr& mod_ptr = n_term_mod_ptrs_[i];
        if (mod_ptr->getOriResiduePtr()->getAminoAcidPtr() != acid_ptr) {
          continue;
        }
        GraphNTermForm form;
        form.prot_mod_ptr = plain.prot_mod_ptr;
        form.trunc_n_term_mod_ptr = mod_ptr;
        form.takes_first_shift = true;
        form.remain_mass = first_shift_ptr->getMassShift() - mod_ptr->getShift();
        mod_forms.push_back(form);
      }
    }
  }
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  ProtModPtrVec prot_mod_ptrs;
  ResiduePtrVec residues;
  if (start <= 1) {
    prot_mod_ptrs = prsm_para_ptr->getProtModPtrVec();
    ProteoformPtr db_form_ptr = proteoform_factory::geneDbProteoformPtr(
        form_ptr->getFastaSeqPtr(), prsm_para_ptr->getFixModPtrVec());
    residues = db_form_ptr->getResSeqPtr()->getResidues();
  }
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    const ProtModPtr& mod_ptr = prot_mod_ptrs[i];
    if (!prot_mod_util::allowMod(mod_ptr, residues)) {
      continue;
    }
    std::string type = mod_ptr->getType();
    if (type == ProtModBase::getType_NME()) {
      if (start == 1) {
        plain.prot_mod_ptr = mod_ptr;
      }
      continue;
    }
    if (first_shift_ptr == nullptr) {
      continue;
    }
    double acetyl_mass = mod_ptr->getModPtr()->getShift();
    GraphNTermForm form;
    form.prot_mod_ptr = mod_ptr;
    form.takes_first_shift = true;
    if (type == ProtModBase::getType_M_ACETYLATION() && start == 0) {
      form.remain_mass = first_shift_ptr->getMassShift() - acetyl_mass;
    } else if (type == ProtModBase::getType_NME_ACETYLATION() && start == 0) {
      form.start_offset = 1;
      form.remain_mass = first_shift_ptr->getMassShift() -
                         (acetyl_mass - residues[0]->getMass());
    } else if (type == ProtModBase::getType_NME_ACETYLATION() && start == 1) {
      form.remain_mass = first_shift_ptr->getMassShift() - acetyl_mass;
    } else {
      continue;
    }
    mod_forms.push_back(form);
  }
  if (mod_forms.empty()) {
    return plain;
  }

  double tole = mng_ptr_->topmg_align_error_tolerance_;
  double shift_mass = first_shift_ptr->getMassShift();
  PtmPtrVec ptm_vec;
  nearestPtmComb(std::ceil(shift_mass * mng_ptr_->convert_ratio_),
                 mng_ptr_->max_known_mods_, ptm_vec);
  double plain_err = std::abs(shift_mass - ptmMassSum(ptm_vec));
  double max_err = std::max(tole, std::min(plain_err, 2 * tole));
  // a margin for the rounding of masses that are equal
  const double eps = 1e-6;
  GraphNTermForm best = plain;
  int best_ptm_num = static_cast<int>(ptm_vec.size());
  double best_err = plain_err;
  for (size_t i = 0; i < mod_forms.size(); i++) {
    GraphNTermForm form = mod_forms[i];
    double err;
    int ptm_num = 0;
    if (std::abs(form.remain_mass) <= max_err + eps) {
      form.drop_first_shift = true;
      err = std::abs(form.remain_mass);
    } else {
      if (nearestPtmComb(
              std::ceil(form.remain_mass * mng_ptr_->convert_ratio_),
              mng_ptr_->max_known_mods_, ptm_vec) < 0) {
        continue;
      }
      err = std::abs(form.remain_mass - ptmMassSum(ptm_vec));
      ptm_num = static_cast<int>(ptm_vec.size());
    }
    if (err > max_err + eps) {
      continue;
    }
    if (ptm_num < best_ptm_num ||
        (ptm_num == best_ptm_num && err < best_err - eps)) {
      best = form;
      best_ptm_num = ptm_num;
      best_err = err;
    }
  }
  return best;
}

void GraphPostProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();

  geneMassPtmMap();

  int spectrum_num =
      msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  std::string input_file_name =
      file_util::basename(sp_file_name) + "." + input_ext_;
  std::string output_file_name =
      file_util::basename(sp_file_name) + "." + output_ext_;

  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);
  PrsmXmlWriterPtr prsm_writer =
      std::make_shared<PrsmXmlWriter>(output_file_name);

  FastaIndexReaderPtr fasta_reader =
      std::make_shared<FastaIndexReader>(db_file_name);

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(
      fasta_reader, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr(),
      sp_para_ptr->getEnvCnnCutoff(),
        sp_para_ptr->getFragDedupPpo());

  int cnt = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() > 0) {
    cnt += group_spec_num;
    if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->containsPrec()) {
      double prec_mono_mass =
          deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
      SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
          deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        int spec_id = spec_set_ptr->getSpectrumId();
        while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
          DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
          double adjusted_prec_mass = prsm_ptr->getAdjustedPrecMass();
          ExtendMsPtrVec refine_ms_ptr_vec =
              extend_ms_factory::geneMsThreePtrVec(
                  deconv_ms_ptr_vec, sp_para_ptr, adjusted_prec_mass);
          PeakIonPairPtrVec pair_vec = peak_ion_pair_util::genePeakIonPairs(
              prsm_ptr->getProteoformPtr(), refine_ms_ptr_vec[0],
              sp_para_ptr->getMinMass());
          std::sort(pair_vec.begin(), pair_vec.end(),
                    PeakIonPair::cmpTheoPeakPosInc);

          ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
          MassShiftPtrVec shift_vec =
              form_ptr->getMassShiftPtrVec(AlterType::VARIABLE);

          // The protein N-terminal form first: an acetylated form takes its
          // acetylation from the variable shift that begins at the first
          // residue, and the rest of that shift is named below, or the shift
          // is dropped when nothing is left.
          int first_idx = -1;
          for (size_t k = 0; k < shift_vec.size(); k++) {
            if (shift_vec[k]->getLeftBpPos() == 0) {
              first_idx = static_cast<int>(k);
              break;
            }
          }
          MassShiftPtr first_shift_ptr =
              first_idx < 0 ? nullptr : shift_vec[first_idx];
          GraphNTermForm n_form = chooseNTermForm(form_ptr, first_shift_ptr);
          if (n_form.takes_first_shift) {
            if (n_form.drop_first_shift) {
              shift_vec.erase(shift_vec.begin() + first_idx);
            } else {
              MassShiftPtr remain_ptr =
                  std::make_shared<MassShift>(shift_vec[first_idx], 0);
              remain_ptr->setMassShift(n_form.remain_mass);
              shift_vec[first_idx] = remain_ptr;
            }
          }

          std::vector<PtmPtrVec> ptm_combs = choosePtmCombs(shift_vec);
          for (size_t k = 0; k < shift_vec.size(); k++) {
            const PtmPtrVec& ptm_vec = ptm_combs[k];

            MassShiftPtr shift_ptr =
                std::make_shared<toppic::MassShift>(shift_vec[k], 0);

            std::vector<double> mass_vec =
                mass_split(shift_vec[k]->getMassShift(), ptm_vec);

            AlterPtr alter_ptr = shift_vec[k]->getAlterPtr(0);

            AminoAcidPtr acid_ptr =
                alter_ptr->getModPtr()->getModResiduePtr()->getAminoAcidPtr();

            AlterPtrVec alter_vec;
            for (size_t i = 0; i < ptm_vec.size(); i++) {
              ResiduePtr mod_res =
                  std::make_shared<Residue>(acid_ptr, ptm_vec[i]);
              ModPtr mod = std::make_shared<Mod>(
                  alter_ptr->getModPtr()->getOriResiduePtr(), mod_res,
                  ModType::SIDE_CHAIN);
              AlterPtr a = std::make_shared<Alter>(
                  alter_ptr->getLeftBpPos(), alter_ptr->getRightBpPos(),
                  alter_ptr->getTypePtr(), mass_vec[i], mod);
              alter_vec.push_back(a);
            }

            shift_ptr->setAlterPtrVec(alter_vec);
            shift_vec[k] = shift_ptr;
          }

          // it used to be NME whenever the proteoform started at residue 1
          ProtModPtr prot_mod = n_form.prot_mod_ptr;

          MassShiftPtrVec unknown_shift_vec =
              prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(
                  AlterType::UNEXPECTED);

          for (size_t k = 0; k < unknown_shift_vec.size(); k++) {
            if (unknown_shift_vec[k]->getRightBpPos() ==
                unknown_shift_vec[k]->getLeftBpPos()) {
              int right_pos = unknown_shift_vec[k]->getRightBpPos();
              int left_pos = unknown_shift_vec[k]->getLeftBpPos();
              if (pair_vec[pair_vec.size() - 1]
                      ->getTheoPeakPtr()
                      ->getIonPtr()
                      ->getPos() > right_pos) {
                for (size_t i = 0; i < pair_vec.size(); i++) {
                  if (pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos() >
                      right_pos) {
                    right_pos =
                        pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos();
                    break;
                  }
                }
                unknown_shift_vec[k]->setRightBpPos(right_pos);
              } else {
                for (int i = static_cast<int>(pair_vec.size() - 1); i >= 0;
                     i--) {
                  if (pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos() <
                      left_pos) {
                    left_pos =
                        pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos();
                    break;
                  }
                }
                unknown_shift_vec[k]->setLeftBpPos(left_pos);
              }
            }
          }

          shift_vec.insert(shift_vec.end(), unknown_shift_vec.begin(),
                           unknown_shift_vec.end());

          std::sort(shift_vec.begin(), shift_vec.end(), MassShift::cmpPosInc);

          MassShiftPtrVec fix_shift_vec =
              prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(
                  AlterType::FIXED);

          fix_shift_vec.insert(fix_shift_vec.end(), shift_vec.begin(),
                               shift_vec.end());

          // The form as toppic builds it (proteoform_factory::
          // geneProtModProteoform): the residues from the form's start, the
          // acetylation as the N-terminal modification of the residue
          // sequence (its mass) and a PROTEIN_VARIABLE shift on the first
          // residue (the annotation).
          int offset = n_form.start_offset;
          ResSeqPtr res_seq_ptr = form_ptr->getResSeqPtr();
          ModPtr n_term_mod_ptr = prot_mod->getModPtr();
          bool acetylated = !ModBase::isNTermNoneModPtr(n_term_mod_ptr);
          if (offset > 0) {
            for (size_t k = 0; k < fix_shift_vec.size(); k++) {
              fix_shift_vec[k] = moveShift(fix_shift_vec[k], offset);
            }
          }
          if (offset > 0 || acetylated) {
            const ResiduePtrVec& residues = res_seq_ptr->getResidues();
            ResiduePtrVec form_residues(residues.begin() + offset,
                                        residues.end());
            res_seq_ptr =
                std::make_shared<ResidueSeq>(form_residues, n_term_mod_ptr);
          }
          if (acetylated) {
            AlterPtr alter_ptr = std::make_shared<Alter>(
                0, 1, AlterType::PROTEIN_VARIABLE, n_term_mod_ptr->getShift(),
                n_term_mod_ptr);
            fix_shift_vec.push_back(std::make_shared<MassShift>(alter_ptr));
          }

          ProteoformPtr new_form = std::make_shared<Proteoform>(
              form_ptr->getFastaSeqPtr(), prot_mod,
              form_ptr->getStartPos() + offset, form_ptr->getEndPos(),
              res_seq_ptr, fix_shift_vec);
          // a truncation N-terminal modification: the residue sequence's
          // N-terminal modification, as toppic's variable PTM search sets it
          if (n_form.trunc_n_term_mod_ptr != nullptr &&
              !new_form->setTruncNTermModPtr(n_form.trunc_n_term_mod_ptr)) {
            LOG_WARN("The truncation N-terminal modification of a proteoform "
                     "of " << new_form->getSeqName() << " cannot be added!");
          }

          // new_form->setVariablePtmNum(prsm_ptr->getProteoformPtr()->getVariablePtmNum());

          PrsmPtr new_prsm = std::make_shared<Prsm>(
              new_form, spec_set_ptr->getDeconvMsPtrVec(), adjusted_prec_mass,
              prsm_para_ptr->getSpParaPtr());
          // a PrSM rebuilt from a graph alignment is still TopMG's
          new_prsm->setTool(prsm_ptr->getTool());

          if (new_prsm->getMatchTheoFragNum() > 0) prsm_writer->write(new_prsm);

          prsm_ptr = prsm_reader->readOnePrsm(
              fasta_reader, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
        }
      }
    }
    std::cout << std::flush << "Mass graph - post-processing " << cnt << " of "
              << spectrum_num << " spectra.\r";
    deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  }
  prsm_reader->close();
  prsm_writer->close();
  std::cout << std::endl;
}

}  // namespace toppic
