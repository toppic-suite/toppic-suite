// Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

#include "base/ptm_util.hpp"
#include "base/prot_mod_base.hpp"
#include "base/file_util.hpp"
#include "base/mass_shift.hpp"

#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_util.hpp"
#include "spec/msalign_reader.hpp"

#include "prsm/prsm_reader.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_xml_writer.hpp"

#include "graphalign/graph_align_mng.hpp"
#include "graphalign/graph_post_processor.hpp"

namespace prot {

std::vector<double> mass_split(double mass, const std::vector<PtmPtr> & ptm_vec) {
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
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  std::string var_mod_file_name = mng_ptr_->var_mod_file_name_;

  std::vector<PtmPtr> ptm_vec = prot::ptm_util::readPtmTxt(var_mod_file_name);

  ProtModPtrVec prot_mod_vec = prsm_para_ptr->getProtModPtrVec();
  for (size_t i = 0; i < prot_mod_vec.size(); i++) {
    if (prot_mod_vec[i]->getType() == ProtModBase::getType_M_ACETYLATION()
        || prot_mod_vec[i]->getType() == ProtModBase::getType_NME_ACETYLATION()) {
      ptm_vec.push_back(PtmBase::getPtmPtrByAbbrName("Acetyl"));
      break;
    }
  }

  std::sort(ptm_vec.begin(), ptm_vec.end(),
            [](PtmPtr a, PtmPtr b) {return a->getName() < b->getName();});

  std::vector<PtmPtr>::iterator it
      = std::unique(ptm_vec.begin(), ptm_vec.end(),
                    [](PtmPtr a, PtmPtr b) {return a->getName() == b->getName();});
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
      int mass = std::ceil(ptm_vec[i]->getMonoMass() * mng_ptr_->convert_ratio_);
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
          mass_ptm_map_[new_mass].insert(mass_ptm_map_[new_mass].end(),
                                         mass_ptm_map_[cur_map_mass[cnt]].begin(),
                                         mass_ptm_map_[cur_map_mass[cnt]].end());
          mass_ptm_map_[new_mass].push_back(ptm_vec[i]);
        }
      }
    }
  }
}

void GraphPostProcessor::process() {
  PrsmParaPtr prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  std::string db_file_name = prsm_para_ptr->getSearchDbFileName();
  LOG_DEBUG("Search db file name " << db_file_name);
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();

  geneMassPtmMap();

  int spectrum_num = msalign_util::getSpNum(prsm_para_ptr->getSpectrumFileName());
  std::string input_file_name = file_util::basename(sp_file_name)+ "." + input_ext_;
  std::string output_file_name = file_util::basename(sp_file_name) + "." + output_ext_;

  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);
  PrsmXmlWriterPtr prsm_writer = std::make_shared<PrsmXmlWriter>(output_file_name);

  FastaIndexReaderPtr fasta_reader
      = std::make_shared<FastaIndexReader>(prsm_para_ptr->getSearchDbFileName());

  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(fasta_reader,
                                              mng_ptr_->prsm_para_ptr_->getFixModPtrVec());

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          sp_para_ptr->getActivationPtr(),
                          sp_para_ptr->getSkipList());

  int cnt = 0;

  SpectrumSetPtr spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];

  while (spec_set_ptr != nullptr) {
    cnt += group_spec_num;
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        double adjusted_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec refine_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, adjusted_prec_mass);
        PeakIonPairPtrVec pair_vec
            = peak_ion_pair_util::genePeakIonPairs(prsm_ptr->getProteoformPtr(),
                                                   refine_ms_ptr_vec[0], 
                                                   sp_para_ptr->getMinMass());
        std::sort(pair_vec.begin(), pair_vec.end(), PeakIonPair::cmpTheoPeakPosInc);

        MassShiftPtrVec shift_vec
            = prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(MassShiftType::VARIABLE);

        for (size_t k = 0; k < shift_vec.size(); k++) {
          int mass = std::ceil(shift_vec[k]->getMassShift() * mng_ptr_->convert_ratio_);
          PtmPtrVec ptm_vec;
          int err = std::abs(mass - mass_ptm_map_.begin()->first);
          for (auto it = mass_ptm_map_.begin(); it != mass_ptm_map_.end(); it++) {
            if (std::abs(mass - it->first) < err) {
              err = std::abs(mass - it->first);
              ptm_vec = it->second;
            }
          }

          MassShiftPtr shfit_ptr
              = std::make_shared<prot::MassShift>(shift_vec[k]->getLeftBpPos(),
                                                  shift_vec[k]->getRightBpPos(),
                                                  shift_vec[k]->getTypePtr());

          std::vector<double> mass_vec = mass_split(shift_vec[k]->getMassShift(), ptm_vec);

          ChangePtr change_ptr = shift_vec[k]->getChangePtr(0);

          AminoAcidPtr acid_ptr = change_ptr->getModPtr()->getModResiduePtr()->getAminoAcidPtr();

          for (size_t i = 0; i < ptm_vec.size(); i++) {
            ResiduePtr mod_res = std::make_shared<Residue>(acid_ptr, ptm_vec[i]);
            ModPtr mod = std::make_shared<Mod>(change_ptr->getModPtr()->getOriResiduePtr(), mod_res);
            ChangePtr c = std::make_shared<Change>(change_ptr->getLeftBpPos(),
                                                   change_ptr->getRightBpPos(),
                                                   change_ptr->getTypePtr(),
                                                   mass_vec[i],
                                                   mod);
            shfit_ptr->setChangePtr(c);
          }

          shift_vec[k] = shfit_ptr;
        }

        ProtModPtr prot_mod = prsm_ptr->getProteoformPtr()->getProtModPtr();

        if (prsm_ptr->getProteoformPtr()->getStartPos() == 1) {
          prot_mod = ProtModBase::getProtModPtrByName(ProtModBase::getType_NME());
        }

        MassShiftPtrVec unknown_shift_vec
            = prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(MassShiftType::UNEXPECTED);

        for (size_t k = 0; k < unknown_shift_vec.size(); k++) {
          if (unknown_shift_vec[k]->getRightBpPos() == unknown_shift_vec[k]->getLeftBpPos()) {
            int right_pos = unknown_shift_vec[k]->getRightBpPos();
            int left_pos = unknown_shift_vec[k]->getLeftBpPos();
            if (pair_vec[pair_vec.size() - 1]->getTheoPeakPtr()->getIonPtr()->getPos() > right_pos) {
              for (size_t i = 0; i < pair_vec.size(); i++) {
                if (pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos() > right_pos) {
                  right_pos = pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos();
                  break; 
                }
              }
              unknown_shift_vec[k]->setRightBpPos(right_pos);
            } else {
              for (int i = static_cast<int>(pair_vec.size() - 1); i >= 0; i--) {
                if (pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos() < left_pos) {
                  left_pos = pair_vec[i]->getTheoPeakPtr()->getIonPtr()->getPos();
                  break; 
                }
              }
              unknown_shift_vec[k]->setLeftBpPos(left_pos);
            }
          }
        }

        shift_vec.insert(shift_vec.end(), unknown_shift_vec.begin(), unknown_shift_vec.end());

        std::sort(shift_vec.begin(), shift_vec.end(), MassShift::cmpPosInc);

        MassShiftPtrVec fix_shift_vec
            = prsm_ptr->getProteoformPtr()->getMassShiftPtrVec(MassShiftType::FIXED);

        fix_shift_vec.insert(fix_shift_vec.end(), shift_vec.begin(), shift_vec.end());

        ProteoformPtr new_form
            = std::make_shared<Proteoform>(prsm_ptr->getProteoformPtr()->getFastaSeqPtr(),
                                           prot_mod,
                                           prsm_ptr->getProteoformPtr()->getStartPos(),
                                           prsm_ptr->getProteoformPtr()->getEndPos(),
                                           prsm_ptr->getProteoformPtr()->getResSeqPtr(),
                                           fix_shift_vec);

        new_form->setVariablePtmNum(prsm_ptr->getProteoformPtr()->getVariablePtmNum());

        PrsmPtr new_prsm = std::make_shared<Prsm>(new_form, spec_set_ptr->getDeconvMsPtrVec(),
                                                  adjusted_prec_mass,
                                                  prsm_para_ptr->getSpParaPtr());

        if (new_prsm->getMatchFragNum() > 0) prsm_writer->write(new_prsm);

        prsm_ptr = prsm_reader->readOnePrsm(fasta_reader,
                                            mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
      }
    }
    std::cout << std::flush <<  "Mass graph - post-processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
    spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr)[0];
  }
  sp_reader.close();
  prsm_reader->close();
  prsm_writer->close();
  std::cout << std::endl;
}

}  // namespace prot
