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

#include "stat/local_align/local_align_processor.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "common/base/mod_util.hpp"
#include "common/base/ptm_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "seq/alter_type.hpp"
#include "seq/local_anno.hpp"

namespace toppic {

LocalAlignProcessor::LocalAlignProcessor(const LocalMngPtr& mng_ptr)
    : mng_ptr_(mng_ptr) {
  init();
}

void LocalAlignProcessor::init() {
  para_.ptm_ptr_vec = ptm_util::readPtmTxt(mng_ptr_->residueModFileName_);
  std::vector<ModPtrVec> mod_ptr_vec2d =
      mod_util::readModTxt(mng_ptr_->residueModFileName_);
  para_.mod_list_N = mod_ptr_vec2d[0];
  para_.mod_list_C = mod_ptr_vec2d[1];
  para_.mod_list_any = mod_ptr_vec2d[2];
  readPara();
  const char* log_name = std::getenv("TOPPIC_LOCAL_ALIGN_LOG");
  if (log_name != nullptr && log_name[0] != '\0') {
    log_file_.open(log_name);
    log_file_ << "spectrum_id\tprsm_id\tprotein\tshift_idx\tshift_mass\tpiece"
              << "\texplanations\tdecision\tsites" << std::endl;
  }
}

void LocalAlignProcessor::readPara() {
  const char* str = std::getenv("TOPPIC_LOCAL_ALIGN_PARA");
  if (str == nullptr) return;
  std::vector<std::string> items = str_util::split(str, ",");
  for (size_t i = 0; i < items.size(); i++) {
    std::vector<std::string> kv = str_util::split(items[i], "=");
    if (kv.size() != 2) continue;
    std::string name = kv[0];
    str_util::trim(name);
    double value = std::stod(kv[1]);
    if (name == "iso") {
      para_.iso = value != 0;
    } else if (name == "max_ptm") {
      para_.max_ptm = static_cast<int>(value);
    } else if (name == "split") {
      para_.split = value != 0;
    } else {
      LOG_ERROR("Unknown TOPPIC_LOCAL_ALIGN_PARA switch " << name);
    }
  }
}

void LocalAlignProcessor::log(const PrsmPtr& prsm_ptr, int shift_idx,
                              double shift_mass, const std::string& piece_name,
                              const std::string& expl_names,
                              const std::string& decision,
                              const std::string& sites) {
  if (!log_file_.is_open()) return;
  log_file_ << prsm_ptr->getSpectrumId() << "\t" << prsm_ptr->getPrsmId()
            << "\t" << prsm_ptr->getProteoformPtr()->getSeqName() << "\t"
            << shift_idx << "\t" << std::fixed << std::setprecision(4)
            << shift_mass << "\t" << piece_name << "\t" << expl_names << "\t"
            << decision << "\t" << sites << std::endl;
}

// LocalProcessor::process with this class's PrSM handling
void LocalAlignProcessor::process() {
  std::string spec_file_name = mng_ptr_->prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name =
      file_util::basename(spec_file_name) + "." + mng_ptr_->input_file_ext_;
  std::string output_file_name =
      file_util::basename(spec_file_name) + "." + mng_ptr_->output_file_ext_;
  PrsmXmlWriterPtr prsm_writer =
      std::make_shared<PrsmXmlWriter>(output_file_name);
  std::string db_file_name =
      mng_ptr_->prsm_para_ptr_->getSearchDbFileNameWithFolder();

  FastaIndexReaderPtr seq_reader =
      std::make_shared<FastaIndexReader>(db_file_name);
  PrsmReaderPtr prsm_reader = std::make_shared<PrsmReader>(input_file_name);
  ModPtrVec fix_mod_list = mng_ptr_->prsm_para_ptr_->getFixModPtrVec();
  PrsmPtr prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_list);
  int group_spec_num = mng_ptr_->prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();

  MsAlignReaderPtr reader_ptr = std::make_shared<MsAlignReader>(
      spec_file_name, group_spec_num, sp_para_ptr->getActivationPtr());

  int spectrum_num = msalign_util::getSpNum(spec_file_name);

  int cnt = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() != 0) {
    cnt += group_spec_num;
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
    if (header_ptr->containsPrec()) {
      double prec_mono_mass = header_ptr->getFirstPrecMonoMass();
      SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
          deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        int spec_id = spec_set_ptr->getSpectrumId();
        while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
          DeconvMsPtrVec set_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
          prsm_ptr->setDeconvMsPtrVec(set_ms_ptr_vec);
          double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
          ExtendMsPtrVec extend_ms_ptr_vec =
              extend_ms_factory::geneMsThreePtrVec(set_ms_ptr_vec, sp_para_ptr,
                                                   new_prec_mass);
          prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
          if (prsm_ptr->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED) >
              0) {
            prsm_ptr = processOnePrsm(prsm_ptr);
          }
          prsm_writer->write(prsm_ptr);
          prsm_ptr = prsm_reader->readOnePrsm(seq_reader, fix_mod_list);
        }
      }
      std::cout << std::flush << "PTM characterization is processing " << cnt
                << " of " << spectrum_num << " spectra.\r";
    }
    deconv_ms_ptr_vec = reader_ptr->getNextMsPtrVec();
  }
  prsm_reader->close();
  prsm_writer->close();
  std::cout << std::endl;
}

namespace {

std::string shiftSites(const MassShiftPtrVec& shifts, int start_pos) {
  std::ostringstream out;
  for (size_t k = 0; k < shifts.size(); k++) {
    LocalAnnoPtr anno = shifts[k]->getAlterPtr(0)->getLocalAnno();
    if (k > 0) out << ";";
    // the residue numbers of the protein, 1-based
    out << anno->getPtmPtr()->getAbbrName() << "@"
        << (anno->getLeftBpPos() + start_pos + 1) << "-"
        << (anno->getRightBpPos() + start_pos + 1) << ":" << std::fixed
        << std::setprecision(3) << anno->getConf();
  }
  return out.str();
}

}  // namespace

bool LocalAlignProcessor::processPiece(const PrsmPtr& prsm_ptr,
                                       const local_align::Piece& piece,
                                       const std::string& piece_name,
                                       int shift_idx,
                                       MassShiftPtrVec& new_shifts) {
  double err_tole = mng_ptr_->prec_tole_ptr_->compStrictErrorTole(
      prsm_ptr->getAdjustedPrecMass());
  std::vector<local_align::Explanation> expls =
      local_align::findExplanations(piece, err_tole, para_);
  std::string names;
  for (size_t k = 0; k < expls.size(); k++) {
    if (k > 0) names += ";";
    names += expls[k].getName();
  }
  if (expls.empty()) {
    log(prsm_ptr, shift_idx, piece.shift_mass, piece_name, "-", "none", "-");
    return false;
  }
  if (expls.size() > 1) {
    log(prsm_ptr, shift_idx, piece.shift_mass, piece_name, names, "ambiguous",
        "-");
    return false;
  }
  int best_score = 0;
  MassShiftPtrVec shifts = local_align::localize(
      piece, expls[0], prsm_ptr->getRefineMsPtrVec(), para_, mng_ptr_,
      best_score);
  if (shifts.empty()) {
    log(prsm_ptr, shift_idx, piece.shift_mass, piece_name, names, "no_sites",
        "-");
    return false;
  }
  log(prsm_ptr, shift_idx, piece.shift_mass, piece_name, names, "localized",
      shiftSites(shifts, prsm_ptr->getProteoformPtr()->getStartPos()) +
          " score=" + std::to_string(best_score));
  new_shifts.insert(new_shifts.end(), shifts.begin(), shifts.end());
  return true;
}

PrsmPtr LocalAlignProcessor::processOnePrsm(const PrsmPtr& prsm_ptr) {
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  MassShiftPtrVec unexp_shifts =
      form_ptr->getMassShiftPtrVec(AlterType::UNEXPECTED);
  std::sort(unexp_shifts.begin(), unexp_shifts.end(), MassShift::cmpPosInc);
  int len = form_ptr->getLen();
  double err_tole = mng_ptr_->prec_tole_ptr_->compStrictErrorTole(
      prsm_ptr->getAdjustedPrecMass());
  MassShiftPtrVec explained;
  MassShiftPtrVec new_shifts;
  if (unexp_shifts.size() == 1) {
    double mass = unexp_shifts[0]->getMassShift();
    // a shift within [-1, 1] Da is not characterized (LocalProcessor)
    if (std::abs(mass) <= 1 + err_tole) {
      log(prsm_ptr, 1, mass, "whole", "-", "small", "-");
      return prsm_ptr;
    }
    local_align::Piece piece =
        local_align::makePiece(form_ptr, 0, len, unexp_shifts[0]);
    if (processPiece(prsm_ptr, piece, "whole", 1, new_shifts)) {
      explained.push_back(unexp_shifts[0]);
    }
  } else if (unexp_shifts.size() == 2) {
    if (!para_.split) return prsm_ptr;
    if (unexp_shifts[0]->getRightBpPos() > unexp_shifts[1]->getLeftBpPos()) {
      log(prsm_ptr, 0, 0, "whole", "-", "overlap", "-");
      return prsm_ptr;
    }
    std::vector<int> matched = local_align::compMatchedSites(
        form_ptr, prsm_ptr->getRefineMsPtrVec(), mng_ptr_);
    int matched_num = 0;
    int break_site = local_align::findBreakSite(unexp_shifts[0], unexp_shifts[1],
                                                matched, matched_num);
    std::string break_str = "break=" + std::to_string(break_site) +
                            "/matched=" + std::to_string(matched_num);
    local_align::Piece piece_1 =
        local_align::makePiece(form_ptr, 0, break_site, unexp_shifts[0]);
    if (processPiece(prsm_ptr, piece_1, "1:0-" + std::to_string(break_site) +
                     "," + break_str, 1, new_shifts)) {
      explained.push_back(unexp_shifts[0]);
    }
    local_align::Piece piece_2 =
        local_align::makePiece(form_ptr, break_site, len, unexp_shifts[1]);
    if (processPiece(prsm_ptr, piece_2, "2:" + std::to_string(break_site) + "-" +
                     std::to_string(len) + "," + break_str, 2, new_shifts)) {
      explained.push_back(unexp_shifts[1]);
    }
  } else {
    return prsm_ptr;
  }
  if (explained.empty()) return prsm_ptr;

  MassShiftPtrVec all_shifts;
  const MassShiftPtrVec& ori_shifts = form_ptr->getMassShiftPtrVec();
  for (size_t k = 0; k < ori_shifts.size(); k++) {
    if (std::find(explained.begin(), explained.end(), ori_shifts[k]) ==
        explained.end()) {
      all_shifts.push_back(ori_shifts[k]);
    }
  }
  all_shifts.insert(all_shifts.end(), new_shifts.begin(), new_shifts.end());
  std::sort(all_shifts.begin(), all_shifts.end(), MassShift::cmpPosInc);
  ProteoformPtr new_form_ptr = std::make_shared<Proteoform>(
      form_ptr->getFastaSeqPtr(), form_ptr->getProtModPtr(),
      form_ptr->getStartPos(), form_ptr->getEndPos(), form_ptr->getResSeqPtr(),
      all_shifts);
  if (form_ptr->hasTruncNTermMod()) {
    new_form_ptr->setTruncNTermModPtr(form_ptr->getTruncNTermModPtr());
  }
  new_form_ptr->setProteoClusterId(form_ptr->getProteoClusterId());
  new_form_ptr->setProtId(form_ptr->getProtId());
  new_form_ptr->setProtClusterId(form_ptr->getProtClusterId());
  prsm_ptr->setProteoformPtr(new_form_ptr,
                             mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  prsm_ptr->setAdjustedPrecMass(new_form_ptr->getMass());
  return prsm_ptr;
}

}  // namespace toppic
