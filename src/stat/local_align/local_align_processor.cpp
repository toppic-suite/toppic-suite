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
#include "common/base/prot_mod_base.hpp"
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
#include "stat/local/local_proteoform.hpp"

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
              << "\texplanations\tdecision\tsites\trule\tall_explanations"
              << std::endl;
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
    } else if (name == "term") {
      para_.term = value != 0;
    } else if (name == "term_range") {
      para_.term_range = static_cast<int>(value);
    } else if (name == "term_max_matched") {
      para_.term_max_matched = static_cast<int>(value);
    } else if (name == "tie") {
      para_.tie = value != 0;
    } else {
      LOG_ERROR("Unknown TOPPIC_LOCAL_ALIGN_PARA switch " << name);
    }
  }
}

void LocalAlignProcessor::log(const PrsmPtr& prsm_ptr, int shift_idx,
                              double shift_mass, const std::string& piece_name,
                              const std::string& expl_names,
                              const std::string& decision,
                              const std::string& sites,
                              const std::string& rule,
                              const std::string& all_names) {
  if (!log_file_.is_open()) return;
  log_file_ << prsm_ptr->getSpectrumId() << "\t" << prsm_ptr->getPrsmId()
            << "\t" << prsm_ptr->getProteoformPtr()->getSeqName() << "\t"
            << shift_idx << "\t" << std::fixed << std::setprecision(4)
            << shift_mass << "\t" << piece_name << "\t" << expl_names << "\t"
            << decision << "\t" << sites << "\t" << (rule.empty() ? "-" : rule)
            << "\t" << (all_names.empty() ? "-" : all_names) << std::endl;
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

namespace {

bool isAcetylForm(const std::string& type) {
  return type == ProtModBase::getType_M_ACETYLATION() ||
         type == ProtModBase::getType_NME_ACETYLATION();
}

}  // namespace

int LocalAlignProcessor::choose(const PrsmPtr& prsm_ptr,
                                const std::vector<Cand>& cands,
                                MassShiftPtrVec& shifts, int& best_score,
                                std::string& rule) {
  rule = "";
  std::vector<size_t> kept;
  for (size_t k = 0; k < cands.size(); k++) kept.push_back(k);
  if (cands.size() > 1 && para_.tie) {
    std::vector<local_align::Explanation> pool;
    std::vector<bool> acetyl_forms;
    for (size_t k = 0; k < cands.size(); k++) {
      pool.push_back(cands[k].expl);
      acetyl_forms.push_back(cands[k].acetyl_form);
    }
    kept = local_align::breakTies(pool, acetyl_forms, rule);
  }
  if (kept.size() > 1 && para_.tie) {
    // Acetyl vs Trimethyl: the explanations left are on the same
    // proteoform and the same up to Acetyl <-> Trimethyl
    auto key = [&](size_t i) {
      std::vector<std::string> abbrs;
      abbrs.push_back(cands[i].expl.first->getAbbrName());
      if (cands[i].expl.second != nullptr) {
        abbrs.push_back(cands[i].expl.second->getAbbrName());
      }
      for (size_t k = 0; k < abbrs.size(); k++) {
        if (abbrs[k] == "Trimethyl") abbrs[k] = "Acetyl";
      }
      std::sort(abbrs.begin(), abbrs.end());
      return abbrs[0] + (abbrs.size() > 1 ? "+" + abbrs[1] : "");
    };
    bool eligible = true;
    bool any_acetyl = false;
    bool any_trimethyl = false;
    for (size_t k = 0; k < kept.size(); k++) {
      if (key(kept[k]) != key(kept[0]) ||
          cands[kept[k]].form != cands[kept[0]].form) {
        eligible = false;
      }
      if (cands[kept[k]].expl.hasPtm("Acetyl")) any_acetyl = true;
      if (cands[kept[k]].expl.hasPtm("Trimethyl")) any_trimethyl = true;
    }
    if (eligible && any_acetyl && any_trimethyl) {
      int best = -1;
      int best_scr = -1;
      double best_conf = -1;
      MassShiftPtrVec best_shifts;
      std::string detail;
      for (size_t k = 0; k < kept.size(); k++) {
        const Cand& cand = cands[kept[k]];
        int scr = 0;
        MassShiftPtrVec cur = local_align::localize(
            cand.piece, cand.expl, prsm_ptr->getRefineMsPtrVec(), para_,
            mng_ptr_, scr);
        if (cur.empty()) scr = -1;
        double conf = 0;
        for (size_t m = 0; m < cur.size(); m++) {
          conf += cur[m]->getAlterPtr(0)->getLocalAnno()->getConf();
        }
        detail += (detail.empty() ? "" : "/") + cand.expl.getName() + ":" +
                  std::to_string(scr);
        bool better = scr > best_scr || (scr == best_scr && conf > best_conf) ||
                      (scr == best_scr && conf == best_conf &&
                       cand.expl.hasPtm("Acetyl") && best >= 0 &&
                       !cands[best].expl.hasPtm("Acetyl"));
        if (best < 0 || better) {
          best = static_cast<int>(kept[k]);
          best_scr = scr;
          best_conf = conf;
          best_shifts = cur;
        }
      }
      rule += (rule.empty() ? "" : "+") + std::string("site_scores(") + detail + ")";
      shifts = best_shifts;
      best_score = best_scr;
      return best;
    }
  }
  if (kept.size() != 1) return -1;
  const Cand& cand = cands[kept[0]];
  shifts = local_align::localize(cand.piece, cand.expl,
                                 prsm_ptr->getRefineMsPtrVec(), para_,
                                 mng_ptr_, best_score);
  return static_cast<int>(kept[0]);
}

LocalAlignProcessor::PieceResult LocalAlignProcessor::processPiece(
    const PrsmPtr& prsm_ptr, const local_align::Piece& piece,
    const std::string& piece_name, int shift_idx,
    const MassShiftPtrVec& kept_shifts, bool n_forms) {
  PieceResult result;
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  double adj_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->prec_tole_ptr_->compStrictErrorTole(adj_prec_mass);
  // the candidates: the explanations of the piece on the proteoform as it
  // is (form nullptr) and, with n_forms, on its other N-terminal forms of
  // the same start (the acetylation dropped or added), as main's base
  // candidate forms
  std::vector<Cand> cands;
  std::string ori_type = form_ptr->getProtModPtr()->getType();
  std::vector<local_align::Explanation> expls =
      local_align::findExplanations(piece, err_tole, para_);
  for (size_t k = 0; k < expls.size(); k++) {
    cands.push_back({nullptr, piece, expls[k], isAcetylForm(ori_type),
                     expls[k].getName()});
  }
  if (n_forms) {
    MassShiftPtrVec kept = kept_shifts;
    ProteoformPtrVec forms = local_proteoform::getCandidateForm(
        form_ptr->getFastaSeqPtr(), form_ptr->getStartPos(),
        form_ptr->getStartPos(), form_ptr->getEndPos(), kept, mng_ptr_);
    for (size_t f = 0; f < forms.size(); f++) {
      std::string type = forms[f]->getProtModPtr()->getType();
      if (type == ori_type) continue;
      local_align::Piece alt =
          local_align::makePiece(forms[f], piece.bgn, piece.end, nullptr);
      alt.shift_mass = adj_prec_mass - forms[f]->getMass();
      std::vector<local_align::Explanation> alt_expls =
          local_align::findExplanations(alt, err_tole, para_);
      for (size_t k = 0; k < alt_expls.size(); k++) {
        cands.push_back({forms[f], alt, alt_expls[k], isAcetylForm(type),
                         alt_expls[k].getName() + "/" + type});
      }
    }
  }
  std::string names;
  for (size_t k = 0; k < cands.size(); k++) {
    if (k > 0) names += ";";
    names += cands[k].name;
  }
  if (cands.empty()) {
    result.decision = "none";
    log(prsm_ptr, shift_idx, piece.shift_mass, piece_name, "-", "none", "-");
    return result;
  }
  // the tie-break rules and the site scores, then the one-explanation rule
  std::string rule;
  MassShiftPtrVec shifts;
  int best_score = 0;
  int chosen = choose(prsm_ptr, cands, shifts, best_score, rule);
  if (chosen < 0) {
    result.decision = "ambiguous";
    log(prsm_ptr, shift_idx, piece.shift_mass, piece_name, names, "ambiguous",
        "-", rule);
    return result;
  }
  const Cand& cand = cands[chosen];
  std::string all_names = cands.size() > 1 ? names : "";
  if (cands.size() > 1) names = cand.name;
  result.shifts = shifts;
  if (result.shifts.empty()) {
    result.decision = "no_sites";
    log(prsm_ptr, shift_idx, cand.piece.shift_mass, piece_name, names,
        "no_sites", "-", rule, all_names);
    return result;
  }
  result.ok = true;
  result.cand_form = cand.form;
  result.decision = "localized";
  log(prsm_ptr, shift_idx, cand.piece.shift_mass, piece_name, names,
      "localized",
      shiftSites(result.shifts, form_ptr->getStartPos()) +
          " score=" + std::to_string(best_score),
      rule, all_names);
  return result;
}

LocalAlignProcessor::PieceResult LocalAlignProcessor::processTerminal(
    const PrsmPtr& prsm_ptr, const MassShiftPtr& shift_ptr, bool n_side,
    int piece_bgn, int piece_end, const MassShiftPtrVec& kept_shifts,
    const std::string& piece_name, int shift_idx) {
  PieceResult result;
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  int start = form_ptr->getStartPos();
  int end = form_ptr->getEndPos();
  int len = form_ptr->getLen();
  int prot_len = static_cast<int>(form_ptr->getFastaSeqPtr()->getRawSeq().length());
  double adj_prec_mass = prsm_ptr->getAdjustedPrecMass();
  double err_tole = mng_ptr_->prec_tole_ptr_->compStrictErrorTole(adj_prec_mass);
  // the candidates on the moved termini
  std::vector<Cand> cands;
  std::string names;
  for (int d = -para_.term_range; d <= para_.term_range; d++) {
    // the unmoved terminus, with its N-terminal forms, is the standard
    // attempt (processPiece)
    if (d == 0) continue;
    int new_start = start;
    int new_end = end;
    // The piece keeps at least one residue. A truncation may cut into the
    // shift's range: the shift's residues near the terminus are the ones a
    // truncation explains (a -128.07 Da shift on the residues 1-9 is a
    // Dimethyl once the first residue, an R, is dropped), and the PTM
    // sites are searched over the whole piece of the candidate.
    if (n_side) {
      new_start = start + d;
      if (new_start < 0 || d >= piece_end) continue;
    } else {
      new_end = end + d;
      if (new_end > prot_len - 1 || len + d <= piece_bgn) continue;
    }
    MassShiftPtrVec kept = kept_shifts;
    ProteoformPtrVec forms = local_proteoform::getCandidateForm(
        form_ptr->getFastaSeqPtr(), start, new_start, new_end, kept, mng_ptr_);
    for (size_t f = 0; f < forms.size(); f++) {
      int cand_len = forms[f]->getLen();
      int pb = n_side ? 0 : piece_bgn;
      int pe = n_side ? piece_end - d : cand_len;
      if (pe <= pb) continue;
      // the other shifts of the candidate are all known
      local_align::Piece piece =
          local_align::makePiece(forms[f], pb, pe, nullptr);
      piece.shift_mass = adj_prec_mass - forms[f]->getMass();
      std::vector<local_align::Explanation> expls =
          local_align::findExplanations(piece, err_tole, para_);
      for (size_t k = 0; k < expls.size(); k++) {
        std::string type = forms[f]->getProtModPtr()->getType();
        std::string name = expls[k].getName() + (n_side ? "@start" : "@end") +
                           std::to_string((n_side ? new_start : new_end) + 1) +
                           (type == "NONE" ? "" : "/" + type);
        if (!names.empty()) names += ";";
        names += name;
        cands.push_back({forms[f], piece, expls[k], isAcetylForm(type), name});
      }
    }
  }
  std::string side = n_side ? "N" : "C";
  if (cands.empty()) {
    result.decision = "term_none";
    log(prsm_ptr, shift_idx, shift_ptr->getMassShift(), piece_name + "," + side,
        "-", "term_none", "-");
    return result;
  }
  // the tie-break rules and the site scores over all (candidate,
  // explanation) pairs, then the one-explanation rule
  std::string rule;
  MassShiftPtrVec shifts;
  int best_score = 0;
  int chosen = choose(prsm_ptr, cands, shifts, best_score, rule);
  if (chosen < 0) {
    result.decision = "term_ambiguous";
    log(prsm_ptr, shift_idx, shift_ptr->getMassShift(), piece_name + "," + side,
        names, "term_ambiguous", "-", rule);
    return result;
  }
  const Cand& cand = cands[chosen];
  std::string all_names = cands.size() > 1 ? names : "";
  if (cands.size() > 1) names = cand.name;
  result.shifts = shifts;
  if (result.shifts.empty()) {
    result.decision = "term_no_sites";
    log(prsm_ptr, shift_idx, shift_ptr->getMassShift(), piece_name + "," + side,
        names, "term_no_sites", "-", rule, all_names);
    return result;
  }
  result.ok = true;
  result.cand_form = cand.form;
  result.decision = "term_localized";
  log(prsm_ptr, shift_idx, cand.piece.shift_mass, piece_name + "," + side,
      names, "term_localized",
      shiftSites(result.shifts, cand.form->getStartPos()) +
          " score=" + std::to_string(best_score),
      rule, all_names);
  return result;
}

PrsmPtr LocalAlignProcessor::processOnePrsm(const PrsmPtr& prsm_ptr) {
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  MassShiftPtrVec unexp_shifts =
      form_ptr->getMassShiftPtrVec(AlterType::UNEXPECTED);
  std::sort(unexp_shifts.begin(), unexp_shifts.end(), MassShift::cmpPosInc);
  int len = form_ptr->getLen();
  int start = form_ptr->getStartPos();
  int end = form_ptr->getEndPos();
  double err_tole = mng_ptr_->prec_tole_ptr_->compStrictErrorTole(
      prsm_ptr->getAdjustedPrecMass());
  if (unexp_shifts.size() != 1 && unexp_shifts.size() != 2) return prsm_ptr;
  if (unexp_shifts.size() == 2 && !para_.split) return prsm_ptr;

  // the shifts of the proteoform other than the explained one(s) and the
  // protein N-terminal ones, which the candidate forms add themselves
  auto keptShifts = [&](const MassShiftPtrVec& explained) {
    MassShiftPtrVec kept;
    const MassShiftPtrVec& all = form_ptr->getMassShiftPtrVec();
    for (size_t k = 0; k < all.size(); k++) {
      if (all[k]->getTypePtr() == AlterType::PROTEIN_VARIABLE) continue;
      if (std::find(explained.begin(), explained.end(), all[k]) !=
          explained.end()) {
        continue;
      }
      kept.push_back(all[k]);
    }
    return kept;
  };
  std::vector<int> matched = local_align::compMatchedSites(
      form_ptr, prsm_ptr->getRefineMsPtrVec(), mng_ptr_);
  auto fewOnN = [&](const MassShiftPtr& s) {
    return para_.term &&
           local_align::countMatchedSites(matched, 1, s->getLeftBpPos()) <=
               para_.term_max_matched;
  };
  auto fewOnC = [&](const MassShiftPtr& s) {
    return para_.term &&
           local_align::countMatchedSites(matched, s->getRightBpPos(), len - 1) <=
               para_.term_max_matched;
  };

  MassShiftPtrVec explained;
  std::vector<PieceResult> results;
  if (unexp_shifts.size() == 1) {
    MassShiftPtr s = unexp_shifts[0];
    double mass = s->getMassShift();
    // a shift within [-1, 1] Da is not characterized (LocalProcessor)
    if (std::abs(mass) <= 1 + err_tole) {
      log(prsm_ptr, 1, mass, "whole", "-", "small", "-");
      return prsm_ptr;
    }
    local_align::Piece piece = local_align::makePiece(form_ptr, 0, len, s);
    MassShiftPtrVec kept = keptShifts({s});
    PieceResult r = processPiece(prsm_ptr, piece, "whole", 1, kept, fewOnN(s));
    if (!r.ok && r.decision == "none") {
      if (fewOnN(s)) {
        r = processTerminal(prsm_ptr, s, true, 0, len, kept, "whole", 1);
      }
      if (!r.ok && fewOnC(s)) {
        r = processTerminal(prsm_ptr, s, false, 0, len, kept, "whole", 1);
      }
    }
    if (r.ok) {
      explained.push_back(s);
      results.push_back(r);
    }
  } else {
    MassShiftPtr s1 = unexp_shifts[0];
    MassShiftPtr s2 = unexp_shifts[1];
    if (s1->getRightBpPos() > s2->getLeftBpPos()) {
      log(prsm_ptr, 0, 0, "whole", "-", "overlap", "-");
      return prsm_ptr;
    }
    int matched_num = 0;
    int b = local_align::findBreakSite(s1, s2, matched, matched_num);
    std::string break_str =
        "break=" + std::to_string(b) + "/matched=" + std::to_string(matched_num);
    std::string name_1 = "1:0-" + std::to_string(b) + "," + break_str;
    std::string name_2 = "2:" + std::to_string(b) + "-" + std::to_string(len) +
                         "," + break_str;
    local_align::Piece piece_1 = local_align::makePiece(form_ptr, 0, b, s1);
    PieceResult r1 = processPiece(prsm_ptr, piece_1, name_1, 1, keptShifts({s1}),
                                  fewOnN(s1));
    if (!r1.ok && r1.decision == "none" && fewOnN(s1)) {
      r1 = processTerminal(prsm_ptr, s1, true, 0, b, keptShifts({s1}), name_1,
                           1);
    }
    local_align::Piece piece_2 = local_align::makePiece(form_ptr, b, len, s2);
    PieceResult r2 = processPiece(prsm_ptr, piece_2, name_2, 2, keptShifts({s2}),
                                  false);
    if (!r2.ok && r2.decision == "none" && fewOnC(s2)) {
      r2 = processTerminal(prsm_ptr, s2, false, b, len, keptShifts({s2}),
                           name_2, 2);
    }
    if (r1.ok) {
      explained.push_back(s1);
      results.push_back(r1);
    }
    if (r2.ok) {
      explained.push_back(s2);
      results.push_back(r2);
    }
  }
  if (explained.empty()) return prsm_ptr;

  // the result proteoform: the original one, or the one with the moved
  // start / end (or the other N-terminal form) of the terminal candidates
  int final_start = start;
  int final_end = end;
  ProtModPtr n_prot_mod = form_ptr->getProtModPtr();
  bool n_changed = false;
  for (size_t k = 0; k < results.size(); k++) {
    if (results[k].cand_form == nullptr) continue;
    if (results[k].cand_form->getStartPos() != start ||
        results[k].cand_form->getProtModPtr()->getType() !=
            n_prot_mod->getType()) {
      final_start = results[k].cand_form->getStartPos();
      n_prot_mod = results[k].cand_form->getProtModPtr();
      n_changed = true;
    }
    if (results[k].cand_form->getEndPos() != end) {
      final_end = results[k].cand_form->getEndPos();
    }
  }
  ProteoformPtr base_ptr;
  MassShiftPtrVec all_shifts;
  if (!n_changed && final_end == end) {
    base_ptr = form_ptr;
    const MassShiftPtrVec& ori_shifts = form_ptr->getMassShiftPtrVec();
    for (size_t k = 0; k < ori_shifts.size(); k++) {
      if (std::find(explained.begin(), explained.end(), ori_shifts[k]) ==
          explained.end()) {
        all_shifts.push_back(ori_shifts[k]);
      }
    }
  } else {
    MassShiftPtrVec kept = keptShifts(explained);
    ProteoformPtrVec forms = local_proteoform::getCandidateForm(
        form_ptr->getFastaSeqPtr(), start, final_start, final_end, kept,
        mng_ptr_);
    base_ptr = forms[0];
    for (size_t f = 0; f < forms.size(); f++) {
      if (forms[f]->getProtModPtr()->getType() == n_prot_mod->getType()) {
        base_ptr = forms[f];
      }
    }
    all_shifts = base_ptr->getMassShiftPtrVec();
  }
  for (size_t k = 0; k < results.size(); k++) {
    int cand_start = results[k].cand_form == nullptr
                         ? start
                         : results[k].cand_form->getStartPos();
    for (size_t m = 0; m < results[k].shifts.size(); m++) {
      all_shifts.push_back(local_align::shiftPositions(
          results[k].shifts[m], final_start - cand_start));
    }
  }
  std::sort(all_shifts.begin(), all_shifts.end(), MassShift::cmpPosInc);
  ProteoformPtr new_form_ptr = std::make_shared<Proteoform>(
      form_ptr->getFastaSeqPtr(), base_ptr->getProtModPtr(), final_start,
      final_end, base_ptr->getResSeqPtr(), all_shifts);
  if (form_ptr->hasTruncNTermMod() && final_start == start) {
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
