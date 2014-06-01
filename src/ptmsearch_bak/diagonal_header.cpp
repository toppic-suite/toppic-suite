/*
 * diagonal_header.cpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#include <cmath>

#include "base/prot_mod.hpp"
#include "base/change.hpp"
#include "base/algorithm.hpp"
#include "ptmsearch/diagonal_header.hpp"

namespace prot {

DiagonalHeader::DiagonalHeader(double n_term_shift, bool n_strict,
                               bool c_strict, bool n_trunc, bool c_trunc) {
  prot_N_term_shift_ = n_term_shift;
  n_strict_ = n_strict;
  c_strict_ = c_strict;
  n_trunc_ = n_trunc;
  c_trunc_ = c_trunc;
}
DiagonalHeaderPtr DiagonalHeader::clone() {
  DiagonalHeaderPtr cloned = DiagonalHeaderPtr(
      new DiagonalHeader(prot_N_term_shift_, n_strict_, c_strict_, 
                         n_trunc_, c_trunc_));
  cloned->setId(id_);
  cloned->setTruncFirstResPos(trunc_first_res_pos_);
  cloned->setMatchFirstBpPos(match_first_bp_pos_);
  cloned->setPepNTermShift(pep_N_term_shift_);
  cloned->setPepNTermMatch(pep_N_term_match_);
  cloned->setProtNTermShift(prot_N_term_shift_);
  cloned->setProtNTermMatch(prot_N_term_match_);
  cloned->setAlignPrefix(is_align_prefix_);
  cloned->setTruncLastResPos(trunc_last_res_pos_);
  cloned->setMatchLastBpPos(match_last_bp_pos_);
  cloned->setPepCTermShift(pep_C_term_shift_);
  cloned->setPepCTermMatch(pep_C_term_match_);
  cloned->setProtCTermShift(prot_C_term_shift_);
  cloned->setProtCTermMatch(prot_C_term_match_);
  cloned->setAlignSuffix(is_align_suffix_);
  return cloned;
}

DiagonalHeaderPtr getDiagonalHeaderPtr(int bgn, int end, DiagonalHeaderPtr header) {
  DiagonalHeaderPtr new_header = header->clone();
  new_header->setMatchFirstBpPos(bgn);
  new_header->setMatchLastBpPos(end);
  return new_header;
}

DiagonalHeaderPtrVec getNTermShiftListCommon(std::vector<double> best_shifts) {
  DiagonalHeaderPtrVec headers;
  for (unsigned int i = 0; i < best_shifts.size(); i++) {
    headers.push_back(
        DiagonalHeaderPtr(
            new DiagonalHeader(best_shifts[i], true, false, false, false)));
  }
  return headers;
}

DiagonalHeaderPtr getNTermShiftListCompLeft(ProteoformPtr seq,
                                               PtmMngPtr mng) {
  double shift = 0;
  return DiagonalHeaderPtr(
      new DiagonalHeader(shift, true, false, true, false));
}

DiagonalHeaderPtr getNTermShiftListCompRight(ProteoformPtr seq,
                                             PrmMsPtr ms_six) {
  double prec_mass = ms_six->getHeaderPtr()->getPrecMonoMass();
  double mole_mass = seq->getResSeqPtr()->getSeqMass();
  double shift = prec_mass - mole_mass;
  return DiagonalHeaderPtr(
      new DiagonalHeader(shift, false, true, false, true));
}

void setPrefixSuffix(DiagonalHeaderPtr &header, double c_shift,
                     ProteoformPtr seq, double term_error_tolerance,
                     PtmMngPtr mng) {
  // set protein c term shift
  header->setProtCTermShift(c_shift);

  std::vector<double> seq_b_masses = seq->getBpSpecPtr()->getPrmMasses();
  double prot_n_term_shift = header->getProtNTermShift();
  double prot_c_term_shift = header->getProtCTermShift();
  int trunc_first_res_pos = getFirstResPos(prot_n_term_shift,
                                           seq_b_masses);

  header->setTruncFirstResPos(trunc_first_res_pos);
  int trunc_last_res_pos = getLastResPos(prot_c_term_shift, seq_b_masses);
  header->setTruncLastResPos(trunc_last_res_pos);
  double pep_n_term_shift = prot_n_term_shift
      + seq_b_masses[trunc_first_res_pos];
  header->setPepNTermShift(pep_n_term_shift);
  double pep_c_term_shift = prot_c_term_shift
      + seq_b_masses[seq_b_masses.size() - 1]
      - seq_b_masses[trunc_last_res_pos + 1];
  header->setPepCTermShift(pep_c_term_shift);
  
  header->setProtTermMatch(term_error_tolerance);

  header->setPepTermMatch(term_error_tolerance);

  header->setAlignPrefixSuffix(mng->align_prefix_suffix_shift_thresh_);
}

void DiagonalHeader::setProtTermMatch(double term_error_tolerance) {
  if (std::abs(prot_N_term_shift_) <= term_error_tolerance) {
    setProtNTermMatch(true);
  }
  else {
    setProtNTermMatch(false);
  }

  if (std::abs(prot_C_term_shift_) <= term_error_tolerance) {
    setProtCTermMatch(true);
  }
  else {
    setProtCTermMatch(false);
  }
}

void DiagonalHeader::setPepTermMatch(double term_error_tolerance) {
  if (std::abs(pep_N_term_shift_) <= term_error_tolerance) {
    setPepNTermMatch(true);
  }
  else {
    setPepNTermMatch(false);
  }

  if (std::abs(pep_C_term_shift_) <= term_error_tolerance) {
    setPepCTermMatch(true);
  }
  else {
    setPepCTermMatch(false);
  }
}

void DiagonalHeader::setAlignPrefixSuffix(double term_error_tolerance) {
  if (std::abs(prot_N_term_shift_) <= term_error_tolerance) {
    setAlignPrefix(true);
  }
  else {
    setAlignPrefix(false);
  }

  if (std::abs(prot_C_term_shift_) <= term_error_tolerance) {
    setAlignSuffix(true);
  }
  else {
    setAlignSuffix(false);
  }
}

DiagonalHeaderPtrVec getNTermShiftListTruncPrefix(ProteoformPtr seq) {
  DiagonalHeaderPtrVec extend_n_term_shift;
  std::vector<double> seq_masses = seq->getBpSpecPtr()->getPrmMasses();
  double shift;
  for (unsigned int i = 1; i < seq_masses.size(); i++) {
    shift = -seq_masses[i];
    extend_n_term_shift.push_back(
        DiagonalHeaderPtr(new DiagonalHeader(shift, true, false, true, false)));
  }

  return extend_n_term_shift;
}

DiagonalHeaderPtrVec getNTermShiftListTruncsuffix(PrmMsPtr ms,
                                                  ProteoformPtr seq) {
  DiagonalHeaderPtrVec extend_n_term_shift;
  std::vector<double> ms_masses = prot::getMassList(ms);
  std::vector<double> seq_masses = seq->getBpSpecPtr()->getPrmMasses();
  double shift;
  for (unsigned int i = 1; i < seq_masses.size(); i++) {
    shift = ms_masses[ms_masses.size() - 1] - seq_masses[i];
    extend_n_term_shift.push_back(
        DiagonalHeaderPtr(new DiagonalHeader(shift, true, false, true, false)));
  }
  return extend_n_term_shift;
}
DiagonalHeaderPtrVec get1dHeaders(DiagonalHeaderPtrVec2D headers) {
  DiagonalHeaderPtrVec header_list;
  for (unsigned int i = 0; i < headers.size(); i++) {
    for (unsigned int j = 0; j < headers[i].size(); j++) {
      header_list.push_back(headers[i][j]);
    }
  }
  return header_list;
}

ChangePtrVec getUnexpectedChanges(DiagonalHeaderPtrVec headers, int first, int last) {
  ChangePtrVec change_list;
  if (!headers[0]->isPepNTermMatch()) {
    change_list.push_back(
        ChangePtr(new Change(0, headers[0]->getMatchFirstBpPos()-first,
                             UNEXPECTED_CHANGE, headers[0]->getPepNTermShift(), nullptr)));
  }
  for (unsigned int i = 0; i < headers.size() - 1; i++) {
    change_list.push_back(
        ChangePtr(
            new Change(
                headers[i]->getMatchLastBpPos()-first,
                headers[i + 1]->getMatchFirstBpPos()-first,
                UNEXPECTED_CHANGE,
                headers[i + 1]->getProtNTermShift()
                    - headers[i]->getProtNTermShift(),
                nullptr)));
  }
  DiagonalHeaderPtr lastHeader = headers[headers.size() - 1];
  if (!lastHeader->isPepCTermMatch()) {
    change_list.push_back(
        ChangePtr(new Change(lastHeader->getMatchLastBpPos()-first, (last + 1) -first,
        UNEXPECTED_CHANGE, lastHeader->getPepCTermShift(), nullptr)));
  }
  return change_list;
}

} /* namespace prot */
