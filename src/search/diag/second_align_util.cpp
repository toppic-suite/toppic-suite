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

#include "search/diag/second_align_util.hpp"

#include <cstdlib>
#include <fstream>
#include <limits>
#include <mutex>
#include <sstream>

#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/prm_ms_factory.hpp"
#include "ms/factory/prm_ms_util.hpp"
#include "search/diag/diag_header_util.hpp"
#include "search/diag/diag_pair_util.hpp"
#include "search/diag/diagonal_util.hpp"
#include "seq/alter_type.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace second_align_util {

PrsmPtr genePrsm(const ProteoformPtr& proteo_ptr,
                 const DiagHeaderPtrVec& header_ptrs,
                 DeconvMsPtrVec& deconv_ms_ptr_vec,
                 const ExtendMsPtrVec& ms_three_ptr_vec,
                 const PrsmParaPtr& prsm_para_ptr, double& refine_prec_mass,
                 ExtendMsPtrVec& refine_ms_ptr_vec, double clamp_prec_mass) {
  int first_pos = header_ptrs[0]->getTruncFirstResPos();
  int last_pos = header_ptrs[header_ptrs.size() - 1]->getTruncLastResPos();
  ProteoformPtr sub_proteo_ptr = proteoform_factory::geneSubProteoform(
      proteo_ptr, proteo_ptr->getFastaSeqPtr(), first_pos, last_pos);

  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  double min_mass = sp_para_ptr->getMinMass();
  PeakTolerancePtr prec_tole_ptr = sp_para_ptr->getPrecTolerancePtr();

  // the refinement moves the headers' shifts
  DiagHeaderPtrVec headers = header_ptrs;
  double delta_min = -std::numeric_limits<double>::infinity();
  double delta_max = std::numeric_limits<double>::infinity();
  if (clamp_prec_mass > 0) {
    // the delta that keeps the refined precursor within the tolerance of
    // the observed one
    double prec_mass =
        ms_three_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
    double tole = prec_tole_ptr->compPrecErrorTole(clamp_prec_mass);
    delta_min = clamp_prec_mass - tole - prec_mass;
    delta_max = clamp_prec_mass + tole - prec_mass;
  }
  refine_prec_mass = diagonal_util::refinePrecursorByMatchedIons(
      proteo_ptr, ms_three_ptr_vec, headers, prec_tole_ptr, min_mass,
      delta_min, delta_max);

  refine_ms_ptr_vec = extend_ms_factory::geneMsThreePtrVec(
      deconv_ms_ptr_vec, sp_para_ptr, refine_prec_mass);

  DiagHeaderPtrVec refined_header_ptrs = diagonal_util::refineHeadersBgnEnd(
      proteo_ptr, refine_ms_ptr_vec, headers, min_mass);

  if (refined_header_ptrs.size() == 0) {
    return nullptr;
  }

  MassShiftPtrVec shifts = diag_header_util::getDiagonalMassChanges(
      refined_header_ptrs, first_pos, last_pos, AlterType::UNEXPECTED);
  sub_proteo_ptr->addMassShiftPtrVec(shifts);

  return std::make_shared<Prsm>(sub_proteo_ptr, deconv_ms_ptr_vec,
                                refine_prec_mass, sp_para_ptr);
}

PrmPeakPtrVec genePrmPeaks(const DeconvMsPtrVec& deconv_ms_ptr_vec,
                           const SpParaPtr& sp_para_ptr,
                           double refine_prec_mass, bool strict) {
  // the PRM factory takes the observed (labelled) precursor mass and the
  // label, and builds the peaks for the label-free mass
  double n_term_label_mass = sp_para_ptr->getNTermLabelMass();
  PrmMsPtrVec ms_six_ptr_vec = prm_ms_factory::geneMsSixPtrVec(
      deconv_ms_ptr_vec, sp_para_ptr, refine_prec_mass + n_term_label_mass,
      n_term_label_mass);
  if (strict) {
    for (size_t i = 0; i < ms_six_ptr_vec.size(); i++) {
      for (size_t j = 0; j < ms_six_ptr_vec[i]->size(); j++) {
        PrmPeakPtr peak_ptr = ms_six_ptr_vec[i]->getPeakPtr(j);
        double tole = peak_ptr->getStrictTolerance();
        peak_ptr->setNStrictCRelaxTolerance(tole);
        peak_ptr->setNRelaxCStrictTolerance(tole);
      }
    }
  }
  // the zero and precursor peaks, sized on the refined precursor
  return prm_ms_util::getPrmPeakPtrs(ms_six_ptr_vec,
                                     sp_para_ptr->getPrecTolerancePtr());
}

DiagHeaderPtrVec moveHeaders(const DiagonalPtrVec& diagonal_ptrs, double delta,
                             double refine_prec_mass,
                             const ProteoformPtr& proteo_ptr) {
  double seq_mass = proteo_ptr->getResSeqPtr()->getSeqMass();
  DiagHeaderPtrVec result;
  for (size_t i = 0; i < diagonal_ptrs.size(); i++) {
    DiagHeaderPtr header_ptr = diagonal_ptrs[i]->getHeader();
    double n_shift = header_ptr->getProtNTermShift();
    if (!header_ptr->isNStrict()) {
      // anchored at the C-terminus: the C-terminal shift is kept
      n_shift += delta;
    }
    DiagHeaderPtr moved_ptr = std::make_shared<DiagHeader>(
        n_shift, header_ptr->isNStrict(), header_ptr->isCStrict(),
        header_ptr->isProtNTermMatch(), header_ptr->isProtCTermMatch(),
        header_ptr->isPepNTermMatch(), header_ptr->isPepCTermMatch());
    moved_ptr->initHeader(refine_prec_mass - seq_mass - n_shift, proteo_ptr);
    result.push_back(moved_ptr);
  }
  return result;
}

DiagonalPtrVec geneDiagonals(const DiagHeaderPtrVec& header_ptrs,
                             const PrmPeakPtrVec& prm_peaks,
                             int group_spec_num,
                             const ProteoformPtr& proteo_ptr,
                             std::vector<double>& ms_masses) {
  ms_masses.resize(prm_peaks.size());
  for (size_t i = 0; i < prm_peaks.size(); i++) {
    ms_masses[i] = prm_peaks[i]->getPosition();
  }
  return diag_pair_util::geneDiagonalsWithoutEmptyList(
      header_ptrs, prm_peaks, group_spec_num, proteo_ptr);
}

namespace {

std::mutex log_mutex;

const char* logFileName() {
  static const char* name = std::getenv("TOPPIC_SECOND_ALIGN_LOG");
  return name;
}

}  // namespace

bool logEnabled() { return logFileName() != nullptr; }

std::string prsmInfo(const PrsmPtr& prsm_ptr) {
  if (prsm_ptr == nullptr) {
    return "-\t-\t-\t-";
  }
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();
  std::ostringstream out;
  out.precision(4);
  out << std::fixed;
  out << form_ptr->getStartPos() << "\t" << form_ptr->getEndPos() << "\t"
      << prsm_ptr->getMatchTheoFragNum() << "\t";
  MassShiftPtrVec shifts = form_ptr->getMassShiftPtrVec(AlterType::UNEXPECTED);
  if (shifts.size() == 0) {
    out << "-";
  }
  for (size_t i = 0; i < shifts.size(); i++) {
    if (i > 0) {
      out << ";";
    }
    out << shifts[i]->getLeftBpPos() << ":" << shifts[i]->getRightBpPos()
        << ":" << shifts[i]->getMassShift();
  }
  return out.str();
}

void log(const std::string& line) {
  const char* name = logFileName();
  if (name == nullptr) {
    return;
  }
  std::lock_guard<std::mutex> guard(log_mutex);
  std::ofstream out(name, std::ios::app);
  out << line << "\n";
}

}  // namespace second_align_util

}  // namespace toppic
