//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "seq/proteoform.hpp"
#include "ms/spec/rm_break_type.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_coverage.hpp"

namespace toppic {

PrsmCoverage::PrsmCoverage(PrsmParaPtr prsm_para_ptr, 
                           const std::string &input_file_ext,
                           const std::string &output_file_ext):
    prsm_para_ptr_(prsm_para_ptr), 
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext){}

void PrsmCoverage::processSingleCoverage() {
  std::string sp_file_name = prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = file_util::basename(sp_file_name)+"." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileName();
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr_->getFixModPtrVec();
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  // init variables
  int spectrum_num = msalign_util::getSpNum(sp_file_name);
  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();
  SimpleMsAlignReaderPtr ms_reader_ptr = std::make_shared<SimpleMsAlignReader>(sp_file_name, group_spec_num,
                                                                               prsm_para_ptr_->getSpParaPtr()->getActivationPtr());
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  std::string output_file_name = file_util::basename(sp_file_name)+"."+output_file_ext_;
  std::ofstream out_stream;
  out_stream.open(output_file_name.c_str());
  // write title
  printTitle(out_stream);
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  while ((spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr)) !=  nullptr) {
    cnt+= group_spec_num;
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        processOnePrsm(out_stream, prsm_ptr, prsm_para_ptr_);
        outputMatchPeaks(prsm_ptr, prsm_para_ptr_);
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
      }
    }
    std::cout << std::flush <<  "PrSM coverage is processing " << cnt
        << " of " << spectrum_num << " spectra.\r";
  }
  LOG_DEBUG("Search completed");
  prsm_reader.close();
  out_stream.close();
  std::cout << std::endl;
}

void PrsmCoverage::printTitle(std::ofstream &file) {
  file << "Data_file_name" << "\t"
      << "Prsm_ID" << "\t"
      << "Spectrum_ID"<< "\t"
      << "Activation_type" << "\t"
      << "Scan(s)" << "\t"
      << "#peaks"<< "\t"
      << "Charge" << "\t"
      << "Precursor_mass" << "\t"
      << "Adjusted_precursor_mass" << "\t"
      << "Species_ID" << "\t"
      << "Protein_name" << "\t"
      << "First_residue" << "\t"
      << "Last_residue" << "\t"
      << "Peptide" << "\t"
      << "#unexpected_modifications" << "\t"
      << "#matched_peaks" << "\t"
      << "#matched_fragment_ions" << "\t"
      << "P-Value" << "\t"
      << "E-Value" << "\t"
      << "One_Protein_probabilty"<< "\t"
      << "FDR" << "\t"

      << "N_term_ion_coverage" << "\t"
      << "C_term_ion_coverage" << "\t"
      << "Both_term_ion_coverage" << "\t"
      << "left_N_term_ion_coverage" << "\t"
      << "left_C_term_ion_coverage" << "\t"
      << "left_Both_term_ion_coverage" << "\t"
      << "middle_N_term_ion_coverage" << "\t"
      << "middle_C_term_ion_coverage" << "\t"
      << "middle_Both_term_ion_coverage" << "\t"
      << "right_N_term_ion_coverage" << "\t"
      << "right_C_term_ion_coverage" << "\t"
      << "right_Both_term_ion_coverage" << "\t"
      << std::endl;
}

void PrsmCoverage::printTwoTitle(std::ofstream &file) {
  file << "Data_file_name" << "\t"
      << "Prsm_ID" << "\t"
      << "Spectrum_ID"<< "\t"
      << "Activation_type" << "\t"
      << "Scan(s)" << "\t"
      << "#peaks"<< "\t"
      << "Charge" << "\t"
      << "Precursor_mass" << "\t"
      << "Adjusted_precursor_mass" << "\t"
      << "Species_ID" << "\t"
      << "Protein_name" << "\t"
      << "First_residue" << "\t"
      << "Last_residue" << "\t"
      << "Peptide" << "\t"
      << "#unexpected_modifications" << "\t"
      << "#matched_peaks" << "\t"
      << "#matched_fragment_ions" << "\t"
      << "P-Value" << "\t"
      << "E-Value" << "\t"
      << "One_Protein_probabilty"<< "\t"
      << "FDR" << "\t"

      << "N_term_ion_coverage_1" << "\t"
      << "C_term_ion_coverage_1" << "\t"
      << "Both_term_ion_coverage_1" << "\t"
      << "left_N_term_ion_coverage_1" << "\t"
      << "left_C_term_ion_coverage_1" << "\t"
      << "left_Both_term_ion_coverage_1" << "\t"
      << "middle_N_term_ion_coverage_1" << "\t"
      << "middle_C_term_ion_coverage_1" << "\t"
      << "middle_Both_term_ion_coverage_1" << "\t"
      << "right_N_term_ion_coverage_1" << "\t"
      << "right_C_term_ion_coverage_1" << "\t"
      << "right_Both_term_ion_coverage_1" << "\t"

      << "N_term_ion_coverage_2" << "\t"
      << "C_term_ion_coverage_2" << "\t"
      << "Both_term_ion_coverage_2" << "\t"
      << "left_N_term_ion_coverage_2" << "\t"
      << "left_C_term_ion_coverage_2" << "\t"
      << "left_Both_term_ion_coverage_2" << "\t"
      << "middle_N_term_ion_coverage_2" << "\t"
      << "middle_C_term_ion_coverage_2" << "\t"
      << "middle_Both_term_ion_coverage_2" << "\t"
      << "right_N_term_ion_coverage_2" << "\t"
      << "right_C_term_ion_coverage_2" << "\t"
      << "right_Both_term_ion_coverage_2" << "\t"

      << "N_term_ion_coverage_3" << "\t"
      << "C_term_ion_coverage_3" << "\t"
      << "Both_term_ion_coverage_3" << "\t"
      << "left_N_term_ion_coverage_3" << "\t"
      << "left_C_term_ion_coverage_3" << "\t"
      << "left_Both_term_ion_coverage_3" << "\t"
      << "middle_N_term_ion_coverage_3" << "\t"
      << "middle_C_term_ion_coverage_3" << "\t"
      << "middle_Both_term_ion_coverage_3" << "\t"
      << "right_N_term_ion_coverage_3" << "\t"
      << "right_C_term_ion_coverage_3" << "\t"
      << "right_Both_term_ion_coverage_3" << "\t"

      << std::endl;
}

void PrsmCoverage::computeCoverage(std::ofstream &file, PrsmPtr prsm_ptr,
                                   PeakIonPairPtrVec &pair_ptrs, PrsmParaPtr prsm_para_ptr) {
  int len = prsm_ptr->getProteoformPtr()->getResSeqPtr()->getLen() - 1;
  int begin = 1;
  int end = len - 1;
  double n_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double c_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double both_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  int one_third = len /3;
  end = one_third;
  double left_n_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double left_c_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double left_both_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  begin = one_third + 1;
  int two_thirds = len/3 * 2;
  end = two_thirds;
  double middle_n_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double middle_c_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double middle_both_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  begin = two_thirds + 1;
  end = len -1;
  double right_n_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double right_c_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double right_both_full_coverage = peak_ion_pair_util::computePairCoverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  file << n_full_coverage << "\t"
      << c_full_coverage << "\t"
      << both_full_coverage << "\t"
      << left_n_full_coverage << "\t"
      << left_c_full_coverage << "\t"
      << left_both_full_coverage << "\t"
      << middle_n_full_coverage << "\t"
      << middle_c_full_coverage << "\t"
      << middle_both_full_coverage << "\t"
      << right_n_full_coverage << "\t"
      << right_c_full_coverage << "\t"
      << right_both_full_coverage << "\t";
}

void PrsmCoverage::compOneCoverage(std::ofstream &file, PrsmPtr prsm_ptr,
                                   PeakIonPairPtrVec &pair_ptrs, PrsmParaPtr prsm_para_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  int peak_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
  }
  str_util::trim(spec_ids);
  str_util::trim(spec_activations);
  str_util::trim(spec_scans);
  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm_ptr->getPrsmId() << "\t"
      << spec_ids << "\t"
      << spec_activations << "\t"
      << spec_scans << "\t"
      << peak_num << "\t"
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"  // "Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteoClusterId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSeqName() << "\t"
      << prsm_ptr->getProteoformPtr()->getStartPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getEndPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm_ptr->getProteoformPtr()->getMassShiftNum(AlterType::UNEXPECTED) << "\t"
      << prsm_ptr->getMatchPeakNum() << "\t"
      << prsm_ptr->getMatchFragNum() << "\t"
      << prsm_ptr->getPValue() << "\t"
      << prsm_ptr->getEValue() << "\t"
      << prsm_ptr->getOneProtProb()<< "\t"
      << prsm_ptr->getFdr() << "\t";
  computeCoverage(file, prsm_ptr, pair_ptrs, prsm_para_ptr);
  file << std::endl;
}

void PrsmCoverage::compTwoCoverage(std::ofstream &file, PrsmPtr prsm_ptr,
                                   PeakIonPairPtrVec &pair_ptrs_1,
                                   PeakIonPairPtrVec &pair_ptrs_2,
                                   PeakIonPairPtrVec &pair_ptrs_3,
                                   PrsmParaPtr prsm_para_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  int peak_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
  }
  str_util::trim(spec_ids);
  str_util::trim(spec_activations);
  str_util::trim(spec_scans);

  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm_ptr->getPrsmId() << "\t"
      << spec_ids << "\t"
      << spec_activations << "\t"
      << spec_scans << "\t"
      << peak_num << "\t"
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"  // "Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteoClusterId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSeqName() << "\t"
      << prsm_ptr->getProteoformPtr()->getStartPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getEndPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm_ptr->getProteoformPtr()->getMassShiftNum(AlterType::UNEXPECTED) << "\t"
      << prsm_ptr->getMatchPeakNum() << "\t"
      << prsm_ptr->getMatchFragNum() << "\t"
      << prsm_ptr->getPValue() << "\t"
      << prsm_ptr->getEValue() << "\t"
      << prsm_ptr->getOneProtProb()<< "\t"
      << prsm_ptr->getFdr() << "\t";
  computeCoverage(file, prsm_ptr, pair_ptrs_1, prsm_para_ptr);
  computeCoverage(file, prsm_ptr, pair_ptrs_2, prsm_para_ptr);
  computeCoverage(file, prsm_ptr, pair_ptrs_3, prsm_para_ptr);
  file << std::endl;
}

void PrsmCoverage::outputMatchPeaks(PrsmPtr prsm_ptr, PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs = peak_ion_pair_util::genePeakIonPairs(prsm_ptr->getProteoformPtr(),
                                                                     prsm_ptr->getRefineMsPtrVec(),
                                                                     min_mass);
  std::string file_name = "scan_"
      + prsm_ptr->getDeconvMsPtrVec()[0]->getMsHeaderPtr()->getScansString() + ".peak";
  std::ofstream ps(file_name, std::ofstream::out);
  for (size_t i = 0; i < pair_ptrs.size(); i++) {
    ExtendPeakPtr ex_peak = pair_ptrs[i]->getRealPeakPtr();
    TheoPeakPtr theo_peak = pair_ptrs[i]->getTheoPeakPtr();
    ps << ex_peak->getMonoMass() << "\t"
        << ex_peak->getBasePeakPtr()->getMonoMass() << "\t"
        << theo_peak->getModMass() << std::endl;
  }
  ps.close();
}

void PrsmCoverage::processOnePrsm(std::ofstream &file, PrsmPtr prsm_ptr,
                                  PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs = peak_ion_pair_util::genePeakIonPairs(prsm_ptr->getProteoformPtr(),
                                                                     prsm_ptr->getRefineMsPtrVec(),
                                                                     min_mass);
  compOneCoverage(file, prsm_ptr, pair_ptrs, prsm_para_ptr);
}

void PrsmCoverage::processTwoPrsms(std::ofstream &file, PrsmPtr prsm_ptr_1, PrsmPtr prsm_ptr_2,
                                   PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs_11;
  if (prsm_ptr_1 != nullptr) {
    pair_ptrs_11 = peak_ion_pair_util::genePeakIonPairs(prsm_ptr_1->getProteoformPtr(),
                                                        prsm_ptr_1->getRefineMsPtrVec(), min_mass);
  }

  PeakIonPairPtrVec pair_ptrs_12;
  if (prsm_ptr_1 != nullptr && prsm_ptr_2 != nullptr) {
    pair_ptrs_12 = peak_ion_pair_util::genePeakIonPairs(prsm_ptr_1->getProteoformPtr(),
                                                        prsm_ptr_2->getRefineMsPtrVec(), min_mass);
  }

  PeakIonPairPtrVec pair_ptrs_1;
  pair_ptrs_1.insert(pair_ptrs_1.begin(), pair_ptrs_11.begin(), pair_ptrs_11.end());
  pair_ptrs_1.insert(pair_ptrs_1.begin(), pair_ptrs_12.begin(), pair_ptrs_12.end());
  if (prsm_ptr_1 != nullptr) {
    compTwoCoverage(file, prsm_ptr_1, pair_ptrs_11, pair_ptrs_12, pair_ptrs_1, prsm_para_ptr);
  }

  PeakIonPairPtrVec pair_ptrs_21;
  if (prsm_ptr_1 != nullptr && prsm_ptr_2 != nullptr) {
    pair_ptrs_21 = peak_ion_pair_util::genePeakIonPairs(prsm_ptr_2->getProteoformPtr(),
                                                        prsm_ptr_1->getRefineMsPtrVec(), min_mass);
  }

  PeakIonPairPtrVec pair_ptrs_22;
  if (prsm_ptr_2 != nullptr) {
    pair_ptrs_22 = peak_ion_pair_util::genePeakIonPairs(prsm_ptr_2->getProteoformPtr(),
                                                        prsm_ptr_2->getRefineMsPtrVec(), min_mass);
  }

  PeakIonPairPtrVec pair_ptrs_2;
  pair_ptrs_2.insert(pair_ptrs_2.begin(), pair_ptrs_21.begin(), pair_ptrs_21.end());
  pair_ptrs_2.insert(pair_ptrs_2.begin(), pair_ptrs_22.begin(), pair_ptrs_22.end());
  if (prsm_ptr_2 != nullptr) {
    compTwoCoverage(file, prsm_ptr_2, pair_ptrs_21, pair_ptrs_22, pair_ptrs_2, prsm_para_ptr);
  }
}
}  // namespace toppic

