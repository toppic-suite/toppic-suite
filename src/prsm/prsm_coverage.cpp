#include <iostream>
#include <fstream>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_coverage.hpp"

namespace prot {

PrsmCoverage::PrsmCoverage(PrsmParaPtr prsm_para_ptr,
                           const std::string &input_file_ext,
                           const std::string &output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
}

void PrsmCoverage::processSingleCoverage(){

  ProteoformPtrVec raw_form_ptrs 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string base_name = basename(prsm_para_ptr_->getSpectrumFileName());
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmPtrVec prsm_ptrs = readPrsm(input_file_name, raw_form_ptrs);
  LOG_DEBUG("read prsm_ptr complete ");
  addSpectrumPtrsToPrsms(prsm_ptrs, prsm_para_ptr_);
  LOG_DEBUG("prsm_ptrs loaded");

  std::string output_file_name = base_name+"."+output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  //write title
  printTitle(file);

  for(size_t i=0;i<prsm_ptrs.size();i++){
    processOnePrsm(file, prsm_ptrs[i], prsm_para_ptr_);
  }
  file.close();
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
      << "Protein_ID" << "\t"
      << "Species_ID" << "\t"
      << "Protein_name" << "\t"
      << "Protein_mass" << "\t"
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
      << "Protein_ID" << "\t"
      << "Species_ID" << "\t"
      << "Protein_name" << "\t"
      << "Protein_mass" << "\t"
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
  double n_full_coverage = computePairConverage(pair_ptrs, begin, end, N_TERM_COVERAGE);
  double c_full_coverage = computePairConverage(pair_ptrs, begin, end, C_TERM_COVERAGE);
  double both_full_coverage = computePairConverage(pair_ptrs, begin, end, BOTH_TERM_COVERAGE);
  int one_third = len /3;
  end = one_third;
  double left_n_full_coverage = computePairConverage(pair_ptrs, begin, end, N_TERM_COVERAGE);
  double left_c_full_coverage = computePairConverage(pair_ptrs, begin, end, C_TERM_COVERAGE);
  double left_both_full_coverage = computePairConverage(pair_ptrs, begin, end, BOTH_TERM_COVERAGE);
  begin = one_third + 1;
  int two_thirds = len/3 * 2;
  end = two_thirds;
  double middle_n_full_coverage = computePairConverage(pair_ptrs, begin, end, N_TERM_COVERAGE);
  double middle_c_full_coverage = computePairConverage(pair_ptrs, begin, end, C_TERM_COVERAGE);
  double middle_both_full_coverage = computePairConverage(pair_ptrs, begin, end, BOTH_TERM_COVERAGE);
  begin = two_thirds + 1;
  end = len -1;
  double right_n_full_coverage = computePairConverage(pair_ptrs, begin, end, N_TERM_COVERAGE);
  double right_c_full_coverage = computePairConverage(pair_ptrs, begin, end, C_TERM_COVERAGE);
  double right_both_full_coverage = computePairConverage(pair_ptrs, begin, end, BOTH_TERM_COVERAGE);
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

  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm_ptr->getId() << "\t"
      << prsm_ptr->getSpectrumId()<< "\t"
      << prsm_ptr->getDeconvMsPtr()->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
      << prsm_ptr->getSpectrumScan() << "\t"
      << prsm_ptr->getDeconvMsPtr()->size()<< "\t"
      << prsm_ptr->getDeconvMsPtr()->getHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"//"Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSpeciesId() << "\t"
      << prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getName() << "\t"
      << prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getSeqMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getStartPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getEndPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm_ptr->getProteoformPtr()->getUnexpectedChangeNum() << "\t"
      << prsm_ptr->getMatchPeakNum() << "\t"
      << prsm_ptr->getMatchFragNum() << "\t"
      << prsm_ptr->getPValue() << "\t"
      << prsm_ptr->getEValue() << "\t"
      << prsm_ptr->getProbPtr()->getOneProtProb()<< "\t"
      << prsm_ptr->getFdr() << "\t";
  computeCoverage(file, prsm_ptr, pair_ptrs, prsm_para_ptr);
  file << std::endl;
}


void PrsmCoverage::compTwoCoverage(std::ofstream &file, PrsmPtr prsm_ptr,  
                                   PeakIonPairPtrVec &pair_ptrs_1, 
                                   PeakIonPairPtrVec &pair_ptrs_2, 
                                   PeakIonPairPtrVec &pair_ptrs_3, 
                                   PrsmParaPtr prsm_para_ptr) {
  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm_ptr->getId() << "\t"
      << prsm_ptr->getSpectrumId()<< "\t"
      << prsm_ptr->getDeconvMsPtr()->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
      << prsm_ptr->getSpectrumScan() << "\t"
      << prsm_ptr->getDeconvMsPtr()->size()<< "\t"
      << prsm_ptr->getDeconvMsPtr()->getHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"//"Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSpeciesId() << "\t"
      << prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getName() << "\t"
      << prsm_ptr->getProteoformPtr()->getDbResSeqPtr()->getSeqMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getStartPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getEndPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm_ptr->getProteoformPtr()->getUnexpectedChangeNum() << "\t"
      << prsm_ptr->getMatchPeakNum() << "\t"
      << prsm_ptr->getMatchFragNum() << "\t"
      << prsm_ptr->getPValue() << "\t"
      << prsm_ptr->getEValue() << "\t"
      << prsm_ptr->getProbPtr()->getOneProtProb()<< "\t"
      << prsm_ptr->getFdr() << "\t";
  computeCoverage(file, prsm_ptr, pair_ptrs_1, prsm_para_ptr);
  computeCoverage(file, prsm_ptr, pair_ptrs_2, prsm_para_ptr);
  computeCoverage(file, prsm_ptr, pair_ptrs_3, prsm_para_ptr);
  file << std::endl;
}

void PrsmCoverage::processOnePrsm(std::ofstream &file, PrsmPtr prsm_ptr, 
                                  PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs =  getPeakIonPairs (prsm_ptr->getProteoformPtr(), 
                                              prsm_ptr->getRefineMs(),
                                              min_mass);
  compOneCoverage(file, prsm_ptr, pair_ptrs, prsm_para_ptr);
}

void PrsmCoverage::processTwoPrsms(std::ofstream &file, PrsmPtr prsm_ptr_1, PrsmPtr prsm_ptr_2, 
                                  PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs_11;
  if (prsm_ptr_1 != nullptr) {
    pair_ptrs_11 =  getPeakIonPairs (prsm_ptr_1->getProteoformPtr(), 
                                 prsm_ptr_1->getRefineMs(),
                                 min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_12;
  if (prsm_ptr_1 != nullptr && prsm_ptr_2 != nullptr) {
    pair_ptrs_12 =  getPeakIonPairs (prsm_ptr_1->getProteoformPtr(), 
                                 prsm_ptr_2->getRefineMs(),
                                 min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_1;
  pair_ptrs_1.insert(pair_ptrs_1.begin(), pair_ptrs_11.begin(), pair_ptrs_11.end());
  pair_ptrs_1.insert(pair_ptrs_1.begin(), pair_ptrs_12.begin(), pair_ptrs_12.end());
  if (prsm_ptr_1 != nullptr) {
    compTwoCoverage(file, prsm_ptr_1, pair_ptrs_11, pair_ptrs_12, pair_ptrs_1, prsm_para_ptr);
  }

  PeakIonPairPtrVec pair_ptrs_21;
  if (prsm_ptr_1 != nullptr && prsm_ptr_2 != nullptr) {
    pair_ptrs_21 =  getPeakIonPairs (prsm_ptr_2->getProteoformPtr(), 
                                 prsm_ptr_1->getRefineMs(),
                                 min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_22;
  if (prsm_ptr_2 != nullptr) {
    pair_ptrs_22 =  getPeakIonPairs (prsm_ptr_2->getProteoformPtr(), 
                                 prsm_ptr_2->getRefineMs(),
                                 min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_2;
  pair_ptrs_2.insert(pair_ptrs_2.begin(), pair_ptrs_21.begin(), pair_ptrs_21.end());
  pair_ptrs_2.insert(pair_ptrs_2.begin(), pair_ptrs_22.begin(), pair_ptrs_22.end());
  if (prsm_ptr_2 != nullptr) {
    compTwoCoverage(file, prsm_ptr_2, pair_ptrs_21, pair_ptrs_22, pair_ptrs_2, prsm_para_ptr);
  }
}

void PrsmCoverage::processCombineCoverage(){

  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string base_name = basename(prsm_para_ptr_->getSpectrumFileName());
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmPtrVec prsm_ptrs = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("read prsm_ptr complete ");
  addSpectrumPtrsToPrsms(prsm_ptrs, prsm_para_ptr_);
  LOG_DEBUG("prsm_ptrs loaded");

  std::string output_file_name = base_name+".COMBINE_"+output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  printTwoTitle(file);

  std::string spectrum_file_name = prsm_para_ptr_->getSpectrumFileName();
  MsAlignReader reader (spectrum_file_name);

  DeconvMsPtr ms_ptr_1 = reader.getNextMs();
  DeconvMsPtr ms_ptr_2 = reader.getNextMs();
  while (ms_ptr_1.get() != nullptr && ms_ptr_2.get() != nullptr) {
    PrsmPtrVec sele_prsm_ptrs_1;
    PrsmPtrVec sele_prsm_ptrs_2;
    filterPrsms(prsm_ptrs, ms_ptr_1->getHeaderPtr(), sele_prsm_ptrs_1);
    filterPrsms(prsm_ptrs, ms_ptr_2->getHeaderPtr(), sele_prsm_ptrs_2);
    if (sele_prsm_ptrs_1.size() >= 1 || sele_prsm_ptrs_2.size() >= 1) {
      if (sele_prsm_ptrs_1.size() == 0) {
        processTwoPrsms(file, PrsmPtr(nullptr), sele_prsm_ptrs_2[0], prsm_para_ptr_);
      }
      else {
        if (sele_prsm_ptrs_2.size() == 0) {
          processTwoPrsms(file, sele_prsm_ptrs_1[0], PrsmPtr(nullptr), prsm_para_ptr_);
        }
        else {
          processTwoPrsms(file, sele_prsm_ptrs_1[0], sele_prsm_ptrs_2[0], prsm_para_ptr_);
        }
      }
    }
    ms_ptr_1 = reader.getNextMs();
    ms_ptr_2 = reader.getNextMs();
  }
  reader.close();
  file.close();
}

}
    
