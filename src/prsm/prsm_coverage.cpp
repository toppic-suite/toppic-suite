#include <iostream>
#include <fstream>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_coverage.hpp"

namespace prot {

PrsmCoverage::PrsmCoverage(PrsmParaPtr prsm_para_ptr,
                           std::string input_file_ext,
                           std::string output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
}

void PrsmCoverage::processSingleCoverage(){

  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string base_name = basename(prsm_para_ptr_->getSpectrumFileName());
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmPtrVec prsms = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("read prsm complete ");
  addSpectrumPtrsToPrsms(prsms, prsm_para_ptr_);
  LOG_DEBUG("prsms loaded");

  std::string output_file_name = base_name+"."+output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  //write title
  printTitle(file);

  for(unsigned int i=0;i<prsms.size();i++){
    processOnePrsm(file, prsms[i], prsm_para_ptr_);
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

void PrsmCoverage::compCoverage(std::ofstream &file, PrsmPtr prsm, 
                                PeakIonPairPtrVec &pairs, PrsmParaPtr prsm_para_ptr) {
  int len = prsm->getProteoformPtr()->getResSeqPtr()->getLen() - 1;
  int begin = 1;
  int end = len - 1;
  double n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
  double c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
  double both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);
  int one_third = len /3;
  end = one_third;
  double left_n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
  double left_c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
  double left_both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);
  begin = one_third + 1;
  int two_thirds = len/3 * 2;
  end = two_thirds;
  double middle_n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
  double middle_c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
  double middle_both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);
  begin = two_thirds + 1;
  end = len -1;
  double right_n_full_coverage = computePairConverage(pairs, begin, end, N_TERM_COVERAGE);
  double right_c_full_coverage = computePairConverage(pairs, begin, end, C_TERM_COVERAGE);
  double right_both_full_coverage = computePairConverage(pairs, begin, end, BOTH_TERM_COVERAGE);

  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm->getId() << "\t"
      << prsm->getSpectrumId()<< "\t"
      << prsm->getDeconvMsPtr()->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
      << prsm->getSpectrumScan() << "\t"
      << prsm->getDeconvMsPtr()->size()<< "\t"
      << prsm->getDeconvMsPtr()->getHeaderPtr()->getPrecCharge() << "\t"
      << prsm->getOriPrecMass()<< "\t"//"Precursor_mass"
      << prsm->getAdjustedPrecMass() << "\t"
      << prsm->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
      << prsm->getProteoformPtr()->getSpeciesId() << "\t"
      << prsm->getProteoformPtr()->getDbResSeqPtr()->getName() << "\t"
      << prsm->getProteoformPtr()->getDbResSeqPtr()->getSeqMass() << "\t"
      << prsm->getProteoformPtr()->getStartPos() << "\t"
      << prsm->getProteoformPtr()->getEndPos() << "\t"
      << prsm->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm->getProteoformPtr()->getUnexpectedChangeNum() << "\t"
      << prsm->getMatchPeakNum() << "\t"
      << prsm->getMatchFragNum() << "\t"
      << prsm->getPValue() << "\t"
      << prsm->getEValue() << "\t"
      << prsm->getProbPtr()->getOneProtProb()<< "\t"
      << prsm->getFdr() << "\t"

      << n_full_coverage << "\t"
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
      << right_both_full_coverage << "\t"
      << std::endl;
}

void PrsmCoverage::processOnePrsm(std::ofstream &file, PrsmPtr prsm, 
                                  PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pairs =  getPeakIonPairs (prsm->getProteoformPtr(), 
                                              prsm->getRefineMs(),
                                              min_mass);
  compCoverage(file, prsm, pairs, prsm_para_ptr);
}

void PrsmCoverage::processTwoPrsms(std::ofstream &file, PrsmPtr prsm_1, PrsmPtr prsm_2, 
                                  PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pairs_11 =  getPeakIonPairs (prsm_1->getProteoformPtr(), 
                                                prsm_1->getRefineMs(),
                                                min_mass);
  PeakIonPairPtrVec pairs_12 =  getPeakIonPairs (prsm_1->getProteoformPtr(), 
                                                prsm_2->getRefineMs(),
                                                min_mass);
  PeakIonPairPtrVec pairs_1;
  pairs_1.insert(pairs_1.begin(), pairs_11.begin(), pairs_11.end());
  pairs_1.insert(pairs_1.begin(), pairs_12.begin(), pairs_12.end());
  compCoverage(file, prsm_1, pairs_1, prsm_para_ptr);

  PeakIonPairPtrVec pairs_21 =  getPeakIonPairs (prsm_2->getProteoformPtr(), 
                                                 prsm_1->getRefineMs(),
                                                 min_mass);
  PeakIonPairPtrVec pairs_22 =  getPeakIonPairs (prsm_2->getProteoformPtr(), 
                                                 prsm_2->getRefineMs(),
                                                 min_mass);
  PeakIonPairPtrVec pairs_2;
  pairs_2.insert(pairs_2.begin(), pairs_21.begin(), pairs_21.end());
  pairs_2.insert(pairs_2.begin(), pairs_22.begin(), pairs_22.end());
  compCoverage(file, prsm_2, pairs_2, prsm_para_ptr);
}

void PrsmCoverage::processCombineCoverage(){

  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string base_name = basename(prsm_para_ptr_->getSpectrumFileName());
  std::string input_file_name = base_name + "." + input_file_ext_;
  PrsmPtrVec prsms = readPrsm(input_file_name, raw_forms);
  LOG_DEBUG("read prsm complete ");
  addSpectrumPtrsToPrsms(prsms, prsm_para_ptr_);
  LOG_DEBUG("prsms loaded");

  std::string output_file_name = base_name+".COMBINE_"+output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  printTitle(file);

  std::string spectrum_file_name = prsm_para_ptr_->getSpectrumFileName();
  MsAlignReader reader (spectrum_file_name);

  DeconvMsPtr ms_ptr_1 = reader.getNextMs();
  DeconvMsPtr ms_ptr_2 = reader.getNextMs();
  while (ms_ptr_1.get() != nullptr && ms_ptr_2.get() != nullptr) {
    PrsmPtrVec sele_prsms_1;
    PrsmPtrVec sele_prsms_2;
    filterPrsms(prsms, ms_ptr_1->getHeaderPtr(), sele_prsms_1);
    filterPrsms(prsms, ms_ptr_2->getHeaderPtr(), sele_prsms_2);
    if (sele_prsms_1.size() == 0 || sele_prsms_2.size() == 0) {
      if (sele_prsms_1.size() == 1) {
        processOnePrsm(file, sele_prsms_1[0], prsm_para_ptr_);
      }

      if (sele_prsms_2.size() == 1) {
        processOnePrsm(file, sele_prsms_2[0], prsm_para_ptr_);
      }
    }
    else {
      processTwoPrsms(file, sele_prsms_1[0], sele_prsms_2[0], prsm_para_ptr_);
    }

    ms_ptr_1 = reader.getNextMs();
    ms_ptr_2 = reader.getNextMs();
  }
  reader.close();
  file.close();
}

}
    
