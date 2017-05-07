// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "base/proteoform.hpp"
#include "base/file_util.hpp"
#include "spec/rm_break_type.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/extend_ms_factory.hpp"
#include "prsm/peak_ion_pair_util.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_coverage.hpp"

namespace prot {

PrsmCoverage::PrsmCoverage(PrsmParaPtr prsm_para_ptr,
                           const std::string &input_file_ext,
                           const std::string &output_file_ext): 
    prsm_para_ptr_(prsm_para_ptr), 
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext) {
    }

void PrsmCoverage::processSingleCoverage(){
  std::string sp_file_name = prsm_para_ptr_->getSpectrumFileName();
  std::string input_file_name = FileUtil::basename(sp_file_name)+"." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileName();
  FastaIndexReaderPtr seq_reader(new FastaIndexReader(db_file_name));
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr_->getFixModPtrVec();
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  //init variables
  int spectrum_num = MsAlignUtil::getSpNum (sp_file_name);
  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          prsm_para_ptr_->getSpParaPtr()->getActivationPtr());
  int cnt = 0;
  SpectrumSetPtr spec_set_ptr;

  std::string output_file_name = FileUtil::basename(sp_file_name)+"."+output_file_ext_;
  std::ofstream out_stream; 
  out_stream.open(output_file_name.c_str());
  //write title
  printTitle(out_stream);
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    cnt+= group_spec_num;
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec 
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        processOnePrsm(out_stream, prsm_ptr, prsm_para_ptr_);
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
      }
    }
    std::cout << std::flush <<  "PrSM coverage is processing " << cnt 
        << " of " << spectrum_num << " spectra.\r";

  }
  LOG_DEBUG("Search completed");
  sp_reader.close();
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
  double n_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double c_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double both_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  int one_third = len /3;
  end = one_third;
  double left_n_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double left_c_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double left_both_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  begin = one_third + 1;
  int two_thirds = len/3 * 2;
  end = two_thirds;
  double middle_n_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double middle_c_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double middle_both_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::BOTH);
  begin = two_thirds + 1;
  end = len -1;
  double right_n_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::N_TERM);
  double right_c_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::C_TERM);
  double right_both_full_coverage = PeakIonPairUtil::computePairConverage(pair_ptrs, begin, end, RmBreakType::BOTH);
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
    spec_ids = spec_ids + std::to_string(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
  }
  boost::algorithm::trim(spec_ids);
  boost::algorithm::trim(spec_activations);
  boost::algorithm::trim(spec_scans);
  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm_ptr->getPrsmId() << "\t"
      << spec_ids << "\t"
      << spec_activations << "\t"
      << spec_scans << "\t"
      << peak_num << "\t"
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"//"Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getSpeciesId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSeqName() << "\t"
      << prsm_ptr->getProteoformPtr()->getStartPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getEndPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm_ptr->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED) << "\t"
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
    spec_ids = spec_ids + std::to_string(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
  }
  boost::algorithm::trim(spec_ids);
  boost::algorithm::trim(spec_activations);
  boost::algorithm::trim(spec_scans);

  file << prsm_para_ptr_->getSpectrumFileName() << "\t"
      << prsm_ptr->getPrsmId() << "\t"
      << spec_ids << "\t"
      << spec_activations << "\t"
      << spec_scans << "\t"
      << peak_num << "\t"
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"//"Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getSpeciesId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSeqName() << "\t"
      << prsm_ptr->getProteoformPtr()->getStartPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getEndPos() << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << prsm_ptr->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED) << "\t"
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

void PrsmCoverage::processOnePrsm(std::ofstream &file, PrsmPtr prsm_ptr, 
                                  PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs =  PeakIonPairFactory::genePeakIonPairs (prsm_ptr->getProteoformPtr(), 
                                                                       prsm_ptr->getRefineMsPtrVec(), min_mass);
  compOneCoverage(file, prsm_ptr, pair_ptrs, prsm_para_ptr);
}

void PrsmCoverage::processTwoPrsms(std::ofstream &file, PrsmPtr prsm_ptr_1, PrsmPtr prsm_ptr_2, 
                                   PrsmParaPtr prsm_para_ptr) {
  double min_mass = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  PeakIonPairPtrVec pair_ptrs_11;
  if (prsm_ptr_1 != nullptr) {
    pair_ptrs_11 =  PeakIonPairFactory::genePeakIonPairs (prsm_ptr_1->getProteoformPtr(), 
                                                          prsm_ptr_1->getRefineMsPtrVec(), min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_12;
  if (prsm_ptr_1 != nullptr && prsm_ptr_2 != nullptr) {
    pair_ptrs_12 =  PeakIonPairFactory::genePeakIonPairs (prsm_ptr_1->getProteoformPtr(), 
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
    pair_ptrs_21 =  PeakIonPairFactory::genePeakIonPairs (prsm_ptr_2->getProteoformPtr(), 
                                                          prsm_ptr_1->getRefineMsPtrVec(), min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_22;
  if (prsm_ptr_2 != nullptr) {
    pair_ptrs_22 =  PeakIonPairFactory::genePeakIonPairs (prsm_ptr_2->getProteoformPtr(), 
                                                          prsm_ptr_2->getRefineMsPtrVec(), min_mass);
  }
  PeakIonPairPtrVec pair_ptrs_2;
  pair_ptrs_2.insert(pair_ptrs_2.begin(), pair_ptrs_21.begin(), pair_ptrs_21.end());
  pair_ptrs_2.insert(pair_ptrs_2.begin(), pair_ptrs_22.begin(), pair_ptrs_22.end());
  if (prsm_ptr_2 != nullptr) {
    compTwoCoverage(file, prsm_ptr_2, pair_ptrs_21, pair_ptrs_22, pair_ptrs_2, prsm_para_ptr);
  }
}

/*
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
   MsAlignReader reader (spectrum_file_name, prsm_para_ptr_->getGroupSpecNum());

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
   */

}

