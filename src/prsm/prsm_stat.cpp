//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>

#include "base/acid_base.hpp"
#include "base/file_util.hpp"
#include "spec/peak.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "prsm/peak_ion_pair_factory.hpp"
#include "prsm/prsm_stat.hpp"

namespace prot {

PrsmStat::PrsmStat(PrsmParaPtr prsm_para_ptr, 
                   const std::string &input_file_ext, 
                   const std::string &output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  min_mass_ = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
  acid_ptr_vec_ = AcidBase::getBaseAcidPtrVec();
}

int countCoverage(const std::vector<bool> &match_ion_vec, int start, int end) {
  int count = 0;
  for (size_t i = 0; i < match_ion_vec.size(); i++) {
    if ((int)i >= start && (int)i < end && match_ion_vec[i]) {
      count++;
    }
  }
  return count;
}

int countAcid(ResSeqPtr res_seq_ptr, AcidPtr acid_ptr) {
  int count = 0;
  for (int i = 0; i < res_seq_ptr->getLen(); i++) {
    if (res_seq_ptr->getResiduePtr(i)->getAcidPtr() == acid_ptr) {
      count++;
    }
  }
  return count;
}

int countAcidLeftCoverage(ResSeqPtr res_seq_ptr, AcidPtr acid_ptr, 
                          const std::vector<bool> &match_ion_vec) {
  int count = 0;
  for (int i = 0; i < res_seq_ptr->getLen(); i++) {
    if (res_seq_ptr->getResiduePtr(i)->getAcidPtr() == acid_ptr && match_ion_vec[i]) {
      count++;
    }
  }
  return count;
}

int countAcidRightCoverage(ResSeqPtr res_seq_ptr, AcidPtr acid_ptr, 
                          const std::vector<bool> &match_ion_vec) {
  int count = 0;
  for (int i = 0; i < res_seq_ptr->getLen(); i++) {
    if (res_seq_ptr->getResiduePtr(i)->getAcidPtr() == acid_ptr && match_ion_vec[i+1]) {
      count++;
    }
  }
  return count;
}

void PrsmStat::writePrsm(std::ofstream &file, PrsmPtr prsm_ptr) {
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
      << spec_activations<< "\t"
      << spec_scans << "\t"
      << peak_num << "\t"
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << "\t"
      << prsm_ptr->getOriPrecMass()<< "\t"//"Precursor_mass"
      << prsm_ptr->getAdjustedPrecMass() << "\t"
      << prsm_ptr->getProteoformPtr()->getSpeciesId() << "\t"
      << prsm_ptr->getProteoformPtr()->getSeqName() << " "
      << prsm_ptr->getProteoformPtr()->getSeqDesc() << "\t"
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

  ProteoformPtr proteo_ptr = prsm_ptr->getProteoformPtr();
  ExtendMsPtrVec refine_ms_ptr_vec = prsm_ptr->getRefineMsPtrVec();
  int proteo_len = proteo_ptr->getEndPos() - proteo_ptr->getStartPos() + 1;
  int third_start = 1;
  int third_first = proteo_len /3 + 1;
  int third_second = proteo_len * 2 / 3 + 1;
  int third_end = proteo_len;

  file << (proteo_len - 1) << "\t";
  file << (third_first - third_start) << "\t";
  file << (third_second - third_first) << "\t";
  file << (third_end - third_second) << "\t";

  std::vector<bool> comb_n_ion (proteo_len + 1, false);
  std::vector<bool> comb_c_ion (proteo_len + 1, false);
  std::vector<bool> comb_both_ion (proteo_len + 1, false);

  std::vector<std::vector<bool>> n_ion_2d;
  std::vector<std::vector<bool>> c_ion_2d;
  std::vector<std::vector<bool>> both_ion_2d;
  for (size_t s = 0; s < deconv_ms_ptr_vec.size(); s++) {
    //get ion_pair
    PeakIonPairPtrVec pair_ptrs = PeakIonPairFactory::genePeakIonPairs(prsm_ptr->getProteoformPtr(), 
                                                                       refine_ms_ptr_vec[s],
                                                                       min_mass_);
    std::vector<bool> n_ion (proteo_len + 1, false);
    std::vector<bool> c_ion (proteo_len + 1, false);
    std::vector<bool> both_ion (proteo_len + 1, false);
    for(size_t p = 0; p<pair_ptrs.size(); p++){
      int pos = pair_ptrs[p]->getTheoPeakPtr()->getIonPtr()->getPos();
      //LOG_DEBUG("start pos " << prot_ptr->getStartPos() << " pos " << pos);
      if(pair_ptrs[p]->getTheoPeakPtr()->getIonPtr()->getIonTypePtr()->isNTerm()){
        n_ion[pos] = true;
        both_ion[pos] = true;
        comb_n_ion[pos] = true;
        comb_both_ion[pos] = true;
      }
      else{
        c_ion[pos] = true;
        both_ion[pos] = true;
        comb_c_ion[pos] = true;
        comb_both_ion[pos] = true;
      }
    }
    n_ion_2d.push_back(n_ion);
    c_ion_2d.push_back(c_ion);
    both_ion_2d.push_back(both_ion);
  }

  for (size_t s = 0; s < deconv_ms_ptr_vec.size(); s++) {
    file << countCoverage(n_ion_2d[s], third_start, third_end) << "\t";
    file << countCoverage(n_ion_2d[s], third_start, third_first) << "\t";
    file << countCoverage(n_ion_2d[s], third_first, third_second) << "\t";
    file << countCoverage(n_ion_2d[s], third_second, third_end) << "\t";

    file << countCoverage(c_ion_2d[s], third_start, third_end) << "\t";
    file << countCoverage(c_ion_2d[s], third_start, third_first) << "\t";
    file << countCoverage(c_ion_2d[s], third_first, third_second) << "\t";
    file << countCoverage(c_ion_2d[s], third_second, third_end) << "\t";

    file << countCoverage(both_ion_2d[s], third_start, third_end) << "\t";
    file << countCoverage(both_ion_2d[s], third_start, third_first) << "\t";
    file << countCoverage(both_ion_2d[s], third_first, third_second) << "\t";
    file << countCoverage(both_ion_2d[s], third_second, third_end) << "\t";
  }

  file << countCoverage(comb_n_ion, third_start, third_end) << "\t";
  file << countCoverage(comb_n_ion, third_start, third_first) << "\t";
  file << countCoverage(comb_n_ion, third_first, third_second) << "\t";
  file << countCoverage(comb_n_ion, third_second, third_end) << "\t";

  file << countCoverage(comb_c_ion, third_start, third_end) << "\t";
  file << countCoverage(comb_c_ion, third_start, third_first) << "\t";
  file << countCoverage(comb_c_ion, third_first, third_second) << "\t";
  file << countCoverage(comb_c_ion, third_second, third_end) << "\t";

  file << countCoverage(comb_both_ion, third_start, third_end) << "\t";
  file << countCoverage(comb_both_ion, third_start, third_first) << "\t";
  file << countCoverage(comb_both_ion, third_first, third_second) << "\t";
  file << countCoverage(comb_both_ion, third_second, third_end) << "\t";

  ResSeqPtr res_seq_ptr = proteo_ptr->getResSeqPtr();

  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    AcidPtr acid = acid_ptr_vec_[i];
    file << countAcid(res_seq_ptr, acid) << "\t";
    for (size_t j = 0; j < deconv_ms_ptr_vec.size(); j++) {
      file << countAcidLeftCoverage(res_seq_ptr, acid, n_ion_2d[j]) << "\t";
      file << countAcidLeftCoverage(res_seq_ptr, acid, c_ion_2d[j]) << "\t";
      file << countAcidLeftCoverage(res_seq_ptr, acid, both_ion_2d[j]) << "\t";

      file << countAcidRightCoverage(res_seq_ptr, acid, n_ion_2d[j]) << "\t";
      file << countAcidRightCoverage(res_seq_ptr, acid, c_ion_2d[j]) << "\t";
      file << countAcidRightCoverage(res_seq_ptr, acid, both_ion_2d[j]) << "\t";
    }

    file << countAcidLeftCoverage(res_seq_ptr, acid, comb_n_ion) << "\t";
    file << countAcidLeftCoverage(res_seq_ptr, acid, comb_c_ion) << "\t";
    file << countAcidLeftCoverage(res_seq_ptr, acid, comb_both_ion) << "\t";

    file << countAcidRightCoverage(res_seq_ptr, acid, comb_n_ion) << "\t";
    file << countAcidRightCoverage(res_seq_ptr, acid, comb_c_ion) << "\t";
    file << countAcidRightCoverage(res_seq_ptr, acid, comb_both_ion) << "\t";
  }

  file << std::endl;
}

void PrsmStat::process() {
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName(); 
  std::string base_name = FileUtil::basename(spectrum_file_name);
  std::string output_file_name = base_name + "." + output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  //write title
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
      << "FDR" << "\t";
  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();

  file << "Total break point number " << "\t"
      << "Left break point number " <<  "\t"
      << "Middle break point number " <<  "\t"
      << "Right break point number " <<  "\t";
  for (int i = 1; i <= group_spec_num; i++) {
    file << "Spectrum " << i << " N_term_ion_coverage" << "\t"
        << "Spectrum " << i << " left_N_term_ion_coverage" << "\t"
        << "Spectrum " << i << " middle_N_term_ion_coverage" << "\t"
        << "Spectrum " << i << " right_N_term_ion_coverage" << "\t"

        << "Spectrum " << i << " C_term_ion_coverage" << "\t"
        << "Spectrum " << i << " left_C_term_ion_coverage" << "\t"
        << "Spectrum " << i << " middle_C_term_ion_coverage" << "\t"
        << "Spectrum " << i << " right_C_term_ion_coverage" << "\t"

        << "Spectrum " << i << " Both_term_ion_coverage" << "\t"
        << "Spectrum " << i << " left_Both_term_ion_coverage" << "\t"
        << "Spectrum " << i << " middle_both_term_ion_coverage" << "\t"
        << "Spectrum " << i << " right_Both_term_ion_coverage" << "\t";
  }

  file << "Combined spectra N_term_ion_coverage" << "\t"
      << "Combined spectra left_N_term_ion_coverage" << "\t"
      << "Combined spectra middle_N_term_ion_coverage" << "\t"
      << "Combined spectra right_N_term_ion_coverage" << "\t"

      << "Combined spectra C_term_ion_coverage" << "\t"
      << "Combined spectra left_C_term_ion_coverage" << "\t"
      << "Combined spectra middle_C_term_ion_coverage" << "\t"
      << "Combined spectra right_C_term_ion_coverage" << "\t"

      << "Combined spectra both_term_ion_coverage" << "\t"
      << "Combined spectra left_both_term_ion_coverage" << "\t"
      << "Combined spectra middle_both_term_ion_coverage" << "\t"
      << "Combined spectra right_both_term_ion_coverage" << "\t";

  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string acid = acid_ptr_vec_[i]->getOneLetter();
    file << acid << " number" << "\t";
    for (int j = 0; j < group_spec_num; j++) {
      file << "Spectrum " << j << " " << acid << " left N term" << "\t";
      file << "Spectrum " << j << " " << acid << " left C term" << "\t";
      file << "Spectrum " << j << " " << acid << " left both term" << "\t";
      file << "Spectrum " << j << " " << acid << " right N term" << "\t";
      file << "Spectrum " << j << " " << acid << " right C term" << "\t";
      file << "Spectrum " << j << " " << acid << " right both term" << "\t";
    }
    file << "Combined spectra " << acid << " left N term" << "\t";
    file << "Combined spectra " << acid << " left C term" << "\t";
    file << "Combined spectra " << acid << " left both term" << "\t";
    file << "Combined spectra " << acid << " right N term" << "\t";
    file << "Combined spectra " << acid << " right C term" << "\t";
    file << "Combined spectra " << acid << " right both term" << "\t";
  }

  file << std::endl;

  std::string input_file_name = FileUtil::basename(spectrum_file_name) + "." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileName();
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr_->getFixModPtrVec();

  std::string sp_file_name = prsm_para_ptr_->getSpectrumFileName();
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  // init variables
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          prsm_para_ptr_->getSpParaPtr()->getActivationPtr(),
                          prsm_para_ptr_->getSpParaPtr()->getSkipList());

  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  std::vector<SpectrumSetPtr> spec_set_vec = sp_reader.getNextSpectrumSet(sp_para_ptr);

  while (spec_set_vec[0] != nullptr) {
    if(spec_set_vec[0]->isValid()){
      int spec_id = spec_set_vec[0]->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_vec[0]->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec 
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        writePrsm(file, prsm_ptr);
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
      }
    }
    spec_set_vec = sp_reader.getNextSpectrumSet(sp_para_ptr);
  }

  sp_reader.close();
  prsm_reader.close();
  //write end;
  file.close();
}

} /* namespace prot */
