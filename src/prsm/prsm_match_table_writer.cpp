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

#include <iomanip>
#include <sstream>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_match_table_writer.hpp"

namespace toppic {

PrsmMatchTableWriter::PrsmMatchTableWriter(PrsmParaPtr prsm_para_ptr, 
                                 std::string argu_str,
                                 const std::string &input_file_ext, 
                                 const std::string &output_file_ext, 
                                 bool write_multiple_matches):
    prsm_para_ptr_(prsm_para_ptr),
    input_file_ext_(input_file_ext),
    argu_str_(argu_str),
    output_file_ext_(output_file_ext),
    write_multiple_matches_(write_multiple_matches) {
      std::string db_file_name = prsm_para_ptr_->getSearchDbFileNameWithFolder();
      search_match_ptr_ = std::make_shared<SearchFastaMatch>(db_file_name);
    }

void PrsmMatchTableWriter::write() {
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName();
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string output_file_name = base_name + output_file_ext_;
  std::ofstream file;
  file.open(output_file_name.c_str());
  file << argu_str_;
  file << "\n";
  // write title
  std::string delim = "\t";
  
  file << "Data file name" << delim
      << "Prsm ID" << delim
      << "Spectrum ID"<< delim
      << "Fragmentation" << delim
      << "Scan(s)" << delim
      << "Retention time" << delim
      << "#peaks"<< delim
      << "Charge" << delim
      << "Precursor mass" << delim
      << "Adjusted precursor mass" << delim
      << "Proteoform ID" << delim
      << "Feature intensity" << delim
      << "Feature score" << delim
      << "Feature apex time" << delim
      << "#Protein hits" << delim
      << "Protein accession" << delim
      << "Protein description" << delim
      << "First residue" << delim
      << "Last residue" << delim
      << "Special amino acids" << delim
      << "Database protein sequence" << delim
      << "Proteoform" << delim
      << "Proteoform mass" << delim
      << "Protein N-terminal form" << delim
      << "Fixed PTMs" << delim
      << "#unexpected modifications" << delim
      << "unexpected modifications" << delim
      << "#variable PTMs" << delim
      << "variable PTMs" << delim
      << "MIScore" << delim
      << "#matched peaks" << delim
      << "#matched fragment ions" << delim
      << "E-value" << delim
      << "Spectrum-level Q-value" << delim
      << "Proteoform-level Q-value" << std::endl;

  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;  
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileNameWithFolder();
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr_->getFixModPtrVec();
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  // init variables
  std::string sp_file_name = prsm_para_ptr_->getSpectrumFileName();
  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  SimpleMsAlignReaderPtr ms_reader_ptr 
      = std::make_shared<SimpleMsAlignReader>(sp_file_name, 
                                              group_spec_num,
                                              sp_para_ptr->getActivationPtr());
  SpectrumSetPtr spec_set_ptr;
  while ((spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr)) 
         != nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        //writePrsm(file, prsm_ptr);
        writePrsmStandardFormat(file, prsm_ptr);
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
      }
    }
  }
  prsm_reader.close();
  // write end;
  file.close();
}


void PrsmMatchTableWriter::writePrsm(std::ofstream &file, PrsmPtr prsm_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  std::string retention_time;
  std::string delim = "\t";
  std::string empty_str = "-";
  std::vector<std::pair<FastaSeqPtr,int>> matches;

  int peak_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations 
        + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
    retention_time = retention_time 
        + str_util::fixedToString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getRetentionTime(), 2) + " ";
  }

  str_util::trim(spec_ids);
  str_util::trim(spec_activations);
  str_util::trim(spec_scans);
  str_util::trim(retention_time);

  if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getRetentionTime() <= 0.0) {
    retention_time = empty_str;
  }

  matches = search_match_ptr_->process(prsm_ptr);

  file << std::setprecision(10);
  LOG_DEBUG("start output prsm ");
  file << prsm_ptr->getFileName() << delim
      << prsm_ptr->getPrsmId() << delim
      << spec_ids << delim
      << spec_activations<< delim
      << spec_scans << delim
      << retention_time << delim
      << peak_num << delim
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << delim
      << prsm_ptr->getOriPrecMass()<< delim
      << prsm_ptr->getAdjustedPrecMass() << delim
      << prsm_ptr->getProteoformPtr()->getProteoClusterId() << delim;

  if (prsm_ptr->getSampleFeatureInte() > 0) {
    std::ostringstream str_stream;
    str_stream << std::scientific << std::setprecision(3);
    str_stream << prsm_ptr->getSampleFeatureInte();
    file << str_stream.str() << delim;
  } else {
    file << empty_str << delim;
  }

  file << prsm_ptr->getFracFeatureScore() << delim;
  file << prsm_ptr->getTimeApex() << delim;
  file << matches.size() << delim;
  
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();

  int start_pos = form_ptr->getStartPos();
  int end_pos = form_ptr->getEndPos();
  file << form_ptr->getSeqName() << delim
      << form_ptr->getSeqDesc() << delim
      << (start_pos + 1) << delim
      << (end_pos + 1) << delim
      << form_ptr->getFastaSeqPtr()->getAcidReplaceStr(start_pos, end_pos) << delim
      << form_ptr->getFastaSeqPtr()->getSubSeq(start_pos, end_pos) << delim
      << form_ptr->getProteoformMatchSeq() << delim
      << form_ptr->getMass() << delim
      << form_ptr->getProtModPtr()->getType() << delim
      << form_ptr->getAlterStr(AlterType::FIXED) << delim
      << form_ptr->getAlterNum(AlterType::UNEXPECTED) << delim
      << form_ptr->getAlterStr(AlterType::UNEXPECTED) << delim
      << form_ptr->getAlterNum(AlterType::VARIABLE) << delim
      << form_ptr->getAlterStr(AlterType::VARIABLE) << delim
      << form_ptr->getMIScore() << delim
      << prsm_ptr->getMatchPeakNum() << delim
      << prsm_ptr->getMatchFragNum() << delim
      << prsm_ptr->getEValue() << delim;

  double fdr = prsm_ptr->getFdr();
  if (fdr >= 0) {
    file << fdr << delim;
  } else {
    file << empty_str << delim;
  }

  double proteoform_fdr = prsm_ptr->getProteoformFdr();
  if (proteoform_fdr >= 0) {
    file << proteoform_fdr << std::endl;
  } else {
    file << empty_str << std::endl;
  }

  if (write_multiple_matches_) {
    // print out other matches
    for (size_t i = 0; i < matches.size(); i++) {
      FastaSeqPtr seq_ptr = matches[i].first;
      int seq_pos = matches[i].second;
      if (seq_ptr->getName() == form_ptr->getSeqName()) {
        continue;
      }

    file << prsm_ptr->getFileName() << delim
        << prsm_ptr->getPrsmId() << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim;

      // feature
      file << delim
        << delim
        << delim
        << delim;

      file << seq_ptr->getName() << delim
        << seq_ptr->getDesc() << delim
        << (seq_pos + 1) << delim
        << (seq_pos + form_ptr->getLen()) << delim
        << seq_ptr->getAcidReplaceStr(form_ptr->getStartPos(), form_ptr->getEndPos()) << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim;

      // fdr
      file << delim
        << std::endl;
    }
  }
}

void PrsmMatchTableWriter::writePrsmStandardFormat(std::ofstream &file, PrsmPtr prsm_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  std::string retention_time;
  std::string delim = "\t";
  std::string empty_str = "-";
  std::vector<std::pair<FastaSeqPtr,int>> matches;

  int peak_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations 
        + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
    retention_time = retention_time 
        + str_util::fixedToString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getRetentionTime(), 2) + " ";
  }

  str_util::trim(spec_ids);
  str_util::trim(spec_activations);
  str_util::trim(spec_scans);
  str_util::trim(retention_time);

  if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getRetentionTime() <= 0.0) {
    retention_time = empty_str;
  }

  matches = search_match_ptr_->process(prsm_ptr);

  file << std::setprecision(10);
  LOG_DEBUG("start output prsm ");
  file << prsm_ptr->getFileName() << delim
      << prsm_ptr->getPrsmId() << delim
      << spec_ids << delim
      << spec_activations<< delim
      << spec_scans << delim
      << retention_time << delim
      << peak_num << delim
      << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecCharge() << delim
      << prsm_ptr->getOriPrecMass()<< delim
      << prsm_ptr->getAdjustedPrecMass() << delim
      << prsm_ptr->getProteoformPtr()->getProteoClusterId() << delim;

  if (prsm_ptr->getSampleFeatureInte() > 0) {
    std::ostringstream str_stream;
    str_stream << std::scientific << std::setprecision(3);
    str_stream << prsm_ptr->getSampleFeatureInte();
    file << str_stream.str() << delim;
  } else {
    file << empty_str << delim;
  }

  file << prsm_ptr->getFracFeatureScore() << delim;
  file << prsm_ptr->getTimeApex() << delim;
  file << matches.size() << delim;
  
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();

  std::string seq = form_ptr->getProteoformMatchSeq();

  int start_pos = form_ptr->getStartPos();
  int end_pos = form_ptr->getEndPos();
  file << form_ptr->getSeqName() << delim
      << form_ptr->getSeqDesc() << delim
      << (start_pos + 1) << delim
      << (end_pos + 1) << delim
      << form_ptr->getFastaSeqPtr()->getAcidReplaceStr(start_pos, end_pos) << delim
      << form_ptr->getFastaSeqPtr()->getSubSeq(start_pos, end_pos) << delim
      << form_ptr->getProteoformMatchSeq() << delim
      << form_ptr->getMass() << delim
      << form_ptr->getProtModPtr()->getType() << delim
      << form_ptr->getAlterStr(AlterType::FIXED) << delim
      << form_ptr->getAlterNum(AlterType::UNEXPECTED) << delim
      << form_ptr->getAlterStr(AlterType::UNEXPECTED) << delim
      << form_ptr->getAlterNum(AlterType::VARIABLE) << delim
      << form_ptr->getAlterStr(AlterType::VARIABLE) << delim
      << form_ptr->getMIScore() << delim
      << prsm_ptr->getMatchPeakNum() << delim
      << prsm_ptr->getMatchFragNum() << delim
      << prsm_ptr->getEValue() << delim;

  double fdr = prsm_ptr->getFdr();
  if (fdr >= 0) {
    file << fdr << delim;
  } else {
    file << empty_str << delim;
  }

  double proteoform_fdr = prsm_ptr->getProteoformFdr();
  if (proteoform_fdr >= 0) {
    file << proteoform_fdr << std::endl;
  } else {
    file << empty_str << std::endl;
  }

  if (write_multiple_matches_) {
    // print out other matches
    for (size_t i = 0; i < matches.size(); i++) {
      FastaSeqPtr seq_ptr = matches[i].first;
      int seq_pos = matches[i].second;
      if (seq_ptr->getName() == form_ptr->getSeqName()) {
        continue;
      }

    file << prsm_ptr->getFileName() << delim
        << prsm_ptr->getPrsmId() << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim;

      // feature
      file << delim
        << delim
        << delim
        << delim;

      file << seq_ptr->getName() << delim
        << seq_ptr->getDesc() << delim
        << (seq_pos + 1) << delim
        << (seq_pos + form_ptr->getLen()) << delim
        << seq_ptr->getAcidReplaceStr(form_ptr->getStartPos(), form_ptr->getEndPos()) << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim
        << delim;

      // fdr
      file << delim
        << std::endl;
    }
  }
}

}  // namespace toppic
