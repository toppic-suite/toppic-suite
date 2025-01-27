//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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
#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_match_table_writer.hpp"

namespace toppic {

PrsmMatchTableWriter::PrsmMatchTableWriter(PrsmParaPtr prsm_para_ptr, 
                                 std::string argu_str,
                                 const std::string &input_file_ext):
    prsm_para_ptr_(prsm_para_ptr),
    input_file_ext_(input_file_ext),
    argu_str_(argu_str) {
      std::string db_file_name = prsm_para_ptr_->getSearchDbFileNameWithFolder();
      search_match_ptr_ = std::make_shared<SearchFastaMatch>(db_file_name);
    }

void PrsmMatchTableWriter::writeHeader(std::ofstream &file) {
  // write title
  std::string delim = "\t";
  file << "Data file name" << delim
      << "Prsm ID" << delim
      << "Spectrum ID"<< delim
      << "Fragmentation" << delim
      << "Scan(s)" << delim
      << "Retention time" << delim
      << "# Masses"<< delim
      << "Charge" << delim
      << "Precursor mass" << delim
      << "Adjusted precursor mass" << delim
      << "Proteoform ID" << delim
      << "Proteoform intensity" << delim
      << "Feature ID" << delim
      << "Feature intensity" << delim
      << "Feature score" << delim
      << "Feature apex time" << delim
      << "# Protein hits" << delim
      << "Protein accession" << delim
      << "Protein description" << delim
      << "First residue" << delim
      << "Last residue" << delim
      << "Special amino acids" << delim
      << "Database protein sequence" << delim
      << "Previous amino acid" << delim
      << "Proteoform" << delim
      << "Next amino acid" << delim
      << "Proteoform mass" << delim
      << "Protein N-terminal form" << delim
      << "Fixed PTMs" << delim
      << "# Unexpected modifications" << delim
      << "Unexpected modifications" << delim
      << "# Variable PTMs" << delim
      << "Variable PTMs" << delim
      << "MIScore" << delim
      << "# Matched masses" << delim
      << "# Matched fragments" << delim
      << "E-value" << delim
      << "Spectrum-level Q-value" << delim
      << "Proteoform-level Q-value" << std::endl;
}

void PrsmMatchTableWriter::write(const std::string &output_file_ext, 
                                 bool write_multiple_matches) {
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName();
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string output_file_name = base_name + output_file_ext;
  std::ofstream file;
  file.open(output_file_name.c_str());
  file << argu_str_ << std::endl;

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
  MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(sp_file_name, 
                                                                   group_spec_num,
                                                                   sp_para_ptr->getActivationPtr());
  DeconvMsPtrVec deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  PrsmPtrVec prsm_list;
  int max_proteo_id = -1;
  while (deconv_ms_ptr_vec.size() != 0) {
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
    if (header_ptr->containsPrec()) {
      double prec_mono_mass = header_ptr->getFirstPrecMonoMass() - sp_para_ptr->getNTermLabelMass();
      SpectrumSetPtr spec_set_ptr  
        = spectrum_set_factory::geneSpectrumSetPtr(deconv_ms_ptr_vec,
                                                   sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        int spec_id = spec_set_ptr->getSpectrumId();
        while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
          DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
          prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
          double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
          ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
          prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
          prsm_list.push_back(prsm_ptr);
          if (prsm_ptr->getProteoformPtr()->getProteoClusterId() > max_proteo_id) {
            max_proteo_id = prsm_ptr->getProteoformPtr()->getProteoClusterId();
          }
          prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
        }
      }
    }
    deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  }
  prsm_reader.close();

  file << "Number of identified PrSMs: " << prsm_list.size() << std::endl;
  PrsmPtrVec2D prsm_2d; 
  for (size_t i = 0; i < max_proteo_id + 1; i++) {
    PrsmPtrVec prsm_list; 
    prsm_2d.push_back(prsm_list);
  }
  for (size_t i = 0; i < prsm_list.size(); i++) {
    int proteo_id = prsm_list[i]->getProteoformPtr()->getProteoClusterId();
    prsm_2d[proteo_id].push_back(prsm_list[i]);
  }
  int proteo_num = 0;
  std::set<std::string> prot_id_set;
  for (size_t i = 0; i < max_proteo_id + 1; i++) {
    if (prsm_2d[i].size() > 0) {
      proteo_num++;
      std::sort(prsm_2d[i].begin(), prsm_2d[i].end(), Prsm::cmpEValueIncProtInc);
      prot_id_set.insert(prsm_2d[i][0]->getProteoformPtr()->getSeqName());
    }
  }
  file << "Number of identified proteoforms: " << proteo_num << std::endl;
  file << "Number of identified proteins: " << prot_id_set.size() << std::endl << std::endl;
  writeHeader(file);
  for (size_t i = 0; i < prsm_list.size(); i++) {
    writePrsmStandardFormat(file, prsm_list[i], write_multiple_matches);
  }
  // write end;
  file.close();
}

void PrsmMatchTableWriter::writePrsmStandardFormat(std::ofstream &output_file, PrsmPtr prsm_ptr, 
                                                   bool write_multiple_matches) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  std::string retention_time;
  std::string delim = "\t";
  std::string empty_str = "-";
  std::vector<std::pair<FastaSeqPtr,int>> matches;
  
  
  int mass_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getSpecId()) + " ";
    spec_activations = spec_activations 
        + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    mass_num += deconv_ms_ptr_vec[i]->size();
    retention_time = retention_time 
        + str_util::fixedToString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getRetentionTime()/60, 3) + " ";
  }

  str_util::trim(spec_ids);
  str_util::trim(spec_activations);
  str_util::trim(spec_scans);
  str_util::trim(retention_time);

  if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getRetentionTime() <= 0.0) {
    retention_time = empty_str;
  }

  matches = search_match_ptr_->process(prsm_ptr);

  std::stringstream line_str;
  line_str << std::setprecision(10);
  LOG_DEBUG("start output prsm ");
  line_str << prsm_ptr->getFileName() << delim
    << prsm_ptr->getPrsmId() << delim
    << spec_ids << delim
    << spec_activations << delim
    << spec_scans << delim
    << retention_time << delim
    << mass_num << delim
    << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecCharge() << delim
    << prsm_ptr->getOriPrecMass()<< delim
    << prsm_ptr->getAdjustedPrecMass() << delim
    << prsm_ptr->getProteoformPtr()->getProteoClusterId() << delim
    << prsm_ptr->getProteoformPtr()->getProteoInte() << delim;

  if (prsm_ptr->getFracFeatureInte() > 0) {
    line_str << prsm_ptr->getFracFeatureId() << delim;
    std::ostringstream str_stream;
    str_stream << std::scientific << std::setprecision(3);
    str_stream << prsm_ptr->getFracFeatureInte();
    line_str << str_stream.str() << delim;
    line_str << prsm_ptr->getFracFeatureScore() << delim;
    line_str << str_util::fixedToString(prsm_ptr->getFracFeatureApexTime()/60,3) << delim;
  } else {
    line_str << empty_str << delim;
    line_str << empty_str << delim;
    line_str << empty_str << delim;
    line_str << empty_str << delim;
  }

  line_str << matches.size() << delim;
  
  ProteoformPtr form_ptr = prsm_ptr->getProteoformPtr();

  std::string seq = form_ptr->getProteoformMatchSeq();

  int start_pos = form_ptr->getStartPos();
  int end_pos = form_ptr->getEndPos();
  line_str << form_ptr->getSeqName() << delim
      << form_ptr->getSeqDesc() << delim
      << (start_pos + 1) << delim
      << (end_pos + 1) << delim
      << form_ptr->getFastaSeqPtr()->getAcidReplaceStr(start_pos, end_pos) << delim
      << form_ptr->getFastaSeqPtr()->getSubSeq(start_pos, end_pos) << delim
      << form_ptr->getPrevAminoAcid() << delim
      << form_ptr->getProteoformMatchSeq() << delim
      << form_ptr->getNextAminoAcid() << delim
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
    line_str << fdr << delim;
  } else {
    line_str << empty_str << delim;
  }

  double proteoform_fdr = prsm_ptr->getProteoformFdr();
  if (proteoform_fdr >= 0) {
    line_str << proteoform_fdr << std::endl;
  } else {
    line_str << empty_str << std::endl;
  }

  if (write_multiple_matches) {
    // print out other matches
    for (size_t i = 0; i < matches.size(); i++) {
      FastaSeqPtr seq_ptr = matches[i].first;
      int seq_pos = matches[i].second;
      if (seq_ptr->getName() == form_ptr->getSeqName()) {
        continue;
      }

    line_str << prsm_ptr->getFileName() << delim
        << prsm_ptr->getPrsmId() << delim
        << delim
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
      line_str << delim
        << delim
        << delim
        << delim
        << delim;

      line_str << seq_ptr->getName() << delim
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
        << delim
        << delim
        << delim;

      // fdr
      line_str << delim
        << std::endl;
    }
  }

  output_file << line_str.str();
}

}  // namespace toppic
