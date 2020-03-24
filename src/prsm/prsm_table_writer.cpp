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
#include "prsm/prsm_table_writer.hpp"

namespace toppic {

PrsmTableWriter::PrsmTableWriter(PrsmParaPtr prsm_para_ptr, 
                                 std::string argu_str,
                                 const std::string &input_file_ext, 
                                 const std::string &output_file_ext):
    prsm_para_ptr_(prsm_para_ptr),
    input_file_ext_(input_file_ext),
    argu_str_(argu_str),
    output_file_ext_(output_file_ext) {}

void PrsmTableWriter::write() {
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName();
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string output_file_name = base_name + output_file_ext_;
  std::ofstream file;
  file.open(output_file_name.c_str());
  file << argu_str_;
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
      << "Protein accession" << delim
      << "Protein description" << delim
      << "First residue" << delim
      << "Last residue" << delim
      << "Proteoform" << delim
      << "#unexpected modifications" << delim
      << "MIScore" << delim
      << "#variable PTMs" << delim
      << "#matched peaks" << delim
      << "#matched fragment ions" << delim
      << "P-value" << delim
      << "E-value" << delim
      << "Q-value (spectral FDR)" << delim
      << "Proteoform FDR" << std::endl;

  std::string input_file_name = file_util::basename(spectrum_file_name) + "." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileName();
  FastaIndexReaderPtr seq_reader = std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr_->getFixModPtrVec();
  PrsmReader prsm_reader(input_file_name);
  // LOG_DEBUG("start read prsm");
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
  // LOG_DEBUG("end read prsm");

  // init variables
  std::string sp_file_name = prsm_para_ptr_->getSpectrumFileName();
  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  SimpleMsAlignReaderPtr ms_reader_ptr = std::make_shared<SimpleMsAlignReader>(sp_file_name, 
                                                                               group_spec_num,
                                                                               sp_para_ptr->getActivationPtr());
  SpectrumSetPtr spec_set_ptr;
  while ((spec_set_ptr = spectrum_set_factory::readNextSpectrumSetPtr(ms_reader_ptr, sp_para_ptr)) != nullptr) {
    if (spec_set_ptr->isValid()) {
      int spec_id = spec_set_ptr->getSpectrumId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec
            = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        writePrsm(file, prsm_ptr);
        // LOG_DEBUG("start read prsm");
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
        // LOG_DEBUG("end read prsm");
      }
    }
  }
  prsm_reader.close();
  // write end;
  file.close();
}

void PrsmTableWriter::writePrsm(std::ofstream &file, PrsmPtr prsm_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  std::string retention_time;
  std::string delim = "\t";

  int ptm_num = prsm_ptr->getProteoformPtr()->getMassShiftNum(AlterType::UNEXPECTED);
  int peak_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = prsm_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + str_util::toString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
    retention_time = retention_time + str_util::fixedToString(deconv_ms_ptr_vec[i]->getMsHeaderPtr()->getRetentionTime(), 2) + " ";
  }

  str_util::trim(spec_ids);
  str_util::trim(spec_activations);
  str_util::trim(spec_scans);
  str_util::trim(retention_time);

  if (deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getRetentionTime() <= 0.0) retention_time = "-";

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
    file << "-" << delim;
  }

  file << prsm_ptr->getFracFeatureScore() << delim;

  file << prsm_ptr->getProteoformPtr()->getSeqName() << delim
      << "\"" << prsm_ptr->getProteoformPtr()->getSeqDesc() << "\"" << delim
      << (prsm_ptr->getProteoformPtr()->getStartPos() + 1) << delim
      << (prsm_ptr->getProteoformPtr()->getEndPos() + 1) << delim
      << "\"" << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\"" << delim
      << ptm_num << delim
      << prsm_ptr->getProteoformPtr()->getMIScore() << delim
      << prsm_ptr->getProteoformPtr()->getVariablePtmNum() << delim
      << prsm_ptr->getMatchPeakNum() << delim
      << prsm_ptr->getMatchFragNum() << delim
      << prsm_ptr->getPValue() << delim
      << prsm_ptr->getEValue() << delim;

  double fdr = prsm_ptr->getFdr();
  if (fdr >= 0) {
    file << fdr << delim;
  } else {
    file << "-" << delim;
  }

  double proteoform_fdr = prsm_ptr->getProteoformFdr();
  if (proteoform_fdr >= 0) {
    file << proteoform_fdr << std::endl;
  } else {
    file << "-" << std::endl;
  }
}

}  // namespace toppic
