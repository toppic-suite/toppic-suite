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


#include <iomanip>
#include <ctime>

#include <boost/algorithm/string.hpp>

#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/extend_ms_factory.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "console/argument.hpp"

namespace prot {

PrsmTableWriter::PrsmTableWriter(PrsmParaPtr prsm_para_ptr, 
                                 std::map<std::string, std::string> arguments,
                                 const std::string &input_file_ext, 
                                 const std::string &output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  arguments_ = arguments;
  output_file_ext_ = output_file_ext;
}


void PrsmTableWriter::write(){
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName(); 
  std::string base_name = FileUtil::basename(spectrum_file_name);
  std::string output_file_name = base_name + "." + output_file_ext_;
  std::ofstream file; 
  file.open(output_file_name.c_str());
  /*time_t ctt = time(0);*/
  //file << "Time: ";
  /*file << asctime(localtime(&ctt));*/
  Argument::outputArguments(file, arguments_);
  //write title
  file << "Data file name" << "\t"
      << "Prsm ID" << "\t"
      << "Spectrum ID"<< "\t"
      << "Fragmentation" << "\t"
      << "Scan(s)" << "\t"
      << "#peaks"<< "\t"
      << "Charge" << "\t"
      << "Precursor mass" << "\t"
      << "Adjusted precursor mass" << "\t"
      << "Proteoform ID" << "\t"
      << "Feature intensity" << "\t"
      << "Protein name" << "\t"
      << "First residue" << "\t"
      << "Last residue" << "\t"
      << "Proteoform" << "\t"
      << "#unexpected modifications" << "\t";

  if (prsm_para_ptr_->doLocaliztion()) {
    file << "MIScore" << "\t";
  }

  file << "#matched peaks" << "\t"
      << "#matched fragment ions" << "\t"
#ifdef TOPPIC
      << "P-value" << "\t"
      << "E-value" << "\t"
      //      << "One Protein probabilty"<< "\t"
      << "Q-value (spectral FDR)" << "\t"
      << "Proteoform FDR"
#endif
#ifdef MASS_GRAPH
      << "#Variable PTMs" << "\t"
#endif
      << std::endl;

  std::string input_file_name 
      = FileUtil::basename(spectrum_file_name) + "." + input_file_ext_;
  std::string db_file_name = prsm_para_ptr_->getSearchDbFileName();
  FastaIndexReaderPtr seq_reader(new FastaIndexReader(db_file_name));
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr_->getFixModPtrVec();
  PrsmReader prsm_reader(input_file_name);
  //LOG_DEBUG("start read prsm");
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
  //LOG_DEBUG("end read prsm");

  //init variables
  std::string sp_file_name = prsm_para_ptr_->getSpectrumFileName();
  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();
  MsAlignReader sp_reader(sp_file_name, group_spec_num,
                          prsm_para_ptr_->getSpParaPtr()->getActivationPtr());
  SpectrumSetPtr spec_set_ptr;
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  while((spec_set_ptr = sp_reader.getNextSpectrumSet(sp_para_ptr))!= nullptr){
    if(spec_set_ptr->isValid()){
      int spec_id = spec_set_ptr->getSpecId();
      while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
        DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
        prsm_ptr->setDeconvMsPtrVec(deconv_ms_ptr_vec);
        double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
        ExtendMsPtrVec extend_ms_ptr_vec 
            = ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, sp_para_ptr, new_prec_mass);
        prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
        writePrsm(file, prsm_ptr);
        //LOG_DEBUG("start read prsm");
        prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
        //LOG_DEBUG("end read prsm");
      }
    }
  }
  sp_reader.close();
  prsm_reader.close();
  //write end;
  file.close();
}

std::string outputChangePtr(ProteoformPtr proteoform_ptr) {
  StringPairVec string_pairs = proteoform_ptr->getFastaSeqPtr()->getAcidPtmPairVec();
  int start_pos = proteoform_ptr->getStartPos();
  ChangePtrVec change_vec = proteoform_ptr->getChangePtrVec(ChangeType::UNEXPECTED);
  std::string res = "";
  for (size_t i = 0; i < change_vec.size(); i++) {
    if (change_vec[i]->getLocalAnno() == nullptr) 
      continue;

    if (change_vec[i]->getLocalAnno()->getPtmPtr() != nullptr) {
      std::vector<double> scr_vec = change_vec[i]->getLocalAnno()->getScrVec();
      int left_db_bp = change_vec[i]->getLeftBpPos() + start_pos;
      int right_db_bp = change_vec[i]->getRightBpPos() + start_pos;
      res = res + change_vec[i]->getLocalAnno()->getPtmPtr()->getAbbrName() + "[";
      for (int j = left_db_bp; j < right_db_bp; j++) {
        //std::string acid_letter = fasta_seq.substr(j, 1);
        std::string acid_letter = string_pairs[j].first;
        double scr = std::floor(scr_vec[j - left_db_bp] * 1000) / 10;
        if (scr == 100) scr = 99.9;
        if (scr == 0) continue;

        res = res + acid_letter + std::to_string(j + 1) + ":";
        std::stringstream ss;
        ss << std::fixed << std::setprecision(1) << scr;
        res  = res + ss.str() + "%";
        if (j != right_db_bp - 1) {
          res = res + "; ";
        }
      }
      res = res + "]";
    }

    if (i != change_vec.size() - 1) {
      res = res + "; ";
    }
  }
  if (res == "") {
    res = "-";
  }
  return res;
}

void PrsmTableWriter::writePrsm(std::ofstream &file, PrsmPtr prsm_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  int ptm_num = prsm_ptr->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED);
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

  file << std::setprecision(10);
  LOG_DEBUG("start output prsm ");
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
      << prsm_ptr->getPrecFeatureInte() << "\t"
      << prsm_ptr->getProteoformPtr()->getSeqName() << " "
      << prsm_ptr->getProteoformPtr()->getSeqDesc() << "\t"
      << (prsm_ptr->getProteoformPtr()->getStartPos() + 1) << "\t"
      << (prsm_ptr->getProteoformPtr()->getEndPos() + 1) << "\t"
      << prsm_ptr->getProteoformPtr()->getProteinMatchSeq() << "\t"
      << ptm_num << "\t";

  if (prsm_para_ptr_->doLocaliztion()) {
    if (ptm_num == 0) {
      file << "-\t";
    } else {
      file << outputChangePtr(prsm_ptr->getProteoformPtr()) << "\t";
    }
  } 

  file << prsm_ptr->getMatchPeakNum() << "\t"
      << prsm_ptr->getMatchFragNum() << "\t";
#ifdef TOPPIC
  file << prsm_ptr->getPValue() << "\t"
      << prsm_ptr->getEValue() << "\t";
  //      << prsm_ptr->getOneProtProb()<< "\t"
  double fdr = prsm_ptr->getFdr();
  if (fdr >= 0) {
    file << fdr << "\t";
  } else { 
    file << "-" << "\t";
  }
  double proteoform_fdr = prsm_ptr->getProteoformFdr();
  if (proteoform_fdr >= 0) {
    file << proteoform_fdr << "\t";
  } else { 
    file << "-" << "\t";
  }
#endif
#ifdef MASS_GRAPH
  file << prsm_ptr->getProteoformPtr()->getVariablePtmNum() << "\t";
#endif
  file << std::endl;
}

} /* namespace prot */
