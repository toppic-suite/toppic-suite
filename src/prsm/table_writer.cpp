/*
 * table_writer.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: xunlikun
 */

#include "base/file_util.hpp"
#include "prsm/table_writer.hpp"

namespace prot {

TableWriter::TableWriter(PrsmParaPtr prsm_para_ptr, 
                         std::string input_file_ext, std::string output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
}


void TableWriter::write(){
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName(); 
  std::string base_name = basename(spectrum_file_name);
  std::string output_file_name = base_name + "." + output_file_ext_;
  std::ofstream file_; 
  file_.open(output_file_name.c_str());
  //write title
  file_ << "Data_file_name" << "\t"
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
      << std::endl;

  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string input_file_name = base_name + "." + input_file_ext_;
  LOG_DEBUG("input file_name " << input_file_name);
  PrsmPtrVec prsms = readPrsm(input_file_name, raw_forms);
  sort(prsms.begin(),prsms.end(),prsmSpectrumIdUpMatchFragUp);
  LOG_DEBUG("read prsm complete ");
  addSpectrumPtrsToPrsms(prsms, prsm_para_ptr_);
  LOG_DEBUG("prsms loaded");

  for(unsigned int i=0;i<prsms.size();i++){
    file_ << spectrum_file_name << "\t"
        << prsms[i]->getId() << "\t"
        << prsms[i]->getSpectrumId()<< "\t"
        << prsms[i]->getDeconvMsPtr()->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
        << prsms[i]->getSpectrumScan() << "\t"
        << prsms[i]->getDeconvMsPtr()->size()<< "\t"
        << prsms[i]->getDeconvMsPtr()->getHeaderPtr()->getPrecCharge() << "\t"
        << prsms[i]->getOriPrecMass()<< "\t"//"Precursor_mass"
        << prsms[i]->getAdjustedPrecMass() << "\t"
        << prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
        << prsms[i]->getProteoformPtr()->getSpeciesId() << "\t"
        << prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getName() << "\t"
        << prsms[i]->getProteoformPtr()->getDbResSeqPtr()->getSeqMass() << "\t"
        << prsms[i]->getProteoformPtr()->getStartPos() << "\t"
        << prsms[i]->getProteoformPtr()->getEndPos() << "\t"
        << prsms[i]->getProteoformPtr()->getProteinMatchSeq() << "\t"
        << prsms[i]->getProteoformPtr()->getUnexpectedChangeNum() << "\t"
        << prsms[i]->getMatchPeakNum() << "\t"
        << prsms[i]->getMatchFragNum() << "\t"
        << prsms[i]->getPValue() << "\t"
        << prsms[i]->getEValue() << "\t"
        << prsms[i]->getProbPtr()->getOneProtProb()<< "\t"
        << prsms[i]->getFdr() << "\t"
        << std::endl;
//    std::cout<<prsms[i]->getSpectrumId()<<std::endl;
//    std::cout<<prsms[i]->getProteoformPtr()->getProteinMatchSeq()<<std::endl;
  }
  //write end;
  file_.close();
}

} /* namespace prot */
