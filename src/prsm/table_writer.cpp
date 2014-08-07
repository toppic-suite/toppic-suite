#include "base/file_util.hpp"
#include "prsm/table_writer.hpp"

namespace prot {

TableWriter::TableWriter(PrsmParaPtr prsm_para_ptr, 
                         const std::string &input_file_ext, 
                         const std::string &output_file_ext) {
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

  ProteoformPtrVec raw_form_ptrs 
      = readFastaToProteoform(prsm_para_ptr_->getSearchDbFileName(), 
                              prsm_para_ptr_->getFixModResiduePtrVec());

  LOG_DEBUG("protein data set loaded");
  std::string input_file_name = base_name + "." + input_file_ext_;
  LOG_DEBUG("input file_name " << input_file_name);
  PrsmPtrVec prsm_ptrs = readPrsm(input_file_name, raw_form_ptrs);
  sort(prsm_ptrs.begin(),prsm_ptrs.end(),prsmSpectrumIdUpPrecursorIdUp);
  LOG_DEBUG("read prsm complete ");
  addSpectrumPtrsToPrsms(prsm_ptrs, prsm_para_ptr_);
  LOG_DEBUG("prsm_ptrs loaded");

  for(size_t i=0;i<prsm_ptrs.size();i++){
    file_ << spectrum_file_name << "\t"
        << prsm_ptrs[i]->getId() << "\t"
        << prsm_ptrs[i]->getSpectrumId()<< "\t"
        << prsm_ptrs[i]->getDeconvMsPtr()->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
        << prsm_ptrs[i]->getSpectrumScan() << "\t"
        << prsm_ptrs[i]->getDeconvMsPtr()->size()<< "\t"
        << prsm_ptrs[i]->getDeconvMsPtr()->getHeaderPtr()->getPrecCharge() << "\t"
        << prsm_ptrs[i]->getOriPrecMass()<< "\t"//"Precursor_mass"
        << prsm_ptrs[i]->getAdjustedPrecMass() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getSpeciesId() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getName() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getSeqMass() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getStartPos() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getEndPos() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getProteinMatchSeq() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangeNum() << "\t"
        << prsm_ptrs[i]->getMatchPeakNum() << "\t"
        << prsm_ptrs[i]->getMatchFragNum() << "\t"
        << prsm_ptrs[i]->getPValue() << "\t"
        << prsm_ptrs[i]->getEValue() << "\t"
        << prsm_ptrs[i]->getProbPtr()->getOneProtProb()<< "\t"
        << prsm_ptrs[i]->getFdr() << "\t"
        << std::endl;
  }
  //write end;
  file_.close();
}

} /* namespace prot */
