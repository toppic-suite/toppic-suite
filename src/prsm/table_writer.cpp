
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_reader.hpp"
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

  MsAlignReader reader (spectrum_file_name);

  std::string input_file_name 
      = basename(spectrum_file_name + "." + input_file_ext_);
  PrsmReader prsm_reader(input_file_name);
  ResiduePtrVec residue_ptr_vec = prsm_para_ptr_->getFixModResiduePtrVec();
  faidx_t *fai = fai_load(prsm_para_ptr_->getSearchDbFileName().c_str());
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fai, residue_ptr_vec);

  DeconvMsPtr ms_ptr = reader.getNextMs();
  while (ms_ptr.get() != nullptr) {
    PrsmPtrVec selected_prsm_ptrs;
    while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == ms_ptr->getHeaderPtr()->getId()) {
      selected_prsm_ptrs.push_back(prsm_ptr);
      prsm_ptr = prsm_reader.readOnePrsm(fai, residue_ptr_vec);
    }
    writeSelectedPrsms(file, selected_prsm_ptrs, ms_ptr);
    ms_ptr = reader.getNextMs();
  }
  fai_destroy(fai);
  reader.close();
  prsm_reader.close();
  //write end;
  file.close();
}

void TableWriter::writeSelectedPrsms(std::ofstream &file, PrsmPtrVec &prsm_ptrs, 
                                     DeconvMsPtr ms_ptr) {
  for(size_t i=0;i<prsm_ptrs.size();i++){
    file << prsm_para_ptr_->getSpectrumFileName() << "\t"
        << prsm_ptrs[i]->getId() << "\t"
        << prsm_ptrs[i]->getSpectrumId()<< "\t"
        << ms_ptr->getHeaderPtr()->getActivationPtr()->getName()<< "\t"
        << prsm_ptrs[i]->getSpectrumScan() << "\t"
        << ms_ptr->size()<< "\t"
        << ms_ptr->getHeaderPtr()->getPrecCharge() << "\t"
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
        << prsm_ptrs[i]->getOneProtProb()<< "\t"
        << prsm_ptrs[i]->getFdr() << "\t"
        << std::endl;
  }
}

} /* namespace prot */
