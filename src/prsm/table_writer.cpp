#include <iomanip>
#include <boost/algorithm/string.hpp>

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

  int group_spec_num = prsm_para_ptr_->getGroupSpecNum();
  SpParaPtr sp_para_ptr = prsm_para_ptr_->getSpParaPtr();
  MsAlignReader reader (spectrum_file_name, group_spec_num);

  std::string input_file_name 
      = basename(spectrum_file_name) + "." + input_file_ext_;
  PrsmReader prsm_reader(input_file_name);
  ResiduePtrVec residue_ptr_vec = prsm_para_ptr_->getFixModResiduePtrVec();
  faidx_t *fai = fai_load(prsm_para_ptr_->getSearchDbFileName().c_str());
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(fai, residue_ptr_vec);

  SpectrumSetPtr spec_set_ptr;
  while ((spec_set_ptr= reader.getNextSpectrumSet(sp_para_ptr)) != nullptr) {
    PrsmPtrVec selected_prsm_ptrs;
    //LOG_DEBUG("spectrum id " << ms_ptr->getHeaderPtr()->getId() << " prsm ptr " << prsm_ptr);
    while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_set_ptr->getSpecId()) {
      selected_prsm_ptrs.push_back(prsm_ptr);
      prsm_ptr = prsm_reader.readOnePrsm(fai, residue_ptr_vec);
    }
    writeSelectedPrsms(file, selected_prsm_ptrs, spec_set_ptr);
  }
  fai_destroy(fai);
  reader.close();
  prsm_reader.close();
  //write end;
  file.close();
}

void TableWriter::writeSelectedPrsms(std::ofstream &file, PrsmPtrVec &prsm_ptrs, 
                                     SpectrumSetPtr spec_set_ptr) {
  std::string spec_ids;
  std::string spec_activations;
  std::string spec_scans;
  int peak_num = 0;
  DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
  for (size_t i = 0; i < deconv_ms_ptr_vec.size(); i++) {
    spec_ids = spec_ids + std::to_string(deconv_ms_ptr_vec[i]->getHeaderPtr()->getId()) + " ";
    spec_activations = spec_activations + deconv_ms_ptr_vec[i]->getHeaderPtr()->getActivationPtr()->getName() + " ";
    spec_scans = spec_scans + deconv_ms_ptr_vec[i]->getHeaderPtr()->getScansString() + " ";
    peak_num += deconv_ms_ptr_vec[i]->size();
  }
  boost::algorithm::trim(spec_ids);
  boost::algorithm::trim(spec_activations);
  boost::algorithm::trim(spec_scans);

  file << std::setprecision(10);
  for(size_t i=0;i<prsm_ptrs.size();i++){
    //LOG_DEBUG("prec mass " << prsm_ptrs[i]->getOriPrecMass());
    file << prsm_para_ptr_->getSpectrumFileName() << "\t"
        << prsm_ptrs[i]->getId() << "\t"
        << spec_ids << "\t"
        << spec_activations<< "\t"
        << spec_scans << "\t"
        << peak_num << "\t"
        << deconv_ms_ptr_vec[0]->getHeaderPtr()->getPrecCharge() << "\t"
        << prsm_ptrs[i]->getOriPrecMass()<< "\t"//"Precursor_mass"
        << prsm_ptrs[i]->getAdjustedPrecMass() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getId() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getSpeciesId() << "\t"
        << prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getName() << " "
        << prsm_ptrs[i]->getProteoformPtr()->getDbResSeqPtr()->getDesc() << "\t"
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
