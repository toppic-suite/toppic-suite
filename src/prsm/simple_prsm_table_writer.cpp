/*
 * simple_table_writer.cpp
 *
 *  Created on: Aug 11, 2014
 *      Author: qkou
 */

#include <iomanip>

#include "base/file_util.hpp"
#include "prsm/simple_prsm_table_writer.hpp"

namespace prot {

SimplePrsmTableWriter::SimplePrsmTableWriter(PrsmParaPtr prsm_para_ptr, 
                                             const std::string &input_file_ext, 
                                             const std::string &output_file_ext) {
  prsm_para_ptr_ = prsm_para_ptr;
  input_file_ext_ = input_file_ext;
  output_file_ext_ = output_file_ext;
}

void SimplePrsmTableWriter::write() {
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName(); 
  std::string base_name = basename(spectrum_file_name);
  std::string output_file_name = base_name + "." + output_file_ext_;
  std::ofstream file_; 
  file_.open(output_file_name.c_str());
  //write title
  file_ << "Spectrum_ID" << "\t"
        << "Scan(s)" << "\t"
        << "Precursor_ID" << "\t"
        << "Precursor_mass" << "\t"
        << "Protein_ID" << "\t"
        << "Species_ID" << "\t"
        << "Score" << "\t"
        << "Protein_name"
        << std::endl;

  std::string input_file_name = base_name + "." + input_file_ext_;
  LOG_DEBUG("input file_name " << input_file_name);
  SimplePrsmPtrVec simple_prsm_ptrs = readSimplePrsms(input_file_name); 
  LOG_DEBUG("read simple prsm complete ");
  file_ << std::setprecision(10);
  for(size_t i = 0; i < simple_prsm_ptrs.size();i++){
    file_ << simple_prsm_ptrs[i]->getSpectrumId() << "\t"
          << simple_prsm_ptrs[i]->getSpectrumScan() << "\t"
          << simple_prsm_ptrs[i]->getPrecursorId() << "\t"
          << simple_prsm_ptrs[i]->getPrecMass() << "\t"
          << simple_prsm_ptrs[i]->getSeqId() << "\t"
          << simple_prsm_ptrs[i]->getScore() << "\t"
          << simple_prsm_ptrs[i]->getSeqName()
          << std::endl;
  }

  //write end;
  file_.close();
}

}
