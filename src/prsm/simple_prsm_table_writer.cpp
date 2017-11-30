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


#include <iomanip>

#include "base/file_util.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_table_writer.hpp"

namespace prot {

SimplePrsmTableWriter::SimplePrsmTableWriter(PrsmParaPtr prsm_para_ptr, 
                                             const std::string &input_file_ext, 
                                             const std::string &output_file_ext):
    prsm_para_ptr_(prsm_para_ptr),
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext) {
    }

void SimplePrsmTableWriter::write() {
  std::string spectrum_file_name  = prsm_para_ptr_->getSpectrumFileName(); 
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string output_file_name = base_name + "." + output_file_ext_;
  std::ofstream file_; 
  file_.open(output_file_name.c_str());
  //write title
  file_ << "Spectrum_ID" << "\t"
      << "Scan(s)" << "\t"
      << "Precursor_ID" << "\t"
      << "Precursor_mass" << "\t"
      << "Score" << "\t"
      << "Protein_name"
      << std::endl;

  std::string input_file_name = base_name + "." + input_file_ext_;
  LOG_DEBUG("input file_name " << input_file_name);
  SimplePrsmReader reader(input_file_name);
  SimplePrsmPtr prsm_ptr = reader.readOnePrsm();
  LOG_DEBUG("read simple prsm complete ");
  file_ << std::setprecision(10);
  while (prsm_ptr != nullptr) {
    file_ << prsm_ptr->getSpectrumId() << "\t"
        << prsm_ptr->getSpectrumScan() << "\t"
        << prsm_ptr->getPrecursorId() << "\t"
        << prsm_ptr->getPrecMass() << "\t"
        << prsm_ptr->getScore() << "\t"
        << prsm_ptr->getSeqName()
        << std::endl;
    prsm_ptr = reader.readOnePrsm();
  }

  //write end;
  file_.close();
}

}
