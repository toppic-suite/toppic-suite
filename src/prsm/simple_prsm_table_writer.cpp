// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


/*
 * simple_table_writer.cpp
 *
 *  Created on: Aug 11, 2014
 *      Author: qkou
 */

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
  std::string base_name = FileUtil::basename(spectrum_file_name);
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
