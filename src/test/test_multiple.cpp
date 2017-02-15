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


#include <iostream>

#include "base/fasta_reader.hpp"
#include "base/base_data.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_coverage.hpp"

#include "console/argument.hpp"

namespace prot {

int process(int argc, char* argv[]) {
  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC test" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    BaseData::init(exe_dir);

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    std::cout << "Finding species starts " << std::endl;
    double ppo;
    std::istringstream (arguments["error_tolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    PrsmSpeciesPtr prsm_species
        = std::make_shared<PrsmSpecies>(db_file_name, sp_file_name, 
                                        "TOP", prsm_para_ptr->getFixModPtrVec(),
                                        "OUTPUT_RESULT", ppo);
    prsm_species->process();
    prsm_species = nullptr;
    std::cout << "Finding species finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "OUTPUT_RESULT", "OUTPUT_TABLE");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;

    std::cout << "Identification End!" << std::endl;

    PrsmCoveragePtr prsm_coverage = std::make_shared<PrsmCoverage>(prsm_para_ptr, "OUTPUT_RESULT", "COVERAGE");
    prsm_coverage->processSingleCoverage();
    //prsm_coverage->processCombineCoverage();
    prsm_coverage = nullptr;


  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "Topid finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  return prot::process(argc, argv);
}

