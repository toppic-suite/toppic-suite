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
#include <iomanip>

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_prob.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "console/argument.hpp"

namespace prot {

int testProcess(int argc, char* argv[]) {
  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC test " << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    BaseData::init(exe_dir);

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    int n_top;
    std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
    int shift_num;
    std::istringstream (arguments["shiftNumber"]) >> shift_num;
    double max_ptm_mass;
    std::istringstream (arguments["maxPtmMass"]) >> max_ptm_mass;

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    if (arguments["searchType"] == "TARGET+DECOY") {
      FastaUtil::dbPreprocess(ori_db_file_name, db_file_name, false, db_block_size);
    }

    double K1 = 0.03125;
    double K2 = 0.03125;
    double pref = 1.0;
    for (double inte = 10; inte > 0.0001; inte /= 2) {
      PrsmProbPtr prob_processor
          = std::make_shared<PrsmProb>(db_file_name, sp_file_name,
                                       prsm_para_ptr->getFixModPtrVec(),
                                       "EVALUE", "EVALUE_ADJUST", K1, K2, pref, inte);
      prob_processor->process();
      prob_processor = nullptr;
      std::cout << "Combining prsms finished." << std::endl;

      if (arguments["searchType"]=="TARGET") { 
        std::cout << "Top prsm selecting starts" << std::endl;
        PrsmTopSelectorPtr selector = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "EVALUE_ADJUST", "TOP", n_top);
        selector->process();
        selector = nullptr;
        std::cout << "Top prsm selecting finished." << std::endl;
      } else {
        std::cout << "Top prsm selecting starts " << std::endl;
        PrsmTopSelectorPtr selector = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "EVALUE_ADJUST", "TOP_PRE", n_top);
        selector->process();
        selector = nullptr;
        std::cout << "Top prsm selecting finished." << std::endl;

        std::cout << "FDR computation starts " << std::endl;
        PrsmFdrPtr fdr = PrsmFdrPtr(new PrsmFdr(db_file_name, sp_file_name, "TOP_PRE", "TOP"));
        fdr->process();
        fdr = nullptr;
        std::cout << "FDR computation finished." << std::endl;
      }

      std::cout << "Prsm cutoff selecting starts " << std::endl;
      std::string cutoff_type = arguments["cutoffValue"];
      double cutoff_value;
      std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
      PrsmCutoffSelectorPtr output_selector
          = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP",
                                                 "OUTPUT_RESULT_" + std::to_string(K1) + "_" + std::to_string(K2) + "_" + std::to_string(pref) + "_" + std::to_string(inte), 
                                                 cutoff_type, cutoff_value);
      output_selector->process();
      output_selector = nullptr;
      std::cout << "Prsm cutoff selecting finished." << std::endl;

      std::cout << "Outputting table starts " << std::endl;
      PrsmTableWriterPtr table_out
          = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments,
                                              "OUTPUT_RESULT_" + std::to_string(K1) + "_" + std::to_string(K2) + "_" + std::to_string(pref) + "_" + std::to_string(inte), 
                                              "OUTPUT_TABLE_" + std::to_string(K1) + "_" + std::to_string(K2) + "_" + std::to_string(pref) + "_" + std::to_string(inte));
      table_out->write();
      table_out = nullptr;
      std::cout << "Outputting table finished." << std::endl;
    }

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  prot::log_level = 2;
  std::cout << std::setprecision(10);
  return prot::testProcess(argc, argv);
}
