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
#include "base/web_logger.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"

#include "zeroptmsearch/zero_ptm_search_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"

#include "oneptmsearch/ptm_search_mng.hpp"
#include "oneptmsearch/one_ptm_search.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "ptmsearch/ptm_search_processor.hpp"

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/argument.hpp"

namespace prot {

int two_ptm_process(int argc, char* argv[]) {
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

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string log_file_name = arguments["logFileName"];

    bool use_gf = false; 
    if (arguments["useGf"] == "true") {
      use_gf = true;
    }
    // initialize log file 
    //WebLog::init(log_file_name, use_gf, ptm_num);
    LOG_DEBUG("web log inited");

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    LOG_DEBUG("prsm para inited");

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    FastaUtil::dbPreprocess (ori_db_file_name, db_file_name, decoy, db_block_size);
    MsAlignUtil::geneSpIndex(sp_file_name);

    time_t start_s;
    time_t stop_s;

    std::vector<std::string> input_exts ;

    std::cout << "Finding protein species started." << std::endl;
    double ppo;
    std::istringstream (arguments["errorTolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    ModPtrVec mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
    PrsmSpeciesPtr prsm_species
        = std::make_shared<PrsmSpecies>(db_file_name, sp_file_name, "CUTOFF_RESULT", mod_ptr_vec, "OUTPUT_RESULT", ppo);
    prsm_species->process();
    prsm_species = nullptr;
    std::cout << "Finding protein species finished." << std::endl;

    time(&start_s);
    std::cout << "Generating view xml files started." << std::endl;
    XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, exe_dir, "OUTPUT_RESULT");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating view xml files finished." << std::endl;

    std::cout << "Converting xml files to html files started." << std::endl;
    translate(arguments);
    std::cout << "Converting xml files to html files finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Html generation running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPIC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  prot::log_level = 2;
  std::cout << std::setprecision(10);
  return prot::two_ptm_process(argc, argv);
}
