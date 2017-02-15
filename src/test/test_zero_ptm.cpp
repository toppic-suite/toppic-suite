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

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"

#include "console/argument.hpp"

namespace prot {

int zero_ptm_process(int argc, char* argv[]) {
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

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    int ptm_num = std::stoi(arguments["ptmNumber"]);
    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    bool use_gf = false; 
    if (arguments["useGf"] == "true") {
      use_gf = true;
    }
    // initialize log file 
  	WebLog::init(log_file_name, use_gf, false, ptm_num);
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

    time(&start_s);
    std::cout << "Zero PTM filtering started." << std::endl;
    ZeroPtmFilterMngPtr zero_filter_mng_ptr = ZeroPtmFilterMngPtr(new ZeroPtmFilterMng (prsm_para_ptr, "ZERO_FILTER"));
    ZeroPtmFilterProcessorPtr zero_filter_processor = ZeroPtmFilterProcessorPtr(new ZeroPtmFilterProcessor(zero_filter_mng_ptr));
    zero_filter_processor->process();
    //WebLog::completeFunction(WebLog::ZeroPtmTime());
    std::cout << "Zero PTM filtering finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Zero PTM filtering running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    time(&start_s);
    std::cout << "Zero PTM search started." << std::endl;
    ZeroPtmSearchMngPtr zero_search_mng_ptr = std::make_shared<ZeroPtmSearchMng>(prsm_para_ptr, "ZERO_FILTER", "ZERO");
    ZeroPtmSearch::process(zero_search_mng_ptr);
    std::cout << "Zero PTM search finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Zero PTM search running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    time(&start_s);
    std::cout << "E-value computation started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(prsm_para_ptr, ptm_num, max_ptm_mass, use_gf,
                                                        variable_ptm, "ZERO", "EVALUE");
    EValueProcessorPtr processor = std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run 
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Computing e-values running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top PRSM selecting started" << std::endl;
      PrsmTopSelectorPtr selector = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "EVALUE", "TOP", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting finished." << std::endl;
    } else {
      std::cout << "Top PRSM selecting started " << std::endl;
      PrsmTopSelectorPtr selector = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting finished." << std::endl;

      std::cout << "FDR computation started. " << std::endl;
      PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation finished." << std::endl;

    }

    std::cout << "PRSM selecting by cutoff started." << std::endl;
    std::string cutoff_type = arguments["cutoffType"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                                               cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PRSM selecting by cutoff finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "CUTOFF_RESULT", "ZERO_TABLE");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "ZERO_COMPLETE", "ZERO_COMPLETE_TABLE");
    table_out->write();
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "ZERO_PREFIX", "ZERO_PREFIX_TABLE");
    table_out->write();
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "ZERO_SUFFIX", "ZERO_SUFFIX_TABLE");
    table_out->write();
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "ZERO_INTERNAL", "ZERO_INTERNAL_TABLE");
    table_out->write();

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
  return prot::zero_ptm_process(argc, argv);
}
