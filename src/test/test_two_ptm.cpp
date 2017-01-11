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
    std::cout << "TopPIC 0.9.2" << std::endl;

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
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);

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

    LOG_DEBUG("block size 1");
    FastaUtil::dbPreprocess (ori_db_file_name, db_file_name, decoy, db_block_size);
    LOG_DEBUG("block size 2");
    MsAlignUtil::geneSpIndex(sp_file_name);

    time_t start_s;
    time_t stop_s;

    std::vector<std::string> input_exts ;

    input_exts.push_back("ZERO_PTM_COMPLETE");
    input_exts.push_back("ZERO_PTM_PREFIX");
    input_exts.push_back("ZERO_PTM_SUFFIX");
    input_exts.push_back("ZERO_PTM_INTERNAL");

    input_exts.push_back("ONE_PTM_COMPLETE");
    input_exts.push_back("ONE_PTM_PREFIX");
    input_exts.push_back("ONE_PTM_SUFFIX");
    input_exts.push_back("ONE_PTM_INTERNAL");

    time(&start_s);
    std::cout << "Two PTM search started." << std::endl;
    int shift_num = 2;
    PtmSearchMngPtr two_search_mng_ptr 
        = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass, shift_num,
                                         1, "DIAG_FILTER", "PTM");
    PtmSearchProcessorPtr processor = std::make_shared<PtmSearchProcessor>(two_search_mng_ptr);
    processor->process();
    std::cout << "Two PTM search finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Two PTM search running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;
    input_exts.push_back("PTM_2_COMPLETE");
    input_exts.push_back("PTM_2_PREFIX");
    input_exts.push_back("PTM_2_SUFFIX");
    input_exts.push_back("PTM_2_INTERNAL");

    time(&start_s);
    std::cout << "Combining PRSMs started." << std::endl;
    ptm_num = 2;
    int prsm_top_num = (ptm_num + 1) * 4;
    PrsmStrCombinePtr combine_ptr = std::make_shared<PrsmStrCombine>(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num);
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Combining prsms search running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;

    time(&start_s);
    std::cout << "E-value computation started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(prsm_para_ptr, ptm_num, max_ptm_mass, use_gf,
                                                        variable_ptm, "RAW_RESULT", "EVALUE");
    EValueProcessorPtr evalue_processor = std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    evalue_processor->init();
    // compute E-value for a set of prsms each run 
    evalue_processor->process(false);
    evalue_processor = nullptr;
    std::cout << "E-value computation finished." << std::endl;
    time(&stop_s);
    std::cout << "Computing e-values running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;

    time(&start_s);
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
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name,
                                               "TOP", "CUTOFF_RESULT", 
                                               cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PRSM selecting by cutoff finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "CUTOFF_RESULT", "TWO_TABLE");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;

    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "PTM_2_COMPLETE", "PTM_2_COMPLETE_TABLE");
    table_out->write();
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "PTM_2_PREFIX", "PTM_2_PREFIX_TABLE");
    table_out->write();
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "PTM_2_SUFFIX", "PTM_2_SUFFIX_TABLE");
    table_out->write();
    table_out = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "PTM_2_INTERNAL", "PTM_2_INTERNAL_TABLE");
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
  return prot::two_ptm_process(argc, argv);
}
