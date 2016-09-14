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

#include "base/fasta_reader.hpp"
#include "base/base_data.hpp"
#include "base/web_logger.hpp"

#include "spec/msalign_reader.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "prsm/output_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/simple_prsm_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "prsm/table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_stat.hpp"

#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

#include "poisson/poisson_processor.hpp"
#include "poisson/poisson_mng.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

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
    std::cout << "TopPC 0.9 " << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    initBaseData(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string log_file_name = arguments["logFileName"];
    
    //int n_top = std::stoi(arguments["numOfTopPrsms"]);
    int ptm_num = std::stoi(arguments["ptmNumber"]);
    //double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    bool use_gf = false; 
    if (arguments["useGf"] == "true") {
      use_gf = true;
    }
    /* initialize log file */
  	WebLog::init(log_file_name, use_gf, ptm_num);


    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));

    /*
    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    */
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    //int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    /*
    dbPreprocess (ori_db_file_name, db_file_name, decoy, db_block_size);
    generateSpIndex(sp_file_name);
    
    std::cout << "Zero PTM search started." << std::endl;
    ZeroPtmMngPtr zero_mng_ptr = ZeroPtmMngPtr(new ZeroPtmMng (prsm_para_ptr, "ZERO"));
    zeroPtmSearchProcess(zero_mng_ptr);
    std::cout << "Zero PTM search finished." << std::endl;

    std::cout << "Diagonal filtering started." << std::endl;
    DiagFilterMngPtr diag_filter_mng_ptr 
        = DiagFilterMngPtr(new DiagFilterMng(prsm_para_ptr, "DIAG_FILTER"));
    DiagFilterProcessorPtr diag_filter_processor = DiagFilterProcessorPtr(new DiagFilterProcessor(diag_filter_mng_ptr));
    diag_filter_processor->process();
    diag_filter_processor = nullptr;
    std::cout << "Diagonal filtering finished." << std::endl;

    std::cout << "One PTM filtering started." << std::endl;
    OnePtmFilterMngPtr one_ptm_filter_mng_ptr 
        = OnePtmFilterMngPtr(new OnePtmFilterMng(prsm_para_ptr, "ONE_PTM_FILTER"));
    OnePtmFilterProcessorPtr one_ptm_filter_processor = OnePtmFilterProcessorPtr(new OnePtmFilterProcessor(one_ptm_filter_mng_ptr));
    one_ptm_filter_processor->process();
    one_ptm_filter_processor = nullptr;
    std::cout << "One PTM filtering finished." << std::endl;

    std::cout << "Combining simple PRSMs started." << std::endl;
    std::vector<std::string> simple_input_exts;
    simple_input_exts.push_back("DIAG_FILTER");
    simple_input_exts.push_back("ONE_PTM_FILTER");
    int top_num = diag_filter_mng_ptr->ptm_fast_filter_result_num_ 
        + one_ptm_filter_mng_ptr->one_ptm_filter_result_num_;
    LOG_DEBUG("top number " << top_num);
    SimplePrsmStrCombinePtr simple_combine_ptr(
        new SimplePrsmStrCombine(sp_file_name, simple_input_exts, "FILTER", top_num));
    simple_combine_ptr->process();
    simple_combine_ptr = nullptr;
    std::cout << "Combining simple PRSMs finished." << std::endl;
    
    std::cout << "PTM search started." << std::endl;
    PtmMngPtr ptm_mng_ptr = PtmMngPtr(new PtmMng(prsm_para_ptr, n_top, shift_num,
                                                 max_ptm_mass, "FILTER", "PTM"));
    PtmProcessorPtr ptm_processor = PtmProcessorPtr(new PtmProcessor(ptm_mng_ptr));
    ptm_processor->process();
    ptm_processor = nullptr;
    std::cout << "PTM search finished" << std::endl;

    std::cout << "Combining PRSMs started." << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("ZERO_COMPLETE");
    input_exts.push_back("ZERO_PREFIX");
    input_exts.push_back("ZERO_SUFFIX");
    input_exts.push_back("ZERO_INTERNAL");
    input_exts.push_back("PTM");
    int prsm_top_num = (shift_num + 1) * 4;
    PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num));
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs finished." << std::endl;

    std::cout << "E-value computation started." << std::endl;
    TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, shift_num, max_ptm_mass,
                                                      "RAW_RESULT", "EVALUE"));
    if (arguments["useTable"] == "false")
      tdgf_mng_ptr->use_table = false;
                                                            
    EValueProcessorPtr processor = EValueProcessorPtr(new EValueProcessor(tdgf_mng_ptr));
    processor->init();
    // compute E-value for a set of prsms each run 
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation finished." << std::endl;

    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top PRSM selecting started" << std::endl;
      PrsmSelectorPtr selector = PrsmSelectorPtr(new PrsmSelector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top));
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting finished." << std::endl;
    }
    else {
      std::cout << "Top PRSM selecting started " << std::endl;
      PrsmSelectorPtr selector = PrsmSelectorPtr(new PrsmSelector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top));
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting finished." << std::endl;

      std::cout << "FDR computation started. " << std::endl;
      PrsmFdrPtr fdr = PrsmFdrPtr(new PrsmFdr(db_file_name, sp_file_name, "TOP_PRE", "TOP"));
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation finished." << std::endl;
    }
    
    std::cout << "PRSM selecting by cutoff started." << std::endl;
    std::string cutoff_type = arguments["cutoffType"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    OutputSelectorPtr output_selector = OutputSelectorPtr(
        new OutputSelector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                           cutoff_type, cutoff_value));
    output_selector->process();
    output_selector = nullptr;
    std::cout << "PRSM selecting by cutoff finished." << std::endl;

    std::cout << "Finding protein species started." << std::endl;
    double ppo;
    std::istringstream (arguments["errorTolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    ResiduePtrVec residue_ptr_vec = prsm_para_ptr->getFixModResiduePtrVec();
    PrsmSpeciesPtr prsm_species = PrsmSpeciesPtr(
        new PrsmSpecies(db_file_name, sp_file_name, "CUTOFF_RESULT", "OUTPUT_RESULT", residue_ptr_vec, ppo));
    prsm_species->process();
    prsm_species = nullptr;
    std::cout << "Finding protein species finished." << std::endl;

    std::cout << "Outputting table started." << std::endl;
    TableWriterPtr table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "OUTPUT_RESULT", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;
    */

    std::cout << "Generating view xml files started." << std::endl;
    XmlGeneratorPtr xml_gene = XmlGeneratorPtr(new XmlGenerator(prsm_para_ptr, exe_dir,"OUTPUT_RESULT"));
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating view xml files finished." << std::endl;

    std::cout << "Converting xml files to html files started." << std::endl;
    translate(arguments);
    std::cout << "Converting xml files to html files finished." << std::endl;
    
    /*
    if (arguments["keepTempFiles"] != "true"){
      std::cout << "Deleting temporary files started." << std::endl;
      delDir(basename(sp_file_name) + "_xml");
      delFile(exe_dir + "/run.log");
      cleanDir(sp_file_name);
      cleanDir(db_file_name);	  
      std::cout << "Deleting temporary files finished." << std::endl;
    }  
    */

    std::cout << "Statistics started." << std::endl;
    PrsmStatPtr stat_ptr = PrsmStatPtr(new PrsmStat(prsm_para_ptr, "OUTPUT_RESULT", "STAT"));
    stat_ptr->process();
    stat_ptr = nullptr;
    std::cout << "Statistics finished." << std::endl;
    
    WebLog::close();

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPIC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  //prot::log_level = 2;
  return prot::process(argc, argv);
}

