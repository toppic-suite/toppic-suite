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
#include "prsm/prsm_feature_species.hpp"
#include "prsm/prsm_form_filter.hpp"

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

#include "local/local_mng.hpp"
#include "local/local_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/argument.hpp"

namespace prot {

int two_base_opt(int argc, char* argv[]) {
  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC 1.1.0 (" << __DATE__ << ")" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    time_t start = time(0);
    char buf[50];
    arguments["start_time"] = std::string(ctime_r(&start, buf));
    Argument::outputArguments(std::cout, arguments);

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
    int thread_num = std::stoi(arguments["threadNumber"]);

    bool use_gf = false;
    if (arguments["useGf"] == "true") {
      use_gf = true;
    }
    bool localization = false;
    if (arguments["residueModFileName"] != "") {
      localization = true;
    }
    // initialize log file
    WebLog::init(log_file_name, use_gf, localization, ptm_num);
    LOG_DEBUG("web log inited");

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    LOG_DEBUG("prsm para inited");

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    FastaUtil::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    MsAlignUtil::geneSpIndex(sp_file_name);

    std::vector<std::string> input_exts;

    std::cout << "Zero PTM filtering - started." << std::endl;
    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, "ZERO_FILTER");
    ZeroPtmFilterProcessorPtr zero_filter_processor
        = std::make_shared<ZeroPtmFilterProcessor>(zero_filter_mng_ptr);
    zero_filter_processor->process();
    WebLog::completeFunction(WebLog::ZeroPtmFilterTime());
    std::cout << "Zero PTM filtering - finished." << std::endl;

    std::cout << "Zero PTM search - started." << std::endl;
    ZeroPtmSearchMngPtr zero_search_mng_ptr
        = std::make_shared<ZeroPtmSearchMng>(prsm_para_ptr, "ZERO_FILTER", "ZERO_PTM");
    ZeroPtmSearch::process(zero_search_mng_ptr);
    WebLog::completeFunction(WebLog::ZeroPtmSearchTime());
    std::cout << "Zero PTM search - finished." << std::endl;

    input_exts.push_back("ZERO_PTM_COMPLETE");
    input_exts.push_back("ZERO_PTM_PREFIX");
    input_exts.push_back("ZERO_PTM_SUFFIX");
    input_exts.push_back("ZERO_PTM_INTERNAL");

    if (ptm_num >= 1) {
      std::cout << "One PTM filtering - started." << std::endl;
      OnePtmFilterMngPtr one_ptm_filter_mng_ptr
          = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "ONE_PTM_FILTER");
      OnePtmFilterProcessorPtr one_filter_processor
          = std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);
      one_filter_processor->process();
      WebLog::completeFunction(WebLog::OnePtmFilterTime());
      std::cout << "One PTM filtering - finished." << std::endl;

      std::cout << "One PTM search - started." << std::endl;
      int shift_num = 1;
      PtmSearchMngPtr one_search_mng_ptr
          = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass,
                                           shift_num, thread_num, "ONE_PTM_FILTER", "ONE_PTM");
      OnePtmSearch::process(one_search_mng_ptr);
      WebLog::completeFunction(WebLog::OnePtmSearchTime());
      std::cout << "One PTM search - finished." << std::endl;

      input_exts.push_back("ONE_PTM_COMPLETE");
      input_exts.push_back("ONE_PTM_PREFIX");
      input_exts.push_back("ONE_PTM_SUFFIX");
      input_exts.push_back("ONE_PTM_INTERNAL");
    }

    if (ptm_num >= 2) {
      std::cout << "Diagonal PTM filtering - started." << std::endl;
      DiagFilterMngPtr diag_filter_mng_ptr
          = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num,
                                            thread_num, "DIAG_FILTER");
      DiagFilterProcessorPtr diag_filter_processor
          = std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      WebLog::completeFunction(WebLog::DiagFilterTime());
      std::cout << "Diagonal filtering - finished." << std::endl;

      std::cout << "Two PTM search - started." << std::endl;
      PtmSearchMngPtr two_search_mng_ptr
          = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass, ptm_num,
                                           thread_num, "DIAG_FILTER", "PTM");
      PtmSearchProcessorPtr processor = std::make_shared<PtmSearchProcessor>(two_search_mng_ptr);
      processor->process();
      WebLog::completeFunction(WebLog::TwoPtmSearchTime());
      std::cout << "Two PTM search - finished." << std::endl;
      input_exts.push_back("PTM");
    }

    std::cout << "Combining PRSMs - started." << std::endl;
    int prsm_top_num = (ptm_num + 1) * 4;
    PrsmStrCombinePtr combine_ptr
        = std::make_shared<PrsmStrCombine>(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num);
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs - finished." << std::endl;

    std::cout << "E-value computation - started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr
        = std::make_shared<TdgfMng>(prsm_para_ptr, ptm_num, max_ptm_mass,
                                    use_gf, variable_ptm, "RAW_RESULT", "EVALUE");
    EValueProcessorPtr processor = std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run
    processor->process(false);
    processor = nullptr;
    if (use_gf) {
      WebLog::completeFunction(WebLog::GfEvalueTime());
    } else {
      WebLog::completeFunction(WebLog::TableEvalueTime());
    }
    std::cout << "E-value computation - finished." << std::endl;

    if (arguments["searchType"] == "TARGET") {
      std::cout << "Top PRSM selecting - started" << std::endl;
      PrsmTopSelectorPtr selector
          = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "EVALUE", "TOP", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting - finished." << std::endl;
    } else {
      std::cout << "Top PRSM selecting - started " << std::endl;
      PrsmTopSelectorPtr selector
          = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name,
                                              "EVALUE", "TOP_PRE", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting - finished." << std::endl;

      std::cout << "FDR computation - started. " << std::endl;
      PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation - finished." << std::endl;
    }

    std::cout << "PRSM selecting by cutoff - started." << std::endl;
    std::string cutoff_type = arguments["cutoffType"];
    double cutoff_value;
    std::istringstream(arguments["cutoffValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP",
                                               "CUTOFF_RESULT", cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PRSM selecting by cutoff - finished." << std::endl;

    std::string suffix = "CUTOFF_RESULT";

    if (localization) {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng = LocalMngPtr(
          new LocalMng(prsm_para_ptr, arguments["local_threshold"],
                       arguments["residueModFileName"], max_ptm_mass,
                       "CUTOFF_RESULT", "LOCAL_RESULT"));
      LocalProcessorPtr local_ptr = LocalProcessorPtr(new LocalProcessor(local_mng));
      local_ptr->process();
      local_ptr = nullptr;
      WebLog::completeFunction(WebLog::LocalizationTime());
      std::cout << "PTM characterization - finished." << std::endl;
      suffix = "LOCAL_RESULT";
    }

    std::cout << "Finding protein species - started." << std::endl;
    double ppo;
    std::istringstream(arguments["errorTolerance"]) >> ppo;
    ppo = ppo / 1000000.0;
    ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
    PrsmFeatureSpeciesPtr prsm_forms = PrsmFeatureSpeciesPtr(
        new PrsmFeatureSpecies(db_file_name, sp_file_name, "CUTOFF_RESULT",
                               "FORMS", fix_mod_list));
    prsm_forms->process();
    prsm_forms = nullptr;
    std::cout << "Finding protein species - finished." << std::endl;
    WebLog::completeFunction(WebLog::SelectingTime());

    time_t end = time(0);
    arguments["end_time"] = std::string(ctime_r(&end, buf));
    arguments["running_time"] = std::to_string(static_cast<int>(difftime(end, start)));

    std::cout << "Outputting the result table - started." << std::endl;
    PrsmTableWriterPtr table_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, arguments, "FORMS", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting the result table - finished." << std::endl;

    std::cout << "Generating xml files started." << std::endl;
    XmlGeneratorPtr xml_gene = XmlGeneratorPtr(new XmlGenerator(prsm_para_ptr, exe_dir, "FORMS"));
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating xml files - finished." << std::endl;

    std::cout << "Converting xml files to html files - started." << std::endl;
    translate(arguments);
    std::cout << "Converting xml files to html files - finished." << std::endl;
    WebLog::completeFunction(WebLog::OutPutTime());

    std::cout << "PRSM proteoform filtering - started." << std::endl;
    PrsmFormFilterPtr form_filter
        = std::make_shared<PrsmFormFilter>(db_file_name, sp_file_name, "FORMS",
                                           "FORM_FILTER_RESULT", "FORM_RESULT");
    form_filter->process();
    form_filter = nullptr;
    std::cout << "PRSM proteoform filtering - finished." << std::endl;

    std::cout << "Outputting the proteoform table - started." << std::endl;
    PrsmTableWriterPtr form_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, arguments, "FORM_RESULT", "FORM_OUTPUT_TABLE"));
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting the proteoform table - finished." << std::endl;

    if (arguments["keepTempFiles"] != "true") {
      std::cout << "Deleting temporary files - started." << std::endl;
      FileUtil::delDir(FileUtil::basename(sp_file_name) + "_xml");
      FileUtil::delFile(exe_dir + "/run.log");
      FileUtil::cleanDir(ori_db_file_name, sp_file_name);
      std::cout << "Deleting temporary files - finished." << std::endl;
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPIC finished." << std::endl;
  return 0;
}

}  // namespace prot

int main(int argc, char* argv[]) {
  // prot::log_level = 2;
  std::cout << std::setprecision(10);
  return prot::two_base_opt(argc, argv);
}
