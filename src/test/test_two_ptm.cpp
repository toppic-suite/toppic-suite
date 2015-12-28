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

#include "oneptmsearch/ptm_search_mng.hpp"

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
  	WebLog::init(log_file_name, use_gf, ptm_num);
    LOG_DEBUG("web log inited");

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));
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

    /*
    time(&start_s);
    std::cout << "Diagonal PTM filtering started." << std::endl;
    DiagFilterMngPtr diag_filter_mng_ptr 
        = DiagFilterMngPtr(new DiagFilterMng (prsm_para_ptr, filter_result_num, "DIAG_FILTER"));
    DiagFilterProcessorPtr diag_filter_processor 
        = DiagFilterProcessorPtr(new DiagFilterProcessor(diag_filter_mng_ptr));
    diag_filter_processor->process();
    //WebLog::completeFunction(WebLog::ZeroPtmTime());
    std::cout << "Diagonal filtering finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Diagonal filtering running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;
    */

    time(&start_s);
    std::cout << "Two PTM search started." << std::endl;
    int shift_num = 2;
    PtmSearchMngPtr two_search_mng_ptr 
        = PtmSearchMngPtr(new PtmSearchMng (prsm_para_ptr, n_top, max_ptm_mass, shift_num,
                                            "DIAG_FILTER", "PTM"));
    PtmSearchProcessorPtr processor = PtmSearchProcessorPtr(new PtmSearchProcessor(two_search_mng_ptr));
    processor->process();
    std::cout << "Two PTM search finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Two PTM search running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    /*
    time(&start_s);
    std::cout << "E-value computation started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, ptm_num, max_ptm_mass, use_gf,
                                                      variable_ptm, "ZERO", "EVALUE"));
    EValueProcessorPtr processor = EValueProcessorPtr(new EValueProcessor(tdgf_mng_ptr));
    processor->init();
    // compute E-value for a set of prsms each run 
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Computing e-values running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    std::cout << "Top PRSM selecting started" << std::endl;
    PrsmTopSelectorPtr selector = PrsmTopSelectorPtr(
        new PrsmTopSelector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top));
    selector->process();
    selector = nullptr;
    std::cout << "Top PRSM selecting finished." << std::endl;

    std::cout << "PRSM selecting by cutoff started." << std::endl;
    std::string cutoff_type = arguments["cutoffType"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector = PrsmCutoffSelectorPtr(
        new PrsmCutoffSelector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                           cutoff_type, cutoff_value));
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PRSM selecting by cutoff finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, "CUTOFF_RESULT", "ZERO_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;
    */
    /*
    table_out = PrsmTableWriterPtr(new PrsmTableWriter(prsm_para_ptr, "ZERO_COMPLETE", "ZERO_COMPLETE_TABLE"));
    table_out->write();
    table_out = PrsmTableWriterPtr(new PrsmTableWriter(prsm_para_ptr, "ZERO_PREFIX", "ZERO_PREFIX_TABLE"));
    table_out->write();
    table_out = PrsmTableWriterPtr(new PrsmTableWriter(prsm_para_ptr, "ZERO_SUFFIX", "ZERO_SUFFIX_TABLE"));
    table_out->write();
    table_out = PrsmTableWriterPtr(new PrsmTableWriter(prsm_para_ptr, "ZERO_INTERNAL", "ZERO_INTERNAL_TABLE"));
    table_out->write();
    */

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
  return prot::two_ptm_process(argc, argv);
}
