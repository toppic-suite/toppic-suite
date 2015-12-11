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
    double filtering_result_num = std::stod(arguments["filteringResultNumber"]);
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

    time(&start_s);
    std::cout << "Zero PTM filtering started." << std::endl;
    ZeroPtmFilterMngPtr zero_filter_mng_ptr = ZeroPtmFilterMngPtr(new ZeroPtmFilterMng (prsm_para_ptr, "ZERO_FILTER"));
    ZeroPtmFilterProcessorPtr zero_filter_processor = ZeroPtmFilterProcessorPtr(new ZeroPtmFilterProcessor(zero_filter_mng_ptr));
    zero_filter_processor->process();
    //WebLog::completeFunction(WebLog::ZeroPtmTime());
    std::cout << "Zero PTM filtering finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Zero PTM filtering running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    /*
    std::cout << "Outputting table starts " << std::endl;
    TableWriterPtr table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "ZERO_COMPLETE", "ZERO_COMPLETE_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "ZERO_PREFIX", "ZERO_PREFIX_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "ZERO_SUFFIX", "ZERO_SUFFIX_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "ZERO_INTERNAL", "ZERO_INTERNAL_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;
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
  return prot::zero_ptm_process(argc, argv);
}
