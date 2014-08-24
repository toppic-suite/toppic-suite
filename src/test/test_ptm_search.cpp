#include <iostream>
#include <iomanip>

#include "base/fasta_reader.hpp"
#include "base/base_data.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/simple_prsm_table_writer.hpp"
#include "prsm/table_writer.hpp"

#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

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
    std::cout << "TopPC 0.5 " << std::endl;

    std::string exe_dir = "";
    exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    initBaseData(exe_dir);

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    int n_top;
    std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
    int shift_num;
    std::istringstream (arguments["shiftNumber"]) >> shift_num;
    double max_ptm_mass;
    std::istringstream (arguments["maxPtmMass"]) >> max_ptm_mass;

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));

    if (arguments["searchType"] == "TARGET+DECOY") {
      generateShuffleDb(ori_db_file_name, db_file_name);
    }

    std::cout << "Fast filtering starts " << std::endl;
    PtmFastFilterMngPtr filter_mng_ptr 
        = PtmFastFilterMngPtr(new PtmFastFilterMng(prsm_para_ptr, "FILTER"));
    PtmFastFilterProcessorPtr filter_processor = PtmFastFilterProcessorPtr(new PtmFastFilterProcessor(filter_mng_ptr));
    filter_processor->process();
    filter_processor = nullptr;

    long start_s = clock();
    std::cout << "Ptm searching starts" << std::endl;
    PtmMngPtr ptm_mng_ptr = PtmMngPtr(new PtmMng(prsm_para_ptr, n_top, shift_num,
                                                 max_ptm_mass, "FILTER_COMBINED", "PTM"));
    PtmProcessorPtr ptm_processor = PtmProcessorPtr(new PtmProcessor(ptm_mng_ptr));
    ptm_processor->process();
    ptm_processor = nullptr;

    long stop_s = clock();
    std::cout << std::endl << "Running time: " << (stop_s-start_s) / double(CLOCKS_PER_SEC)  << " seconds " << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    TableWriterPtr table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_1_COMPLETE", "PTM_1_COMPLETE_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_1_PREFIX", "PTM_1_PREFIX_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_1_SUFFIX", "PTM_1_SUFFIX_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_1_INTERNAL", "PTM_1_INTERNAL_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_2_COMPLETE", "PTM_2_COMPLETE_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_2_PREFIX", "PTM_2_PREFIX_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_2_SUFFIX", "PTM_2_SUFFIX_TABLE"));
    table_out->write();
    table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "PTM_2_INTERNAL", "PTM_2_INTERNAL_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;

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
