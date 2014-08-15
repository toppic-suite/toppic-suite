#include <iostream>
#include <iomanip>

#include "base/fasta_reader.hpp"
#include "base/base_data.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "prsm/output_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/table_writer.hpp"
#include "prsm/prsm_fdr.hpp"

#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

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

    std::string exe_dir = arguments["executiveDir"];
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

    int start_s = clock();

    std::cout << "Zero ptm searching starts " << std::endl;
    ZeroPtmMngPtr zero_mng_ptr = ZeroPtmMngPtr(new ZeroPtmMng (prsm_para_ptr, "ZERO"));
    zeroPtmSearchProcess(zero_mng_ptr);

    int stop_s = clock();
    std::cout << std::endl << "Running time: " << (stop_s-start_s) / double(CLOCKS_PER_SEC)  << " seconds " << std::endl;

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
