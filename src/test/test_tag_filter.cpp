#include <iostream>
#include <iomanip>

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "prsm/prsm_para.hpp"

#include "tagfilter/tag_filter_mng.hpp"
#include "tagfilter/tag_filter_processor.hpp"

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
    std::cout << "TopPIC 1.0" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    argu_processor.outputArguments(std::cout, arguments);

    BaseData::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string log_file_name = arguments["logFileName"];

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    int ptm_num = std::stoi(arguments["ptmNumber"]);
    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    int filter_result_num = 10;

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
    std::cout << "Tag filtering started." << std::endl;
    TagFilterMngPtr diag_filter_mng_ptr 
        = std::make_shared<TagFilterMng>(prsm_para_ptr, filter_result_num, 
                                         arguments["residueModFileName"], 
                                         "TAG_FILTER");

    TagFilterProcessorPtr tag_filter_processor 
        = std::make_shared<TagFilterProcessor>(diag_filter_mng_ptr);
    tag_filter_processor->process();
    std::cout << "Tag filtering finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  //prot::log_level = 2;
  std::cout << std::setprecision(10);
  return prot::two_base_opt(argc, argv);
}
