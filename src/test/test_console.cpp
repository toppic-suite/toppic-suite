#include <iostream>

#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"

#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "prsm/output_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/table_writer.hpp"
#include "prsm/prsm_fdr.hpp"

#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "xpp/xml_generator.hpp"
//#include "xpp/transformer.hpp"


std::map<std::string, std::string> getArgumentMap(const std::string &search_db_file_name,
                                                  const std::string &spectrum_file_name,
                                                  const std::string &search_type,
                                                  const std::string &cutoff_type){
  std::map<std::string, std::string> arguments;
  arguments["oriDatabaseFileName"]=search_db_file_name;
  if (search_type == "TARGET") {
    arguments["databaseFileName"]=search_db_file_name;
  }
  else if (search_type == "TARGET+DECOY") {
    arguments["databaseFileName"]=search_db_file_name + "_target_decoy";
  }
  arguments["spectrumFileName"]=spectrum_file_name;
  arguments["activation"]="FILE";
  arguments["searchType"]=search_type;
  arguments["cysteineProtection"]="C0";
  arguments["shiftNumber"]="2";
  arguments["errorTolerance"]="15";
  arguments["cutoffType"]="EVALUE";
  arguments["cutoff"]="0.01";
  arguments["numOfTopPrsms"]="1";
  return arguments;
}

int main(int argc, char* argv[]) {
  try {
    prot::initBaseData();
    std::cout << "Test msalign+ 0.1 " << std::endl;

    std::map<std::string, std::string> arguments = getArgumentMap(argv[1], argv[2], argv[3], argv[4]);
    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];

    std::string ori_db_file_name = argv[2];
    if (arguments["searchType"] == "TARGET+DECOY") {
      prot::generateShuffleDb(ori_db_file_name, db_file_name);
    }

    std::cout << "Zero ptm search 0.1 " << std::endl;
    prot::ZeroPtmMngPtr zero_mng_ptr
        = prot::ZeroPtmMngPtr(new prot::ZeroPtmMng (arguments));
    zero_mng_ptr->output_file_ext_ = "ZERO";
    prot::zeroPtmSearchProcess(zero_mng_ptr);


    std::cout << "Fast filter 0.1 " << std::endl;
    prot::PtmFastFilterMngPtr filter_mng_ptr 
        = prot::PtmFastFilterMngPtr(new prot::PtmFastFilterMng(arguments));
    filter_mng_ptr->output_file_ext_ = "FILTER";
    prot::PtmFastFilterProcessor filter_processor(filter_mng_ptr);
    filter_processor.process();


    std::cout << "Ptm alignment 0.1 " << std::endl;
    prot::PtmMngPtr ptm_mng_ptr 
        = prot::PtmMngPtr(new prot::PtmMng(arguments));
    ptm_mng_ptr->input_file_ext_ = "FILTER_COMBINED";
    ptm_mng_ptr->output_file_ext_ = "PTM";
    prot::PtmProcessor ptm_processor(ptm_mng_ptr);
    ptm_processor.process();

    std::cout << "Combine prsms " << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("ZERO");
    input_exts.push_back("PTM");
    prot::PrSMCombine combine_processor(db_file_name, sp_file_name, 
                                        input_exts, "RAW_RESULT");
    combine_processor.process();

    std::cout << "EValueConsole 0.1 " << std::endl;
    prot::TdgfMngPtr tdgf_mng_ptr = 
        prot::TdgfMngPtr(new prot::TdgfMng (arguments));
    tdgf_mng_ptr->input_file_ext_="RAW_RESULT";
    tdgf_mng_ptr->output_file_ext_ = "EVALUE";
    prot::EValueProcessor processor(tdgf_mng_ptr);
    processor.init();
    // compute E-value for a set of prsms each run 
    processor.process(false);

    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top selector 0.1 " << std::endl;
      int n_top;
      std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
      prot::PrSMSelector selector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top);
      selector.process();
    }
    else {
      std::cout << "Top selector 0.1 " << std::endl;
      int n_top;
      std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
      std::cout << "argument " << arguments["numOfTopPrsms"] << " n top " << n_top << std::endl;
      prot::PrSMSelector selector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top);
      selector.process();

      std::cout << "FDR selector 0.1 " << std::endl;
      prot::PrSMFdr fdr(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr.process();
    }

    std::cout << "Cutoff selector 0.1 " << std::endl;
    std::string cutoff_type = arguments["cutoffValue"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    prot::OutputSelector output_selector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                                         cutoff_type, cutoff_value);
    output_selector.process();

    std::cout << "Set species 0.1 " << std::endl;
    double ppo;
    std::istringstream (arguments["error_tolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    prot::PrsmSpecies prsm_species(db_file_name, sp_file_name, "CUTOFF_RESULT", 
                                   "OUTPUT_RESULT", ppo);
    prsm_species.process();

    std::cout << "Table output 0.1" << std::endl;
    prot::TableWriter table_out(db_file_name, sp_file_name, "OUTPUT_RESULT",
                                "OUTPUT_TABLE");
    table_out.write();

    /*
    std::cout << "dbfilename " << arguments["databaseFileName"] << std::endl;
    prot::XmlGenerator xml_gene = prot::XmlGenerator(arguments, "OUTPUT_RESULT");
    xml_gene.process();

    TransformerPtr trans = TransformerPtr(new Transformer());
    trans->trans();
    */
  } catch (const char* e) {
    std::cout << "Exception " << e << std::endl;
  }
  return 0;
}
