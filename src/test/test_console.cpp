#include <iostream>

#include "base/fasta_reader.hpp"
#include "base/species.hpp"
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

#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

#include "poisson/poisson_processor.hpp"
#include "poisson/poisson_mng.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "xpp/xml_generator.hpp"
#include "xpp/transformer.hpp"
#include "test/argument.hpp"

namespace prot {

int process(int argc, char* argv[]) {
  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopId 0.5 " << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    initBaseData(exe_dir);

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));
    if (arguments["searchType"] == "TARGET+DECOY") {
      generateShuffleDb(ori_db_file_name, db_file_name);
    }
    /*
    struct timeval start_time; 
    struct timeval end_time; 
    float duration;
    gettimeofday(&start_time, NULL);
    */
    std::cout << "Zero ptm searching starts " << std::endl;
    ZeroPtmMngPtr zero_mng_ptr = ZeroPtmMngPtr(new ZeroPtmMng (prsm_para_ptr, "ZERO"));
    zeroPtmSearchProcess(zero_mng_ptr);


    std::cout << "Fast filtering starts " << std::endl;
    PtmFastFilterMngPtr filter_mng_ptr 
        = PtmFastFilterMngPtr(new PtmFastFilterMng(prsm_para_ptr, "FILTER"));
    PtmFastFilterProcessor filter_processor(filter_mng_ptr);
    filter_processor.process();

    int n_top;
    std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
    int shift_num;
    std::istringstream (arguments["shiftNumber"]) >> shift_num;
    double max_ptm_mass;
    std::istringstream (arguments["maxPtmMass"]) >> max_ptm_mass;

    std::cout << "Ptm searching starts" << std::endl;
    PtmMngPtr ptm_mng_ptr = PtmMngPtr(new PtmMng(prsm_para_ptr, n_top, shift_num,
                                                 max_ptm_mass, "FILTER_COMBINED", "PTM"));
    prot::PtmProcessor ptm_processor(ptm_mng_ptr);
    ptm_processor.process();

    std::cout << "Combining prsms starts" << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("ZERO");
    input_exts.push_back("PTM");
    PrsmCombine combine_processor(db_file_name, sp_file_name, 
                                  input_exts, "RAW_RESULT");
    combine_processor.process();
    std::cout << "Combining prsms finished." << std::endl;

    std::cout << "Poisson computation starts" << std::endl;
    PoissonMngPtr poisson_mng_ptr = PoissonMngPtr(new PoissonMng (prsm_para_ptr, shift_num, max_ptm_mass, 
                                                      "RAW_RESULT", "POISSON_EVALUE"));
    prot::PoissonProcessor poisson(poisson_mng_ptr);
    poisson.init();
    poisson.process();

    std::cout << "E-value computation starts" << std::endl;
    TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, shift_num, max_ptm_mass,
                                                      "POISSON_EVALUE", "EVALUE"));
    prot::EValueProcessor processor(tdgf_mng_ptr);
    processor.init();
    // compute E-value for a set of prsms each run 
    processor.process(false);

    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top prsm selecting starts" << std::endl;
      PrsmSelector selector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top);
      selector.process();
      std::cout << "Top prsm selecting finished." << std::endl;
    }
    else {
      std::cout << "Top prsm selecting starts " << std::endl;
      PrsmSelector selector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top);
      selector.process();
      std::cout << "Top prsm selecting finished." << std::endl;

      std::cout << "FDR computation starts " << std::endl;
      PrsmFdr fdr(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr.process();
      std::cout << "FDR computation finished." << std::endl;
    }

    std::cout << "Prsm cutoff selecting starts " << std::endl;
    std::string cutoff_type = arguments["cutoffValue"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    OutputSelector output_selector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                                         cutoff_type, cutoff_value);
    output_selector.process();
    std::cout << "Prsm cutoff selecting finished." << std::endl;

    std::cout << "Finding species starts " << std::endl;
    double ppo;
    std::istringstream (arguments["error_tolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    PrsmSpecies prsm_species(db_file_name, sp_file_name, "CUTOFF_RESULT", 
                                   "OUTPUT_RESULT", ppo);
    prsm_species.process();
    std::cout << "Finding species finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    TableWriter table_out(prsm_para_ptr, "OUTPUT_RESULT", "OUTPUT_TABLE");
    table_out.write();
    std::cout << "Outputting table finished." << std::endl;

    std::cout << "Generating view xml files starts " << std::endl;
    XmlGenerator xml_gene = XmlGenerator(prsm_para_ptr, exe_dir,"OUTPUT_RESULT");
    xml_gene.process();
    std::cout << "Generating view xml files finished." << std::endl;

    std::cout << "Converting xml files to html files starts " << std::endl;
    prot::translate(arguments);
    std::cout << "Converting xml files to html files finished." << std::endl;

    std::cout << "Identification End!" << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "Topid finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  return prot::process(argc, argv);
}

