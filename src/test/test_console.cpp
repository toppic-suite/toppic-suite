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
    std::cout << "Exectutive directory is: " << exe_dir << std::endl;
    initBaseData(exe_dir);

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];

    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    if (arguments["searchType"] == "TARGET+DECOY") {
      generateShuffleDb(ori_db_file_name, db_file_name);
    }
    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));

    std::cout << "Zero ptm search " << std::endl;
    ZeroPtmMngPtr zero_mng_ptr = ZeroPtmMngPtr(new ZeroPtmMng (prsm_para_ptr, "ZERO"));
    zeroPtmSearchProcess(zero_mng_ptr);


    std::cout << "Fast filter " << std::endl;
    PtmFastFilterMngPtr filter_mng_ptr 
        = PtmFastFilterMngPtr(new PtmFastFilterMng(prsm_para_ptr, "FILTER"));
    PtmFastFilterProcessor filter_processor(filter_mng_ptr);
    filter_processor.process();

    int n_top;
    std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
    int shift_num;
    std::istringstream (arguments["shiftNumber"]) >> shift_num;

    std::cout << "Ptm alignment " << std::endl;
    PtmMngPtr ptm_mng_ptr = PtmMngPtr(new PtmMng(prsm_para_ptr, n_top, shift_num,
                                                 "FILTER_COMBINED", "PTM"));
    prot::PtmProcessor ptm_processor(ptm_mng_ptr);
    ptm_processor.process();

    std::cout << "Combine prsms " << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("ZERO");
    input_exts.push_back("PTM");
    PrSMCombine combine_processor(db_file_name, sp_file_name, 
                                  input_exts, "RAW_RESULT");
    combine_processor.process();

    std::cout << "E-value computation " << std::endl;
    TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, shift_num,
                                                      "RAW_RESULT", "EVALUE"));
    prot::EValueProcessor processor(tdgf_mng_ptr);
    processor.init();
    // compute E-value for a set of prsms each run 
    processor.process(false);

    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top selector " << std::endl;
      PrSMSelector selector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top);
      selector.process();
    }
    else {
      std::cout << "Top selector " << std::endl;
      PrSMSelector selector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top);
      selector.process();

      std::cout << "FDR computation " << std::endl;
      PrSMFdr fdr(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr.process();
    }

    std::cout << "Cutoff selector " << std::endl;
    std::string cutoff_type = arguments["cutoffValue"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    OutputSelector output_selector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                                         cutoff_type, cutoff_value);
    output_selector.process();

    std::cout << "Find species " << std::endl;
    double ppo;
    std::istringstream (arguments["error_tolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    PrsmSpecies prsm_species(db_file_name, sp_file_name, "CUTOFF_RESULT", 
                                   "OUTPUT_RESULT", ppo);
    prsm_species.process();

    std::cout << "Table output " << std::endl;
    TableWriter table_out(db_file_name, sp_file_name, "OUTPUT_RESULT",
                                "OUTPUT_TABLE");
    table_out.write();

    std::cout << "Generate view xml files " << std::endl;
    XmlGenerator xml_gene = XmlGenerator(prsm_para_ptr, exe_dir,"OUTPUT_RESULT");
    xml_gene.process();

    std::cout << "Convert view xml files to html files " << std::endl;
    prot::translate(arguments);

  } catch (const char* e) {
    std::cout << "Exception " << e << std::endl;
  }
  return 0;
}

}

int main(int argc, char* argv[]) {
  return prot::process(argc, argv);
}

