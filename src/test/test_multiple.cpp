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
#include "prsm/prsm_coverage.hpp"

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
    /*
    if (arguments["searchType"] == "TARGET+DECOY") {
      generateShuffleDb(ori_db_file_name, db_file_name);
    }
    
    std::cout << "Zero ptm searching starts " << std::endl;
    ZeroPtmMngPtr zero_mng_ptr = ZeroPtmMngPtr(new ZeroPtmMng (prsm_para_ptr, "ZERO"));
    zeroPtmSearchProcess(zero_mng_ptr);


    std::cout << "Fast filtering starts " << std::endl;
    PtmFastFilterMngPtr filter_mng_ptr 
        = PtmFastFilterMngPtr(new PtmFastFilterMng(prsm_para_ptr, "FILTER"));
    PtmFastFilterProcessorPtr filter_processor = PtmFastFilterProcessorPtr(new PtmFastFilterProcessor(filter_mng_ptr));
    filter_processor->process();
    filter_processor = nullptr;

    int n_top;
    std::istringstream (arguments["numOfTopPrsms"]) >> n_top;
    int shift_num;
    std::istringstream (arguments["shiftNumber"]) >> shift_num;
    double max_ptm_mass;
    std::istringstream (arguments["maxPtmMass"]) >> max_ptm_mass;

    std::cout << "Ptm searching starts" << std::endl;
    PtmMngPtr ptm_mng_ptr = PtmMngPtr(new PtmMng(prsm_para_ptr, n_top, shift_num,
                                                 max_ptm_mass, "FILTER_COMBINED", "PTM"));
    PtmProcessorPtr ptm_processor = PtmProcessorPtr(new PtmProcessor(ptm_mng_ptr));
    ptm_processor->process();
    ptm_processor = nullptr;

    std::cout << "Combining prsms starts" << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("ZERO");
    input_exts.push_back("PTM");
    PrsmCombinePtr combine_processor = PrsmCombinePtr(new PrsmCombine(db_file_name, sp_file_name, 
                                                                    input_exts, "RAW_RESULT"));
    combine_processor->process();
    combine_processor = nullptr;
    std::cout << "Combining prsms finished." << std::endl;

    std::cout << "Poisson computation starts" << std::endl;
    PoissonMngPtr poisson_mng_ptr = PoissonMngPtr(new PoissonMng (prsm_para_ptr, shift_num, max_ptm_mass, 
                                                      "RAW_RESULT", "POISSON_EVALUE"));
    PoissonProcessorPtr poisson = PoissonProcessorPtr(new PoissonProcessor(poisson_mng_ptr));
    poisson->init();
    poisson->process();
    poisson = nullptr;

    std::cout << "E-value computation starts" << std::endl;
    TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, shift_num, max_ptm_mass,
                                                      "POISSON_EVALUE", "EVALUE"));
    EValueProcessorPtr processor = EValueProcessorPtr(new EValueProcessor(tdgf_mng_ptr));
    processor->init();
    // compute E-value for a set of prsms each run 
    processor->process(false);
    processor = nullptr;

    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top prsm selecting starts" << std::endl;
      PrsmSelectorPtr selector = PrsmSelectorPtr(new PrsmSelector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top));
      selector->process();
      selector = nullptr;
      std::cout << "Top prsm selecting finished." << std::endl;
    }
    else {
      std::cout << "Top prsm selecting starts " << std::endl;
      PrsmSelectorPtr selector = PrsmSelectorPtr(new PrsmSelector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top));
      selector->process();
      selector = nullptr;
      std::cout << "Top prsm selecting finished." << std::endl;

      std::cout << "FDR computation starts " << std::endl;
      PrsmFdrPtr fdr = PrsmFdrPtr(new PrsmFdr(db_file_name, sp_file_name, "TOP_PRE", "TOP"));
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation finished." << std::endl;
    }
    */
    std::cout << "Finding species starts " << std::endl;
    double ppo;
    std::istringstream (arguments["error_tolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    PrsmSpeciesPtr prsm_species = PrsmSpeciesPtr(new PrsmSpecies(db_file_name, sp_file_name, 
                                                                 "TOP", "OUTPUT_RESULT", ppo));
    prsm_species->process();
    prsm_species = nullptr;
    std::cout << "Finding species finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    TableWriterPtr table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "OUTPUT_RESULT", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;
    /*
    std::cout << "Generating view xml files starts " << std::endl;
    XmlGeneratorPtr xml_gene = XmlGeneratorPtr(new XmlGenerator(prsm_para_ptr, exe_dir,"OUTPUT_RESULT"));
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating view xml files finished." << std::endl;

    std::cout << "Converting xml files to html files starts " << std::endl;
    translate(arguments);
    std::cout << "Converting xml files to html files finished." << std::endl;
    */

    std::cout << "Identification End!" << std::endl;


    PrsmCoveragePtr prsm_coverage = PrsmCoveragePtr(new PrsmCoverage(prsm_para_ptr, "OUTPUT_RESULT", "COVERAGE"));
    prsm_coverage->processSingleCoverage();
    prsm_coverage->processCombineCoverage();
    prsm_coverage = nullptr;


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

