#include <iostream>

#include "base/fasta_reader.hpp"
#include "base/base_data.hpp"
#include "base/web_logger.hpp"

#include "spec/msalign_reader.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "prsm/output_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/simple_prsm_writer.hpp"
#include "prsm/simple_prsm_str_combine.hpp"
#include "prsm/table_writer.hpp"
#include "prsm/prsm_fdr.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"

#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

#include "poisson/poisson_processor.hpp"
#include "poisson/poisson_mng.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/argument.hpp"
#include "console/summary.hpp"

namespace prot {

int process(int argc, char* argv[]) {
  try {
    boost::posix_time::ptime start_time = boost::posix_time::microsec_clock::local_time();

    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC 0.9.2" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    initBaseData(exe_dir);

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
    /* initialize log file */
  	WebLog::init(log_file_name, use_gf, ptm_num);

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    dbPreprocess (ori_db_file_name, db_file_name, decoy, db_block_size);
    generateSpIndex(sp_file_name);

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

    
    time(&start_s);
    std::cout << "Zero PTM search started." << std::endl;
    ZeroPtmMngPtr zero_mng_ptr = ZeroPtmMngPtr(new ZeroPtmMng (prsm_para_ptr, "ZERO_FILTER", "ZERO"));
    zeroPtmSearchProcess(zero_mng_ptr);
    WebLog::completeFunction(WebLog::ZeroPtmTime());
    std::cout << "Zero PTM search finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Zero PTM search running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    /*
    std::vector<std::string> filtering_result_exts;
    int total_filtering_result_num = 0;
    if (ptm_num >= 1) {
      time(&start_s);
      std::cout << "One PTM filtering started." << std::endl;
      LOG_INFO("filtering result num " << filtering_result_num);
      OnePtmFilterMngPtr one_ptm_filter_mng_ptr 
          = OnePtmFilterMngPtr(new OnePtmFilterMng(prsm_para_ptr, filtering_result_num, "ONE_PTM_FILTER"));
      OnePtmFilterProcessorPtr one_ptm_filter_processor = OnePtmFilterProcessorPtr(new OnePtmFilterProcessor(one_ptm_filter_mng_ptr));
      one_ptm_filter_processor->process();
      one_ptm_filter_processor = nullptr;
      WebLog::completeFunction(WebLog::OnePtmFilterTime());
      std::cout << "One PTM filtering finished." << std::endl;
      //one ptm filtering results information
      filtering_result_exts.push_back("ONE_PTM_FILTER");
      total_filtering_result_num += one_ptm_filter_mng_ptr->one_ptm_filter_result_num_;
      time(&stop_s);
      std::cout <<  "One PTM filtering running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

      if (ptm_num >= 2) {
        time(&start_s);
        std::cout << "Diagonal filtering started." << std::endl;
        DiagFilterMngPtr diag_filter_mng_ptr 
            = DiagFilterMngPtr(new DiagFilterMng(prsm_para_ptr, filtering_result_num, "DIAG_FILTER"));
        DiagFilterProcessorPtr diag_filter_processor = DiagFilterProcessorPtr(new DiagFilterProcessor(diag_filter_mng_ptr));
        diag_filter_processor->process();
        diag_filter_processor = nullptr;
        WebLog::completeFunction(WebLog::DiagFilterTime());
        std::cout << "Diagonal filtering finished." << std::endl;
        filtering_result_exts.push_back("DIAG_FILTER");
        total_filtering_result_num += diag_filter_mng_ptr->ptm_fast_filter_result_num_; 
        time(&stop_s);
        std::cout <<  "Diagonal filtering running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;
      }
      // debug
      // filtering_result_exts.push_back("ONE_PTM_FILTER");
      //filtering_result_exts.push_back("DIAG_FILTER");
      //total_filtering_result_num = 10; 

      time(&start_s);
      std::cout << "Combining simple PRSMs started." << std::endl;
      LOG_DEBUG("filtering result number " << filtering_result_num);
      SimplePrsmStrCombinePtr simple_combine_ptr(
          new SimplePrsmStrCombine(sp_file_name, filtering_result_exts, "FILTER", total_filtering_result_num));
      simple_combine_ptr->process();
      simple_combine_ptr = nullptr;
      std::cout << "Combining simple PRSMs finished." << std::endl;
      time(&stop_s);
      std::cout <<  "Combine simple prsms running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;

      time(&start_s);
      std::cout << "PTM search started." << std::endl;
      PtmMngPtr ptm_mng_ptr = PtmMngPtr(new PtmMng(prsm_para_ptr, n_top, ptm_num,
                                                   max_ptm_mass, "FILTER", "PTM"));
      PtmProcessorPtr ptm_processor = PtmProcessorPtr(new PtmProcessor(ptm_mng_ptr));
      ptm_processor->process();
      ptm_processor = nullptr;
      WebLog::completeFunction(WebLog::PtmTime());
      std::cout << "PTM search finished" << std::endl;
      time(&stop_s);
      std::cout <<  "Ptm search running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;
    }

    time(&start_s);
    std::cout << "Combining PRSMs started." << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("ZERO_COMPLETE");
    input_exts.push_back("ZERO_PREFIX");
    input_exts.push_back("ZERO_SUFFIX");
    input_exts.push_back("ZERO_INTERNAL");
    if (ptm_num >= 1) {
      input_exts.push_back("PTM");
    }
    int prsm_top_num = (ptm_num + 1) * 4;
    PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num));
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Combining prsms search running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;

    time(&start_s);
    std::cout << "E-value computation started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, ptm_num, max_ptm_mass, use_gf,
                                                      variable_ptm, "RAW_RESULT", "EVALUE"));
    EValueProcessorPtr processor = EValueProcessorPtr(new EValueProcessor(tdgf_mng_ptr));
    processor->init();
    // compute E-value for a set of prsms each run 
    processor->process(false);
    processor = nullptr;
    if (use_gf) {
      WebLog::completeFunction(WebLog::GfEvalueTime());
    }
    else {
      WebLog::completeFunction(WebLog::TableEvalueTime());
    }
    std::cout << "E-value computation finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Computing e-values running time: " << difftime(stop_s, start_s)  << " seconds " << std::endl;

    time(&start_s);
    if (arguments["searchType"]=="TARGET") { 
      std::cout << "Top PRSM selecting started" << std::endl;
      PrsmSelectorPtr selector = PrsmSelectorPtr(new PrsmSelector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top));
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting finished." << std::endl;
    }
    else {
      std::cout << "Top PRSM selecting started " << std::endl;
      PrsmSelectorPtr selector = PrsmSelectorPtr(new PrsmSelector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top));
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting finished." << std::endl;

      std::cout << "FDR computation started. " << std::endl;
      PrsmFdrPtr fdr = PrsmFdrPtr(new PrsmFdr(db_file_name, sp_file_name, "TOP_PRE", "TOP"));
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation finished." << std::endl;
    }
    
    std::cout << "PRSM selecting by cutoff started." << std::endl;
    std::string cutoff_type = arguments["cutoffType"];
    double cutoff_value;
    std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    OutputSelectorPtr output_selector = OutputSelectorPtr(
        new OutputSelector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                           cutoff_type, cutoff_value));
    output_selector->process();
    output_selector = nullptr;
    std::cout << "PRSM selecting by cutoff finished." << std::endl;

    std::cout << "Finding protein species started." << std::endl;
    double ppo;
    std::istringstream (arguments["errorTolerance"]) >> ppo;
    ppo = ppo /1000000.0;
    ResiduePtrVec residue_ptr_vec = prsm_para_ptr->getFixModResiduePtrVec();
    PrsmSpeciesPtr prsm_species = PrsmSpeciesPtr(
        new PrsmSpecies(db_file_name, sp_file_name, "CUTOFF_RESULT", "OUTPUT_RESULT", residue_ptr_vec, ppo));
    prsm_species->process();
    prsm_species = nullptr;
    std::cout << "Finding protein species finished." << std::endl;

    std::cout << "Outputting table started." << std::endl;
    TableWriterPtr table_out = TableWriterPtr(new TableWriter(prsm_para_ptr, "OUTPUT_RESULT", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Postprocessing running time: " << (stop_s-start_s)  << " seconds " << std::endl;

    time(&start_s);
    std::cout << "Generating view xml files started." << std::endl;
    XmlGeneratorPtr xml_gene = XmlGeneratorPtr(new XmlGenerator(prsm_para_ptr, exe_dir,"OUTPUT_RESULT"));
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating view xml files finished." << std::endl;

    std::cout << "Converting xml files to html files started." << std::endl;
    translate(arguments);
    std::cout << "Converting xml files to html files finished." << std::endl;
    time(&stop_s);
    std::cout <<  "Html generation running time: " << difftime(stop_s, start_s) << " seconds " << std::endl;
    */
    
    if (arguments["keepTempFiles"] != "true"){
      std::cout << "Deleting temporary files started." << std::endl;
      delDir(basename(sp_file_name) + "_xml");
      delFile(exe_dir + "/run.log");
      cleanDir(sp_file_name);
      cleanDir(db_file_name);	  
      std::cout << "Deleting temporary files finished." << std::endl;
    }  

    boost::posix_time::ptime stop_time = boost::posix_time::microsec_clock::local_time();

    Summary::outputSummary(argu_processor, start_time, stop_time);
    
    WebLog::close();

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPIC finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  prot::log_level = 2;
  return prot::process(argc, argv);
}

