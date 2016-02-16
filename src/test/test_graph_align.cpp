#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"


#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/prsm_table_writer.hpp"

#include "console/argument.hpp"
#include "console/summary.hpp"

#include "graph/graph.hpp"
#include "graph/proteo_graph.hpp"

#include "graphalign/graph_align_mng.hpp"
#include "graphalign/graph_align_processor.hpp"

namespace prot {


int proteoform_graph_test(int argc, char* argv[]) {

  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC 0.9.3" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    std::cout << "Executive file directory is: " << exe_dir << std::endl;
    BaseData::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string log_file_name = arguments["logFileName"];
    std::string residue_mod_file_name = arguments["residueModFileName"];

    int ptm_num = std::stoi(arguments["ptmNumber"]);

    /*
       bool use_gf = false; 
       if (arguments["useGf"] == "true") {
       use_gf = true;
       }
       */

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("Decoy " << decoy);
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    FastaUtil::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    MsAlignUtil::geneSpIndex(sp_file_name);

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));

    std::cout << "Graph alignment started." << std::endl;

    int max_mod_num = 10;
    GraphAlignMngPtr ga_mng_ptr = GraphAlignMngPtr(new GraphAlignMng(prsm_para_ptr, residue_mod_file_name, ptm_num, max_mod_num,  "GRAPH_ALIGN"));
    LOG_DEBUG("shift num " << ptm_num);
    GraphAlignProcessorPtr ga_processor_ptr = GraphAlignProcessorPtr(new GraphAlignProcessor(ga_mng_ptr));
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment finished." << std::endl;

    std::cout << "Combining PRSMs started." << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("GRAPH_ALIGN");
    int prsm_top_num = 1;
    PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num));
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs finished." << std::endl;

    /*
       std::cout << "E-value computation started." << std::endl;
       bool variable_ptm = true;
       TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, ptm_num, max_ptm_mass, use_gf,
       variable_ptm, "GRAPH_ALIGN", "EVALUE"));
       EValueProcessorPtr processor = EValueProcessorPtr(new EValueProcessor(tdgf_mng_ptr));
       processor->init();
       LOG_DEBUG("e value init complete");
    // compute E-value for a set of prsms each run 
    processor->process(false);
    processor = nullptr;

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
    */

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, "RAW_RESULT", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "Proteoform test finished." << std::endl;
  return 0;
}

}

int main(int argc, char* argv[]) {
  prot::log_level = 2;
  return prot::proteoform_graph_test(argc, argv);
}
