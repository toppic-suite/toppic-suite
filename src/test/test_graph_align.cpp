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

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "graph/graph.hpp"
#include "graph/proteo_graph.hpp"

#include "graphalign/graph_align_mng.hpp"
#include "graphalign/graph_align_processor.hpp"

#include <boost/timer.hpp>

namespace prot {


int proteoform_graph_test(int argc, char* argv[]) {

  try {
    Argument argu_processor;
    bool success = argu_processor.parse(argc, argv);
    if (!success) {
      return 1;
    }
    std::map<std::string, std::string> arguments = argu_processor.getArguments();
    std::cout << "TopPIC mass graph" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    argu_processor.outputArguments(std::cout, arguments);

    BaseData::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string log_file_name = arguments["logFileName"];
    std::string residue_mod_file_name = arguments["residueModFileName"];

    int ptm_num = std::stoi(arguments["ptmNumber"]);

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("Decoy " << decoy);
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    FastaUtil::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    MsAlignUtil::geneSpIndex(sp_file_name);

    boost::timer t;

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));

    std::cout << "Graph alignment started." << std::endl;

    int max_mod_num = 10;
    GraphAlignMngPtr ga_mng_ptr = GraphAlignMngPtr(new GraphAlignMng(prsm_para_ptr, residue_mod_file_name, ptm_num, max_mod_num,  "GRAPH_ALIGN"));
    LOG_DEBUG("shift num " << ptm_num);
    GraphAlignProcessorPtr ga_processor_ptr = GraphAlignProcessorPtr(new GraphAlignProcessor(ga_mng_ptr));
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment finished." << std::endl;

    std::cout << "Graph alignment running time: " << t.elapsed()  << " seconds " << std::endl;

    std::cout << "Combining PRSMs started." << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("GRAPH_ALIGN");
    int prsm_top_num = 1;
    PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num));
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs finished." << std::endl;

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, arguments, "RAW_RESULT", "OUTPUT_TABLE"));
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
  //prot::log_level = 2;
  return prot::proteoform_graph_test(argc, argv);
}
