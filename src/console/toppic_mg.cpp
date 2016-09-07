#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"

#include "tdgf/evalue_processor.hpp"
#include "tdgf/tdgf_mng.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/prsm_feature_species.hpp"
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
    time_t start = time(0);
    arguments["start_time"] = std::string(ctime(&start));
    argu_processor.outputArguments(std::cout, arguments);

    BaseData::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string log_file_name = arguments["logFileName"];
    std::string residue_mod_file_name = arguments["residueModFileName"];

    int ptm_num = std::stoi(arguments["ptmNumber"]);

    int thread_num = std::stoi(arguments["threadNumber"]);

    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);

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
    int gap = std::stoi(arguments["proteo_graph_dis"]);
    GraphAlignMngPtr ga_mng_ptr = 
        GraphAlignMngPtr(new GraphAlignMng(prsm_para_ptr, 
                                           residue_mod_file_name, 
                                           ptm_num, max_mod_num, 
                                           gap, max_ptm_mass,
                                           thread_num, "GRAPH_ALIGN"));
    //ga_mng_ptr->prec_error_ = 0;
    LOG_DEBUG("shift num " << ptm_num);
    GraphAlignProcessorPtr ga_processor_ptr = GraphAlignProcessorPtr(new GraphAlignProcessor(ga_mng_ptr));
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment finished." << std::endl;

    std::cout << "Graph alignment running time: " << t.elapsed()  << " seconds " << std::endl;

    std::cout << "Combining PRSMs started." << std::endl;
    std::vector<std::string> input_exts ;
    input_exts.push_back("GRAPH_ALIGN");
    for (int i = 1; i <=100; i++) {
      input_exts.push_back("GRAPH_ALIGN_" + std::to_string(i));
    }
    int prsm_top_num = 1;
    PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num));
    bool normalization = true;
    combine_ptr->process(normalization);
    combine_ptr = nullptr;
    std::cout << "Combining PRSMs finished." << std::endl;

    std::cout << "PRSM selecting by cutoff started." << std::endl;
    std::string cutoff_type = "FRAG";
    double cutoff_value = 10;
    //std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector = PrsmCutoffSelectorPtr(
        new PrsmCutoffSelector(db_file_name, sp_file_name, "RAW_RESULT", "CUTOFF_RESULT", 
                               cutoff_type, cutoff_value));
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PRSM selecting by cutoff finished." << std::endl;

    std::cout << "Finding protein species started." << std::endl;
    double ppo;
    std::istringstream(arguments["errorTolerance"]) >> ppo;
    ppo = ppo / 1000000.0;
    ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
    PrsmFeatureSpeciesPtr prsm_forms = PrsmFeatureSpeciesPtr(
        new PrsmFeatureSpecies(db_file_name, sp_file_name, "CUTOFF_RESULT",
                               "FORMS", fix_mod_list));
    prsm_forms->process();
    prsm_forms = nullptr;
    std::cout << "Finding protein species finished." << std::endl;

    time_t end = time(0);
    arguments["end_time"] = std::string(ctime(&end));
    arguments["running_time"] = std::to_string((int)difftime(end, start));

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, arguments, "FORMS", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;

    std::cout << "PRSM proteoform filtering - started." << std::endl;
    PrsmFormFilterPtr form_filter = PrsmFormFilterPtr(
        new PrsmFormFilter(db_file_name, sp_file_name, "FORMS", "FORM_FILTER_RESULT", "FORM_RESULT"));
    form_filter->process();
    form_filter = nullptr;
    std::cout << "PRSM proteoform filtering - finished." << std::endl;

    std::cout << "Outputting the proteoform table - started." << std::endl;
    PrsmTableWriterPtr form_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, arguments, "FORM_RESULT", "FORM_OUTPUT_TABLE"));
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting the proteoform table - finished." << std::endl;
    /*std::cout << "Generating view xml files started." << std::endl;*/
    //XmlGeneratorPtr xml_gene = XmlGeneratorPtr(new XmlGenerator(prsm_para_ptr, exe_dir, "OUTPUT_RESULT"));
    //xml_gene->process();
    //xml_gene = nullptr;
    //std::cout << "Generating view xml files finished." << std::endl;

    //std::cout << "Converting xml files to html files started." << std::endl;
    //translate(arguments);
    /*std::cout << "Converting xml files to html files finished." << std::endl;*/


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
