// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <map>
#include <string>
#include <algorithm>
#include <vector>

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"
#include "base/file_util.hpp"

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
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "graph/graph.hpp"
#include "graph/proteo_graph.hpp"
#include "graphalign/graph_align_mng.hpp"
#include "graphalign/graph_align_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/topmg_argument.hpp"

namespace prot {

int TopMGProcess(std::map<std::string, std::string> arguments) {
  try {
    std::cout << "TopMG 1.0-beta (" << __DATE__ << ")" << std::endl;

    std::string exe_dir = arguments["executiveDir"];
    time_t start = time(0);
    char buf[50];
    arguments["start_time"] = std::string(ctime_r(&start, buf));
    Argument::outputArguments(std::cout, arguments);

    BaseData::init(exe_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string residue_mod_file_name = arguments["residueModFileName"];

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    int ptm_num = std::stoi(arguments["ptmNumber"]);
    int thread_num = std::stoi(arguments["threadNumber"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);

    bool use_gf = false;
    if (arguments["useGf"] == "true") {
      use_gf = true;
    }

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("Decoy " << decoy);
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    FastaUtil::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    MsAlignUtil::geneSpIndex(sp_file_name, prsm_para_ptr->getSpParaPtr());

    std::vector<std::string> input_exts;

    std::cout << "ASF-One PTM filtering - started." << std::endl;
    OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
        std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "VAR1_ONE_PTM_FILTER",
                                          arguments["residueModFileName"], 1);
    OnePtmFilterProcessorPtr one_filter_processor =
        std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);
    one_filter_processor->process();
    std::cout << "ASF-One PTM filtering - finished." << std::endl;

    input_exts.push_back("VAR1_ONE_PTM_FILTER_COMPLETE");
    input_exts.push_back("VAR1_ONE_PTM_FILTER_PREFIX");
    input_exts.push_back("VAR1_ONE_PTM_FILTER_SUFFIX");
    input_exts.push_back("VAR1_ONE_PTM_FILTER_INTERNAL");

    bool use_asf_diag = false;
    if (arguments["useASFDiag"] == "true") {
      use_asf_diag = true;
    }

    if (use_asf_diag) {
      std::cout << "ASF-Diagonal PTM filtering - started." << std::endl;
      DiagFilterMngPtr diag_filter_mng_ptr1
          = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num,
                                            thread_num, "VAR1_DIAG_FILTER",
                                            arguments["residueModFileName"], 1);
      DiagFilterProcessorPtr diag_filter_processor1
          = std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr1);
      diag_filter_processor1->process();
      std::cout << "ASF-Diagonal filtering - finished." << std::endl;

      input_exts.push_back("VAR1_DIAG_FILTER");
    }

    SimplePrsmPtrVec sim_prsm_ptrs;

    for (size_t i = 0; i < input_exts.size(); i++) {
      std::string input_file_name
          = FileUtil::basename(sp_file_name) + "." + input_exts[i];
      SimplePrsmReaderPtr reader_ptr = std::make_shared<SimplePrsmReader>(input_file_name);
      SimplePrsmPtr str_ptr = reader_ptr->readOnePrsm();
      while (str_ptr != nullptr) {
        sim_prsm_ptrs.push_back(str_ptr);
        str_ptr = reader_ptr->readOnePrsm();
      }
    }

    sim_prsm_ptrs = SimplePrsmUtil::getUniqueMatches(sim_prsm_ptrs);
    std::sort(sim_prsm_ptrs.begin(), sim_prsm_ptrs.end(), SimplePrsm::cmpIdInc);
    SimplePrsmXmlWriter writer(FileUtil::basename(sp_file_name) + "." + "GRAPH_FILTER");
    writer.write(sim_prsm_ptrs);
    writer.close();

    std::cout << "Graph alignment started." << std::endl;

    int max_mod_num = std::stoi(arguments["varPtmNumber"]);
    int gap = std::stoi(arguments["proteo_graph_dis"]);
    int var_ptm_in_gap = std::min(std::stoi(arguments["varPtmNumInGap"]), max_mod_num);
    GraphAlignMngPtr ga_mng_ptr
        = std::make_shared<GraphAlignMng>(prsm_para_ptr,
                                          residue_mod_file_name,
                                          ptm_num, max_mod_num,
                                          gap, var_ptm_in_gap, max_ptm_mass,
                                          thread_num, "GRAPH_FILTER", "GRAPH_ALIGN");
    LOG_DEBUG("shift num " << ptm_num);
    GraphAlignProcessorPtr ga_processor_ptr = std::make_shared<GraphAlignProcessor>(ga_mng_ptr);
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment finished." << std::endl;

    std::cout << "E-value computation - started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr
        = std::make_shared<TdgfMng>(prsm_para_ptr, ptm_num, max_ptm_mass,
                                    use_gf, variable_ptm, thread_num, "GRAPH_ALIGN", "EVALUE");
    EValueProcessorPtr processor = std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation - finished." << std::endl;

    std::cout << "Finding protein species - started." << std::endl;
    if (arguments["featureFileName"] != "") {
      // TopFD msalign file with feature ID
      double prec_error_tole = 1.2;
      ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
      PrsmFeatureSpeciesPtr prsm_forms
          = std::make_shared<PrsmFeatureSpecies>(db_file_name,
                                                 sp_file_name,
                                                 arguments["featureFileName"],
                                                 "EVALUE",
                                                 "FORMS",
                                                 fix_mod_list,
                                                 prec_error_tole,
                                                 prsm_para_ptr);
      prsm_forms->process();
      prsm_forms = nullptr;
    } else {
      double ppo;
      std::istringstream(arguments["errorTolerance"]) >> ppo;
      ppo = ppo / 1000000.0;
      PrsmSpeciesPtr prsm_species
          = std::make_shared<PrsmSpecies>(db_file_name, sp_file_name,
                                          "EVALUE", prsm_para_ptr->getFixModPtrVec(),
                                          "FORMS", ppo);
      prsm_species->process();
      prsm_species = nullptr;
    }
    std::cout << "Finding protein species - finished." << std::endl;

    if (arguments["searchType"] == "TARGET") {
      std::cout << "Top PRSM selecting - started" << std::endl;
      PrsmTopSelectorPtr selector
          = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "FORMS", "TOP", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting - finished." << std::endl;
    } else {
      std::cout << "Top PRSM selecting - started " << std::endl;
      PrsmTopSelectorPtr selector
          = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name,
                                              "FORMS", "TOP_PRE", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PRSM selecting - finished." << std::endl;

      std::cout << "FDR computation - started. " << std::endl;
      PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation - finished." << std::endl;
    }

    std::cout << "PRSM selecting by cutoff - started." << std::endl;
    std::string cutoff_type = arguments["cutoffSpectralType"];
    double cutoff_value = std::stod(arguments["cutoffSpectralValue"]);
    PrsmCutoffSelectorPtr cutoff_selector
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP",
                                               "CUTOFF_RESULT_SPEC", cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PRSM selecting by cutoff - finished." << std::endl;

    time_t end = time(0);
    arguments["end_time"] = std::string(ctime_r(&end, buf));
    arguments["running_time"] = std::to_string(static_cast<int>(difftime(end, start)));

    std::cout << "Outputting table starts " << std::endl;
    PrsmTableWriterPtr table_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, "CUTOFF_RESULT_SPEC", "OUTPUT_TABLE");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting table finished." << std::endl;

    std::cout << "Generating the PRSM xml files - started." << std::endl;
    XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, exe_dir, "CUTOFF_RESULT_SPEC", "prsm_cutoff");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating the PRSM xml files - finished." << std::endl;

    std::cout << "Converting the PRSM xml files to html files - started." << std::endl;
    translate(arguments, "prsm_cutoff");
    std::cout << "Converting the PRSM xml files to html files - finished." << std::endl;

    std::cout << "Proteoform selecting by cutoff - started." << std::endl;
    cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR": "EVALUE";
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    cutoff_selector = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP",
                                                           "CUTOFF_RESULT_FORM", cutoff_type,
                                                           cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "Proteoform selecting by cutoff - finished." << std::endl;

    std::cout << "Proteoform filtering - started." << std::endl;
    PrsmFormFilterPtr form_filter
        = std::make_shared<PrsmFormFilter>(db_file_name, sp_file_name, "CUTOFF_RESULT_FORM",
                                           "FORM_FILTER_RESULT", "FORM_RESULT");
    form_filter->process();
    form_filter = nullptr;
    std::cout << "Proteoform filtering - finished." << std::endl;

    std::cout << "Outputting the proteoform result table - started." << std::endl;
    PrsmTableWriterPtr form_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments,
                                            "FORM_RESULT", "FORM_OUTPUT_TABLE");
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting the proteoform result table - finished." << std::endl;

    std::cout << "Generating the proteoform xml files - started." << std::endl;
    xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, exe_dir, "CUTOFF_RESULT_FORM",
                                              "proteoform_cutoff");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating the proteoform xml files - finished." << std::endl;

    std::cout << "Converting the proteoform xml files to html files - started." << std::endl;
    translate(arguments, "proteoform_cutoff");
    std::cout << "Converting the proteoform xml files to html files - finished." << std::endl;

    if (arguments["keepTempFiles"] != "true") {
      std::cout << "Deleting temporary files - started." << std::endl;
      FileUtil::delDir(FileUtil::basename(sp_file_name) + "_proteoform_cutoff_xml");
      FileUtil::delDir(FileUtil::basename(sp_file_name) + "_prsm_cutoff_xml");
      FileUtil::cleanDir(ori_db_file_name, sp_file_name);
      std::cout << "Deleting temporary files - finished." << std::endl;
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopMG finished." << std::endl;
  return EXIT_SUCCESS;
}

}  // namespace prot
