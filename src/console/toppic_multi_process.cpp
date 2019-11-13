//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include <ctime>

#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"
#include "common/base/base_data.hpp"
#include "common/util/version.hpp"

#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/msalign_frac_merge.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_util.hpp"

#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "search/zeroptmsearch/zero_ptm_search_mng.hpp"
#include "search/zeroptmsearch/zero_ptm_search_processor.hpp"

#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "search/oneptmsearch/ptm_search_mng.hpp"
#include "search/oneptmsearch/one_ptm_search_processor.hpp"

#include "filter/diag/diag_filter_mng.hpp"
#include "filter/diag/diag_filter_processor.hpp"
#include "search/ptmsearch/ptm_search_processor.hpp"

#include "stat/tdgf/tdgf_mng.hpp"
#include "stat/tdgf/evalue_processor.hpp"

#include "stat/local/local_mng.hpp"
#include "stat/local/local_processor.hpp"

#include "visual/xml_generator.hpp"
#include "visual/json_transformer.hpp"

#include "console/toppic_argument.hpp"

namespace toppic {

// proteoform clustering + FDR + HTML generation
int TopPIC_post(std::map<std::string, std::string> & arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];

    base_data::init();
    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    int n_top = std::stoi(arguments["numOfTopPrsms"]);

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    LOG_DEBUG("prsm para inited");

    std::time_t start = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
    std::string start_time = buf;
    arguments["startTime"] = start_time;

    /*
    std::cout << "Finding PrSM clusters - started." << std::endl;
    double tolerance = 1.2; 
    PrsmSimpleClusterPtr prsm_clusters
        = std::make_shared<PrsmSimpleCluster>(db_file_name, sp_file_name,
                                              "toppic_top_pre", 
                                              "toppic_cluster", tolerance);
    prsm_clusters->process();
    prsm_clusters = nullptr;
    std::cout << "Finding PrSM clusters - finished." << std::endl;

    std::cout << "Top PrSM selecting - started " << std::endl;
    PrsmTopSelectorPtr selector
        = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name,
                                            "toppic_cluster", "toppic_merge_pre", n_top);
    selector->process();
    selector = nullptr;
    std::cout << "Top PrSM selecting - finished." << std::endl;

    std::cout << "FDR computation - started. " << std::endl;
    PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "toppic_merge_pre", "toppic_top");
    fdr->process();
    fdr = nullptr;
    std::cout << "FDR computation - finished." << std::endl;

    std::string cutoff_type = arguments["cutoffSpectralType"];
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    double cutoff_value;
    std::istringstream(arguments["cutoffSpectralValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "toppic_top",
                                               "toppic_prsm_cutoff", cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    std::string suffix = "toppic_prsm_cutoff";


    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmTableWriterPtr table_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str, suffix, "_toppic_prsm.csv");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR": "EVALUE";
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    cutoff_selector = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, suffix,
                                                           "toppic_form_cutoff", cutoff_type,
                                                           cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    std::cout << "Selecting top PrSMs for proteoforms - started." << std::endl;
    PrsmFormFilterPtr form_filter
        = std::make_shared<PrsmFormFilter>(db_file_name, sp_file_name, "toppic_form_cutoff",
                                           "toppic_form_cutoff_form");
    form_filter->process();
    form_filter = nullptr;
    std::cout << "Selecting top PrSMs for proteoforms - finished." << std::endl;
    */

    std::time_t end = time(nullptr);
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["endTime"] = buf;

    std::string argu_str = Argument::outputCsvArguments(arguments);

    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmTableWriterPtr form_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str,
                                            "toppic_form_cutoff_form", "_toppic_proteoform.csv");
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

int TopPICMultiProcess(std::map<std::string, std::string> & arguments,
                       const std::vector<std::string> & spec_file_lst) {

  std::string base_path = file_util::absoluteDir(spec_file_lst[0]);
  std::string output_name = arguments["combinedOutputName"];

  std::cout << "TopPIC " << toppic::Version::getVersion() << std::endl;

  std::vector<std::string> base_file_lst;
  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    std::string file_name = file_util::basename(spec_file_lst[i]);
    std::string base_name = file_name.substr(0, file_name.length() - 4);
    std::cout << "file name " << base_name << std::endl;
    base_file_lst.push_back(base_name);
  }
  std::cout << "Combined file name " << output_name << std::endl;

  std::cout << "Merging files - started." << std::endl;
  // merge msalign files
  MsAlignFracMerge merger_processor(base_file_lst, output_name);
  std::string para_str = "";
  merger_processor.process(para_str);

  // merge EVALUE files
  std::vector<std::string> prsm_file_lst(spec_file_lst.size());
  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    prsm_file_lst[i] = file_util::basename(spec_file_lst[i]) + ".toppic_top_pre"; 
    std::cout << "prsm file input " << prsm_file_lst[i] << std::endl;
  }
  int N = 1000000;
  prsm_util::mergePrsmFiles(prsm_file_lst, N, output_name + "_ms2.toppic_top_pre");
  std::cout << "Merging files - finished." << std::endl;

  std::string sp_file_name = output_name + "_ms2.msalign";
  arguments["spectrumFileName"] = sp_file_name;
  std::cout << "sp file name " << sp_file_name << std::endl;
  toppic::TopPIC_post(arguments);

  std::cout << "TopPIC Multi fraction finished." << std::endl << std::flush;

  return 0;
}

}  // namespace toppic

