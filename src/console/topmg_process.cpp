//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
#include <map>
#include <string>
#include <algorithm>
#include <vector>

#include "common/base/mod_util.hpp"
#include "common/base/base_data.hpp"

#include "common/util/file_util.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/version.hpp"
#include "common/util/console_util.hpp"

#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"

#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/deconv_json_merge.hpp"
#include "ms/env/env_base.hpp"
#include "ms/feature/feature_merge.hpp"

#include "para/prsm_para.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_fdr_groups.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_match_table_writer.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_str_merge.hpp"

#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/mng/index_file_name.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"

#include "filter/mng/diag_filter_mng.hpp"
#include "filter/diag/diag_filter_processor.hpp"

#include "search/graph/graph.hpp"
#include "search/graph/proteo_graph.hpp"
#include "search/graphalign/graph_align_mng.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "search/graphalign/graph_post_processor.hpp"

#include "stat/mcmc/mcmc_mng.hpp"
#include "stat/mcmc/mcmc_dpr_processor.hpp"

#include "visual/xml_generator.hpp"
#include "visual/json_transformer.hpp"

#include "console/topmg_argument.hpp"
#include "console/topmg_process.hpp"

namespace toppic {

void cleanTopmgDir(const std::string &fa_name, 
                   const std::string & sp_name,
                   bool keep_temp_files) {
  std::string abs_sp_name = file_util::absoluteName(sp_name);
  std::string sp_base = file_util::basename(abs_sp_name); 
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');

  file_util::delFile(sp_base + "_topmg_proteoform.xml");
  file_util::rename(sp_base + ".topmg_form_cutoff_form", 
                    sp_base + "_topmg_proteoform.xml");

  file_util::delFile(sp_base + "_topmg_prsm_cutoff.xml");
  file_util::rename(sp_base + ".topmg_prsm_cutoff",
                    sp_base + "_topmg_prsm.xml");
  if (!keep_temp_files) {
    file_util::cleanPrefix(sp_name, sp_base + ".msalign_");
    file_util::delFile(abs_sp_name + "_index");
    file_util::delFile(sp_base + ".topmg_one_filter");
    file_util::cleanPrefix(sp_name, sp_base + ".topmg_one_filter_");
    file_util::delFile(sp_base + ".topmg_multi_filter");
    file_util::cleanPrefix(sp_name, sp_base + ".topmg_multi_filter_");
    file_util::delFile(sp_base + ".topmg_graph_filter");
    file_util::cleanPrefix(sp_name, sp_base + ".topmg_graph_filter_");
    file_util::delFile(sp_base + ".topmg_graph_align");
    file_util::cleanPrefix(sp_name, sp_base + ".topmg_graph_align_");
    file_util::delFile(sp_base + ".topmg_graph_post");
    file_util::delFile(sp_base + ".topmg_graph");
    file_util::delFile(sp_base + ".topmg_evalue");
    file_util::cleanPrefix(sp_name, sp_base + ".topmg_evalue_");
    file_util::delFile(sp_base + ".topmg_cluster");
    file_util::delFile(sp_base + ".topmg_cluster_fdr");
    file_util::delFile(sp_base + ".topmg_form_cutoff");
    file_util::delDir(sp_base + "_topmg_proteoform_cutoff_xml");
    file_util::delDir(sp_base + "_topmg_prsm_cutoff_xml");
  }
}

int TopMG_testModFile(std::map<std::string, std::string> & arguments) {
  try {
    base_data::init();
    LOG_DEBUG("Init base data completed");

    // Test arguments
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    if (arguments["varModFileName"] != "") {
      mod_util::readModTxt(arguments["varModFileName"]);
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}

int TopMG_identify(std::map<std::string, std::string> & arguments) {
  try {
    std::time_t start = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

    arguments["startTime"] = buf;
    TopmgArgument::outputArguments(std::cout, " ", arguments);

    std::string resource_dir = arguments["resourceDir"];
    base_data::init();
    EnvBase::initBase(resource_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string var_mod_file_name = arguments["varModFileName"];

    std::string feature_file_name = file_util::basename(sp_file_name) + ".feature";

    if (arguments["useFeatureFile"] == "true") {
      if (!file_util::exists(feature_file_name)) {
        LOG_ERROR("TopFD feature file does not exist!. Please use -x option in command line or select 'Missing MS1 feature file in GUI'.");
        return 1;
      }
    }

    int shift_num = std::stoi(arguments["shiftNumber"]);
    LOG_DEBUG("num of unknown shfit " << shift_num);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
    double max_shift_mass = std::stod(arguments["maxShiftMass"]);

    int thread_num = std::stoi(arguments["threadNumber"]);
    // Filter steps requires a large amount of memory. 
    // We use only one thread to reduce the memory requirement.
    int filter_thread_num = mem_check::getMaxThreads("zero_one_shift_filter");
    if (filter_thread_num > thread_num) {
      filter_thread_num = thread_num;
    }
    LOG_DEBUG("Filter thread number " << filter_thread_num);
    int diag_filter_thread_num = mem_check::getMaxThreads("diag_filter");
    if (diag_filter_thread_num > thread_num) {
      diag_filter_thread_num = thread_num;
    }

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("Decoy " << decoy);
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int max_frag_len = std::stoi(arguments["maxFragmentLength"]);
    int min_block_num = std::stoi(arguments["minBlockNum"]);
    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, 
                             db_block_size, max_frag_len, min_block_num);
    msalign_util::geneSpIndex(sp_file_name);

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    // index file name
    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);


    std::vector<std::string> input_exts;

    std::cout << "ASF-One PTM filtering - started." << std::endl;
    int var_ptm_num_in_filtering = 1;
    OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
      std::make_shared<OnePtmFilterMng>(prsm_para_ptr, index_file_para, 
                                        "topmg_one_filter", 
                                        filter_thread_num,
                                        -max_shift_mass, max_shift_mass, 
                                        var_mod_file_name, 
                                        var_ptm_num_in_filtering);
    one_ptm_filter_mng_ptr->inte_num_ = 4;
    one_ptm_filter_mng_ptr->pref_suff_num_ = 4;
    one_ptm_filter_mng_ptr->comp_num_ = 4;
    one_ptm_filter_processor::process(one_ptm_filter_mng_ptr);
    std::cout << "ASF-One PTM filtering - finished." << std::endl;

    input_exts.push_back("topmg_one_filter_complete");
    input_exts.push_back("topmg_one_filter_prefix");
    input_exts.push_back("topmg_one_filter_suffix");
    input_exts.push_back("topmg_one_filter_internal");


    if (arguments["useAsfDiag"] == "true") {
      std::cout << "ASF-Diagonal PTM filtering - started." << std::endl;
      filter_result_num = 15;
      DiagFilterMngPtr diag_filter_mng_ptr
          = std::make_shared<DiagFilterMng>(prsm_para_ptr, index_file_para, 
                                            filter_result_num,
                                            diag_filter_thread_num, 
                                            "topmg_multi_filter",
                                            var_mod_file_name, 1);
      DiagFilterProcessorPtr diag_filter_processor
          = std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      diag_filter_processor = nullptr;
      std::cout << "ASF-Diagonal filtering - finished." << std::endl;

      input_exts.push_back("topmg_multi_filter");
    }

    std::cout << "Combining filtering results - started." << std::endl;
    SimplePrsmStrMergePtr asf_filter_merger
        = std::make_shared<SimplePrsmStrMerge>(sp_file_name, 
                                               input_exts,
                                               "topmg_graph_filter", 20 * input_exts.size());
    asf_filter_merger->process();
    asf_filter_merger = nullptr;
    std::cout << "Combining filtering results - finished." << std::endl;

    int max_mod_num = std::stoi(arguments["varPtmNumber"]);
    int gap = std::stoi(arguments["proteoGraphGap"]);
    int var_ptm_in_gap = std::min(std::stoi(arguments["varPtmNumInGap"]), max_mod_num);
    bool whole_protein_only = (arguments["wholeProteinOnly"] == "true");
    GraphAlignMngPtr ga_mng_ptr
        = std::make_shared<GraphAlignMng>(prsm_para_ptr,
                                          var_mod_file_name,
                                          shift_num, max_mod_num,
                                          gap, var_ptm_in_gap, max_shift_mass,
                                          thread_num, whole_protein_only, 
                                          "topmg_graph_filter", "topmg_graph_align");
    std::cout << "Graph alignment - started." << std::endl;
    GraphAlignProcessorPtr ga_processor_ptr = std::make_shared<GraphAlignProcessor>(ga_mng_ptr);
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment - finished." << std::endl;

    std::cout << "Graph alignment post-processing - started." << std::endl;
    GraphPostProcessorPtr ga_post_processor_ptr
        = std::make_shared<GraphPostProcessor>(ga_mng_ptr, "topmg_graph_align", "topmg_graph_post");
    ga_post_processor_ptr->process();
    ga_post_processor_ptr = nullptr;
    std::cout << "Graph alignment post-processing - finished." << std::endl;

    std::cout << "E-value computation using MCMC - started." << std::endl;
    MCMCMngPtr mcmc_mng_ptr   
        = std::make_shared<MCMCMng>(prsm_para_ptr, "topmg_graph_post", "topmg_evalue",   
                                    var_mod_file_name, max_mod_num, thread_num);
    DprProcessorPtr processor = std::make_shared<DprProcessor>(mcmc_mng_ptr);   
    processor->process();    
    processor = nullptr;
    std::cout << "E-value computation using MCMC - finished." << std::endl;

    int n_top = std::stoi(arguments["numOfTopPrsms"]);

    std::cout << "Top PrSM selecting - started" << std::endl;
    prsm_top_selector::process(sp_file_name, "topmg_evalue", "topmg_raw_prsm", n_top);
    std::cout << "Top PrSM selecting - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

int TopMG_post(std::map<std::string, std::string> & arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];

    base_data::init();
    LOG_DEBUG("Initialization completed");
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name = ori_db_file_name + "_idx" + file_util::getFileSeparator() + file_util::filenameFromEntirePath(arguments["databaseFileName"]);    
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string var_mod_file_name = arguments["varModFileName"];

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);

    std::cout << "Finding PrSM clusters - started." << std::endl;
    bool is_proteoform_ppm_error = (arguments["proteoformPpmError"] == "true"); 
    double proteoform_error_tole = std::stod(arguments["proteoformErrorTolerance"]);
    if (is_proteoform_ppm_error) {
      proteoform_error_tole = proteoform_error_tole/1000000;
    }

    if (arguments["useFeatureFile"] == "true") {
      // TopFD msalign file with feature ID
      ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
      prsm_feature_cluster::process(sp_file_name,
                                    "topmg_raw_prsm",
                                    "topmg_cluster",
                                    is_proteoform_ppm_error,
                                    proteoform_error_tole);
    } 
    else {
      prsm_simple_cluster::process(db_file_name, 
                                   sp_file_name,
                                   "topmg_raw_prsm", 
                                   prsm_para_ptr->getFixModPtrVec(),
                                   "topmg_cluster", 
                                   is_proteoform_ppm_error,
                                   proteoform_error_tole);
    }
    std::cout << "Finding PrSM clusters - finished." << std::endl;

    std::string cur_suffix = "topmg_cluster";

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      prsm_fdr::process(sp_file_name, "topmg_cluster", "topmg_cluster_fdr", arguments["keepDecoyResults"]);
      //prsm_fdr_groups::process(sp_file_name, "topmg_cluster", "topmg_cluster_fdr", arguments["keepDecoyResults"]);
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "topmg_cluster_fdr";
    }

    std::string cutoff_type = arguments["cutoffSpectralType"];
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    double cutoff_value = std::stod(arguments["cutoffSpectralValue"]);
    prsm_cutoff_selector::process(db_file_name, sp_file_name, cur_suffix,
                                  "topmg_prsm_cutoff", cutoff_type, cutoff_value);
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    std::time_t end = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["endTime"] = buf;

    std::string argu_str = TopmgArgument::outputTsvArguments(arguments);

    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmMatchTableWriterPtr table_out
        = std::make_shared<PrsmMatchTableWriter>(prsm_para_ptr, argu_str, "topmg_prsm_cutoff", "_topmg_prsm_single.tsv", false);
    table_out->write();

    table_out->setOutputName("_topmg_prsm.tsv");
    table_out->setWriteMultiMatches(true);
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, "topmg_prsm_cutoff", "topmg_prsm_cutoff");
    
    if (arguments["geneHTMLFolder"] == "true"){//only when the parameter is set to true
      std::cout << "Generating PrSM xml files - started." << std::endl;
      xml_gene->process();
      xml_gene = nullptr;
      std::cout << "Generating PrSM xml files - finished." << std::endl;

      console_util::copyTopMSV(arguments);

      std::cout << "Converting PrSM xml files to html files - started." << std::endl;
      jsonTranslate(arguments, "topmg_prsm_cutoff");
      std::cout << "Converting PrSM xml files to html files - finished." << std::endl;  
    }
    
    cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR": "EVALUE";
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    prsm_cutoff_selector::process(db_file_name, sp_file_name, cur_suffix,
                                  "topmg_form_cutoff", cutoff_type, cutoff_value);
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    std::cout << "Selecting top PrSMs for proteoforms - started." << std::endl;
    prsm_form_filter::process(db_file_name, sp_file_name, "topmg_form_cutoff",
                              "topmg_form_cutoff_form");
    std::cout << "Selecting top PrSMs for proteoforms - finished." << std::endl;

    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmMatchTableWriterPtr form_out
        = std::make_shared<PrsmMatchTableWriter>(prsm_para_ptr, argu_str,
                                            "topmg_form_cutoff_form", 
                                            "_topmg_proteoform_single.tsv", false);
    form_out->write();

    form_out->setOutputName("_topmg_proteoform.tsv");
    form_out->setWriteMultiMatches(true);
    form_out->write();

    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

    if (arguments["geneHTMLFolder"] == "true"){//only when the parameter is set to true
      std::cout << "Generating proteoform xml files - started." << std::endl;
      xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, "topmg_form_cutoff", "topmg_proteoform_cutoff");

      xml_gene->process();
      xml_gene = nullptr;
      std::cout << "Generating proteoform xml files - finished." << std::endl;

      std::cout << "Converting proteoform xml files to html files - started." << std::endl;
      jsonTranslate(arguments, "topmg_proteoform_cutoff");
      std::cout << "Converting proteoform xml files to html files - finished." << std::endl;
    }
  } catch (const char* e) {
    LOG_ERROR("[Exception]" << e);
  }
  return 0;
}

int TopMGProcess(std::map<std::string, std::string> & arguments) {
  if (TopMG_identify(arguments) != 0) {
    return 1;
  }
  return TopMG_post(arguments);
}

int TopMGProgress_multi_file(std::map<std::string, std::string> & arguments,
                             const std::vector<std::string> & spec_file_lst) {

  std::string base_path = file_util::absoluteDir(spec_file_lst[0]);
  std::string full_combined_name = base_path + file_util::getFileSeparator() 
      +  arguments["combinedOutputName"];

  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
  std::string combined_start_time = buf;

  std::cout << "TopMG " << Version::getVersion() << std::endl;
  arguments["version"] = Version::getVersion();

  xercesc::XMLPlatformUtils::Initialize(); 
  TopMG_testModFile(arguments);

  //check if a combined file name given in -c parameter 
  //is the same as one of the input spectrum file. If so, throw error.
  if (arguments["combinedOutputName"] != "") {
    std::string merged_file_name = arguments["combinedOutputName"] + "_ms2.msalign"; 
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      if (merged_file_name == spec_file_lst[k]) {
        std::string raw_file_name = spec_file_lst[k].substr(0, spec_file_lst[k].find("_ms2.msalign"));
        LOG_ERROR("A combined file name cannot be the same as one of the input file names '" << raw_file_name << "'. Please choose a different name for a combined file and retry.");
        return 1;
      }
    }
  }

  bool keep_temp_files = (arguments["keepTempFiles"] == "true"); 
  std::string ori_db_file_name = arguments["oriDatabaseFileName"];

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
    std::string start_time = buf;
    arguments["startTime"] = start_time;
    arguments["spectrumFileName"] = spec_file_lst[k];
    if (TopMGProcess(arguments) != 0) {
      return 1;
    }
    cleanTopmgDir(ori_db_file_name, spec_file_lst[k], keep_temp_files);
  }

  if (arguments["combinedOutputName"] != "") {
    std::vector<std::string> raw_file_list;
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      std::string raw_file_name = spec_file_lst[k].substr(0, spec_file_lst[k].find("_ms2.msalign"));
      raw_file_list.push_back(raw_file_name);
    }
    std::string para_str = "";
    std::cout << "Merging files started." << std::endl;
    std::cout << "Merging msalign files started." << std::endl;
    msalign_frac_merge::mergeFractions(raw_file_list, full_combined_name, para_str); 
    std::cout << "Merging msalign files finished." << std::endl;

    if (arguments["useFeatureFile"] == "true") {//only when feature files are being used
      std::cout << "Merging feature files started." << std::endl;
      feature_merge::process(raw_file_list, full_combined_name, para_str);
      std::cout << "Merging feature files finished." << std::endl;
    }
    // merge TOP files
    std::cout << "Merging identification files started." << std::endl;
    std::vector<std::string> prsm_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      prsm_file_lst[i] = file_util::basename(spec_file_lst[i]) + ".topmg_raw_prsm"; 
    }
    prsm_util::mergePrsmFiles(prsm_file_lst, SpPara::getMaxSpecNumPerFile(), 
                              SpPara::getMaxFeatureNumPerFile(),
                              full_combined_name + "_ms2.topmg_raw_prsm");
    std::cout << "Merging identification files finished." << std::endl;
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = full_combined_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;
    arguments["startTime"] = combined_start_time;
    // do not generate html files for combined file
    arguments["geneHTMLFolder"] = "false";
    TopMG_post(arguments);
    cleanTopmgDir(ori_db_file_name, sp_file_name, keep_temp_files);
  }


  base_data::release();

  std::cout << "TopMG finished." << std::endl << std::flush;
  return 0; 
}

}  // namespace toppic
