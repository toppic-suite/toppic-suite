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

#include "common/base/base_data.hpp"
#include "common/base/mod_util.hpp"
#include "common/util/version.hpp"

#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"

#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/deconv_json_merge.hpp"
#include "ms/feature/feature_merge.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
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

void copyTopView(std::map<std::string, std::string> &arguments) {
  std::string spectrum_file_name = arguments["spectrumFileName"];
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string base_name_short = base_name.substr(0, base_name.length() - 4);
  std::string topview_dir = base_name_short + "_html" +  file_util::getFileSeparator() + "topview";
  if (file_util::exists(topview_dir)) {
    LOG_WARN("The TopView directory " << topview_dir << " exists!");
    file_util::delDir(topview_dir);
  }

  std::string resource_dir = arguments["resourceDir"];
  // copy resources 
  std::string from_path(resource_dir + file_util::getFileSeparator() + "topview");
  file_util::copyDir(from_path, topview_dir);
}

void cleanToppicDir(const std::string &fa_name, 
                    const std::string & sp_name,
                    bool keep_temp_files) {
  std::string fa_base = file_util::absoluteName(fa_name);
  std::replace(fa_base.begin(), fa_base.end(), '\\', '/');
  std::string abs_sp_name = file_util::absoluteName(sp_name);
  std::string sp_base = file_util::basename(abs_sp_name);
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');

  file_util::rename(sp_base + ".toppic_form_cutoff_form",
                    sp_base + "_toppic_proteoform.xml");
  if (!keep_temp_files) {
    file_util::cleanPrefix(fa_name, fa_base + "_");
    file_util::cleanPrefix(sp_name, sp_base + ".msalign_");
    file_util::delFile(abs_sp_name + "_index");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_zero_filter_");
    file_util::delFile(sp_base + ".toppic_zero_ptm");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_zero_ptm_");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_one_filter_");
    file_util::delFile(sp_base + ".toppic_one_ptm");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_one_ptm_");
    file_util::delFile(sp_base + ".toppic_multi_filter");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_multi_filter_");
    file_util::delFile(sp_base + ".toppic_multi_ptm");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_multi_ptm_");
    file_util::delFile(sp_base + ".toppic_combined");
    file_util::delFile(sp_base + ".toppic_evalue");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_evalue_");
    file_util::delFile(sp_base + ".toppic_top");
    file_util::delFile(sp_base + ".toppic_cluster");
    file_util::delFile(sp_base + ".toppic_cluster_fdr");
    file_util::delFile(sp_base + ".toppic_prsm_cutoff");
    file_util::delFile(sp_base + ".toppic_prsm_cutoff_local");
    file_util::delFile(sp_base + ".toppic_form_cutoff");
    file_util::delDir(sp_base + "_toppic_proteoform_cutoff_xml");
    file_util::delDir(sp_base + "_toppic_prsm_cutoff_xml");
  }
}

// protein filtering + database searching + E-value computation
int TopPIC_testModFile(std::map<std::string, std::string> & arguments) {
  try {
    base_data::init();
    LOG_DEBUG("Init base data completed");

    // Test arguments
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    if (arguments["residueModFileName"] != "") {
      LOG_ERROR("read file");
      mod_util::readModTxt(arguments["residueModFileName"]);
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}

// protein filtering + database searching + E-value computation
int TopPIC_identify(std::map<std::string, std::string> & arguments) {
  try {
    std::time_t start = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

    arguments["startTime"] = buf;
    Argument::outputArguments(std::cout, arguments);

    base_data::init();

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string feature_file_name = file_util::basename(sp_file_name) + ".feature";

    if (arguments["useFeatureFile"] == "true") {
      if (!file_util::exists(feature_file_name)) {
        LOG_ERROR("TopFD feature file does not exist!.");
        LOG_ERROR("Please use -x option in command line or select 'Missing MS1 feature file in GUI'.");
        return 1;
      }
    }

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    int ptm_num = std::stoi(arguments["ptmNumber"]);

    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    double min_ptm_mass = std::stod(arguments["minPtmMass"]);

    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
    int thread_num = std::stoi(arguments["threadNumber"]);

    // Filter steps requires a large amount of memory. 
    // We use only one thread to reduce the memory requirement.
    //int filter_thread_num = 1;

    bool use_gf = true;
    if (arguments["useLookupTable"] == "true") {
      use_gf = false;
    }

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    LOG_DEBUG("prsm para inited");

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);

    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy, db_block_size);
    msalign_util::geneSpIndex(sp_file_name, prsm_para_ptr->getSpParaPtr());

    std::vector<std::string> input_exts;

    std::cout << "Non PTM filtering - started." << std::endl;
    ZeroPtmFilterMngPtr zero_filter_mng_ptr
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, thread_num, "toppic_zero_filter");
    ZeroPtmFilterProcessorPtr zero_filter_processor
        = std::make_shared<ZeroPtmFilterProcessor>(zero_filter_mng_ptr);
    zero_filter_processor->process();
    zero_filter_processor = nullptr;
    std::cout << "Non PTM filtering - finished." << std::endl;

    std::cout << "Non PTM search - started." << std::endl;
    ZeroPtmSearchMngPtr zero_search_mng_ptr
        = std::make_shared<ZeroPtmSearchMng>(prsm_para_ptr, "toppic_zero_filter", "toppic_zero_ptm");
    ZeroPtmSearchProcessorPtr zero_search_processor
        = std::make_shared<ZeroPtmSearchProcessor>(zero_search_mng_ptr);
    zero_search_processor->process();
    zero_search_processor = nullptr;
    std::cout << "Non PTM search - finished." << std::endl;

    input_exts.push_back("toppic_zero_ptm_complete");
    input_exts.push_back("toppic_zero_ptm_prefix");
    input_exts.push_back("toppic_zero_ptm_suffix");
    input_exts.push_back("toppic_zero_ptm_internal");

    if (ptm_num >= 1) {
      std::cout << "One PTM filtering - started." << std::endl;
      OnePtmFilterMngPtr one_ptm_filter_mng_ptr
          = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "toppic_one_filter", thread_num);
      OnePtmFilterProcessorPtr one_filter_processor
          = std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);
      one_filter_processor->process();
      one_filter_processor = nullptr;
      std::cout << "One PTM filtering - finished." << std::endl;

      std::cout << "One PTM search - started." << std::endl;
      int shift_num = 1;
      PtmSearchMngPtr one_search_mng_ptr
          = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass, min_ptm_mass,
                                           shift_num, thread_num, "toppic_one_filter", "toppic_one_ptm");
      OnePtmSearchProcessorPtr one_search_processor
          = std::make_shared<OnePtmSearchProcessor>(one_search_mng_ptr);
      one_search_processor->process();
      one_search_processor = nullptr;
      std::cout << "One PTM search - finished." << std::endl;

      input_exts.push_back("toppic_one_ptm_complete");
      input_exts.push_back("toppic_one_ptm_prefix");
      input_exts.push_back("toppic_one_ptm_suffix");
      input_exts.push_back("toppic_one_ptm_internal");
    }

    if (ptm_num >= 2) {
      std::cout << "Multiple PTM filtering - started." << std::endl;
      // thread number is used because diagonal filter uses only one index
      DiagFilterMngPtr diag_filter_mng_ptr
          = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num,
                                            thread_num, "toppic_multi_filter");
      DiagFilterProcessorPtr diag_filter_processor
          = std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      diag_filter_processor = nullptr;
      std::cout << "Multiple PTM filtering - finished." << std::endl;

      std::cout << "Multiple PTM search - started." << std::endl;
      PtmSearchMngPtr multi_search_mng_ptr
          = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass, min_ptm_mass,
                                           ptm_num, thread_num, "toppic_multi_filter", "toppic_multi_ptm");
      PtmSearchProcessorPtr processor = std::make_shared<PtmSearchProcessor>(multi_search_mng_ptr);
      processor->process();
      processor = nullptr;
      std::cout << "Multiple PTM search - finished." << std::endl;

      input_exts.push_back("toppic_multi_ptm_complete_2");
      input_exts.push_back("toppic_multi_ptm_prefix_2");
      input_exts.push_back("toppic_multi_ptm_suffix_2");
      input_exts.push_back("toppic_multi_ptm_internal_2");
    }

    std::cout << "Merging PrSMs - started." << std::endl;
    int prsm_top_num = (ptm_num + 1) * 4;
    PrsmStrMergePtr merge_ptr
        = std::make_shared<PrsmStrMerge>(sp_file_name, input_exts, "toppic_combined", prsm_top_num);
    merge_ptr->process();
    merge_ptr = nullptr;
    std::cout << "Merging PrSMs - finished." << std::endl;
    
    std::cout << "E-value computation - started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr
        = std::make_shared<TdgfMng>(prsm_para_ptr, ptm_num,
                                    std::max(std::abs(max_ptm_mass), std::abs(min_ptm_mass)),
                                    use_gf, variable_ptm, thread_num, "toppic_combined", "toppic_evalue");
    EValueProcessorPtr processor = std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation - finished." << std::endl;

    std::cout << "Top PrSM selecting - started" << std::endl;
    PrsmTopSelectorPtr selector
        = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, 
                                            "toppic_evalue", "toppic_top", n_top);
    selector->process();
    selector = nullptr;
    std::cout << "Top PrSM selecting - finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

// proteoform clustering + FDR + HTML generation
int TopPIC_post(std::map<std::string, std::string> & arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];

    base_data::init();
    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    double min_ptm_mass = std::stod(arguments["minPtmMass"]);

    bool localization = false;
    if (arguments["residueModFileName"] != "") {
      localization = true;
    }

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name, prsm_para_ptr->getSpParaPtr());
    LOG_DEBUG("prsm para inited");

    std::cout << "Finding PrSM clusters - started." << std::endl;
    double proteoform_error_tole = 1.2;
    if (arguments["useFeatureFile"] == "true") {
      // TopFD msalign file with feature ID
      ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
      PrsmFeatureClusterPtr prsm_clusters
          = std::make_shared<PrsmFeatureCluster>(db_file_name,
                                                 sp_file_name,
                                                 "toppic_top",
                                                 "toppic_cluster",
                                                 fix_mod_list,
                                                 proteoform_error_tole,
                                                 prsm_para_ptr);
      prsm_clusters->process();
      prsm_clusters = nullptr;
    } else {
      PrsmSimpleClusterPtr prsm_clusters
          = std::make_shared<PrsmSimpleCluster>(db_file_name, sp_file_name,
                                                "toppic_top", prsm_para_ptr->getFixModPtrVec(),
                                                "toppic_cluster", proteoform_error_tole);
      prsm_clusters->process();
      prsm_clusters = nullptr;
    }
    std::cout << "Finding PrSM clusters - finished." << std::endl;
    std::string cur_suffix = "toppic_cluster";

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "toppic_cluster", "toppic_cluster_fdr");
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "toppic_cluster_fdr";
    }

    std::string cutoff_type = arguments["cutoffSpectralType"];
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    double cutoff_value;
    std::istringstream(arguments["cutoffSpectralValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, cur_suffix,
                                               "toppic_prsm_cutoff", cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;
    cur_suffix = "toppic_prsm_cutoff";

    if (localization) {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng
          = std::make_shared<LocalMng>(prsm_para_ptr,
                                       std::stod(arguments["localThreshold"]),
                                       arguments["residueModFileName"],
                                       max_ptm_mass,
                                       min_ptm_mass,
                                       "toppic_prsm_cutoff", 
                                       "toppic_prsm_cutoff_local");
      LocalProcessorPtr local_ptr = std::make_shared<LocalProcessor>(local_mng);
      local_ptr->process();
      local_ptr = nullptr;
      std::cout << "PTM characterization - finished." << std::endl;
      cur_suffix = "toppic_prsm_cutoff_local";
    }


    std::time_t end = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["endTime"] = buf;

    std::string argu_str = Argument::outputCsvArguments(arguments);

    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmTableWriterPtr table_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str, cur_suffix, "_toppic_prsm.csv");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    copyTopView(arguments);

    std::cout << "Generating PrSM xml files - started." << std::endl;
    XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, 
                                                              cur_suffix, "toppic_prsm_cutoff");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating PrSM xml files - finished." << std::endl;

    std::cout << "Converting PrSM xml files to json files - started." << std::endl;
    jsonTranslate(arguments, "toppic_prsm_cutoff");
    std::cout << "Converting PrSM xml files to json files - finished." << std::endl;

    cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR": "EVALUE";
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    cutoff_selector = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, cur_suffix,
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

    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmTableWriterPtr form_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str,
                                            "toppic_form_cutoff_form", "_toppic_proteoform.csv");
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

    std::cout << "Generating proteoform xml files - started." << std::endl;
    xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, 
                                              "toppic_form_cutoff", 
                                              "toppic_proteoform_cutoff");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating proteoform xml files - finished." << std::endl;

    std::cout << "Converting proteoform xml files to html files - started." << std::endl;
    jsonTranslate(arguments, "toppic_proteoform_cutoff");
    std::cout << "Converting proteoform xml files to html files - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

int TopPICProgress(std::map<std::string, std::string> & arguments) {
  TopPIC_identify(arguments); 
  TopPIC_post(arguments);
  return 0;
}

int TopPICProgress_multi_file(std::map<std::string, std::string> & arguments,
                              const std::vector<std::string> & spec_file_lst) {

  std::string base_path = file_util::absoluteDir(spec_file_lst[0]);
  std::string base_name = base_path + file_util::getFileSeparator() 
      +  arguments["combinedOutputName"];

  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
  std::string combined_start_time = buf;

  xercesc::XMLPlatformUtils::Initialize(); 

  std::cout << "TopPIC " << toppic::Version::getVersion() << std::endl;

  TopPIC_testModFile(arguments);

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
    std::string start_time = buf;
    arguments["startTime"] = start_time;
    arguments["spectrumFileName"] = spec_file_lst[k];
    if (toppic::TopPICProgress(arguments) != 0) {
      return 1;
    }
  }

  if (spec_file_lst.size() > 1 && arguments["combinedOutputName"] != "") {
    std::string merged_file_name = arguments["combinedOutputName"]; 
    std::string para_str = "";
    std::cout << "Merging files started." << std::endl;
    toppic::MsAlignFracMerge::mergeFiles(spec_file_lst, merged_file_name + "_ms2.msalign", para_str);
    toppic::DeconvJsonMergePtr json_merger 
        = std::make_shared<toppic::DeconvJsonMerge>(spec_file_lst, merged_file_name);
    json_merger->process();
    json_merger = nullptr;
    toppic::FeatureMergePtr feature_merger 
        = std::make_shared<toppic::FeatureMerge>(spec_file_lst, merged_file_name);
    feature_merger->process(para_str);
    feature_merger = nullptr;

    // merge TOP files
    std::vector<std::string> prsm_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      prsm_file_lst[i] = toppic::file_util::basename(spec_file_lst[i]) + ".toppic_top"; 
    }
    int N = 1000000;
    toppic::prsm_util::mergePrsmFiles(prsm_file_lst, N , base_name + "_ms2.toppic_top");
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = base_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;
    arguments["startTime"] = combined_start_time;
    toppic::TopPIC_post(arguments);
  }

  bool keep_temp_files = (arguments["keepTempFiles"] == "true"); 

  std::cout << "Deleting temporary files - started." << std::endl;
  std::string ori_db_file_name = arguments["oriDatabaseFileName"];

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    std::string sp_file_name = spec_file_lst[k];
    cleanToppicDir(ori_db_file_name, sp_file_name, keep_temp_files);
  }

  if (spec_file_lst.size() > 1 && arguments["combinedOutputName"] != "") {
    std::string sp_file_name = base_name + "_ms2.msalign";
    cleanToppicDir(ori_db_file_name, sp_file_name, keep_temp_files);
  }
  std::cout << "Deleting temporary files - finished." << std::endl; 

  std::cout << "TopPIC finished." << std::endl << std::flush;

  return 0;
}

}  // namespace toppic

