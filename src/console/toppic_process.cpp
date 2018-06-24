//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"
#include "base/version.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"
#include "spec/feature_util.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_cluster.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_util.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"
#include "zeroptmsearch/zero_ptm_search_mng.hpp"
#include "zeroptmsearch/zero_ptm_search_processor.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "oneptmsearch/one_ptm_search_processor.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "ptmsearch/ptm_search_processor.hpp"

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"

#include "local/local_mng.hpp"
#include "local/local_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/toppic_argument.hpp"

namespace prot {

// protein filtering + database searching + E-value computation
int TopPIC_identify(std::map<std::string, std::string> & arguments) {
  try {
    std::time_t start = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

    arguments["start_time"] = buf;
    Argument::outputArguments(std::cout, arguments);

    std::string resource_dir = arguments["resourceDir"];

    base_data::init(resource_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string feature_file_name = sp_file_name.substr(0, sp_file_name.length() - 12) + ".feature";

    if (arguments["useFeatureFile"] == "true") {
      if (!boost::filesystem::exists(feature_file_name)) {
        LOG_ERROR("TopFD feature file does not exist!. Please use -x option in command line or select 'Missing MS1 feature file in GUI'.");
        exit(EXIT_FAILURE);
      }
    }

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    int ptm_num = std::stoi(arguments["ptmNumber"]);

    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    double min_ptm_mass = std::stod(arguments["minPtmMass"]);

    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
    int thread_num = std::stoi(arguments["threadNumber"]);

    bool use_gf = false;
    if (arguments["useGf"] == "true") {
      use_gf = true;
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
        = std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, thread_num, "ZERO_FILTER");
    ZeroPtmFilterProcessorPtr zero_filter_processor
        = std::make_shared<ZeroPtmFilterProcessor>(zero_filter_mng_ptr);
    zero_filter_processor->process();
    std::cout << "Non PTM filtering - finished." << std::endl;

    std::cout << "Non PTM search - started." << std::endl;
    ZeroPtmSearchMngPtr zero_search_mng_ptr
        = std::make_shared<ZeroPtmSearchMng>(prsm_para_ptr, "ZERO_FILTER", "ZERO_PTM");
    ZeroPtmSearchProcessorPtr zero_search_processor
        = std::make_shared<ZeroPtmSearchProcessor>(zero_search_mng_ptr);
    zero_search_processor->process();
    zero_search_processor = nullptr;
    std::cout << "Non PTM search - finished." << std::endl;

    input_exts.push_back("ZERO_PTM_COMPLETE");
    input_exts.push_back("ZERO_PTM_PREFIX");
    input_exts.push_back("ZERO_PTM_SUFFIX");
    input_exts.push_back("ZERO_PTM_INTERNAL");

    if (ptm_num >= 1) {
      std::cout << "One PTM filtering - started." << std::endl;
      OnePtmFilterMngPtr one_ptm_filter_mng_ptr
          = std::make_shared<OnePtmFilterMng>(prsm_para_ptr, "ONE_PTM_FILTER", thread_num);
      OnePtmFilterProcessorPtr one_filter_processor
          = std::make_shared<OnePtmFilterProcessor>(one_ptm_filter_mng_ptr);
      one_filter_processor->process();
      std::cout << "One PTM filtering - finished." << std::endl;

      std::cout << "One PTM search - started." << std::endl;
      int shift_num = 1;
      PtmSearchMngPtr one_search_mng_ptr
          = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass, min_ptm_mass,
                                           shift_num, thread_num, "ONE_PTM_FILTER", "ONE_PTM");
      OnePtmSearchProcessorPtr one_search_processor
          = std::make_shared<OnePtmSearchProcessor>(one_search_mng_ptr);
      one_search_processor->process();
      one_search_processor = nullptr;
      std::cout << "One PTM search - finished." << std::endl;

      input_exts.push_back("ONE_PTM_COMPLETE");
      input_exts.push_back("ONE_PTM_PREFIX");
      input_exts.push_back("ONE_PTM_SUFFIX");
      input_exts.push_back("ONE_PTM_INTERNAL");
    }

    if (ptm_num >= 2) {
      std::cout << "Diagonal PTM filtering - started." << std::endl;
      DiagFilterMngPtr diag_filter_mng_ptr
          = std::make_shared<DiagFilterMng>(prsm_para_ptr, filter_result_num,
                                            thread_num, "DIAG_FILTER");
      DiagFilterProcessorPtr diag_filter_processor
          = std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      std::cout << "Diagonal filtering - finished." << std::endl;

      std::cout << "Two PTM search - started." << std::endl;
      PtmSearchMngPtr two_search_mng_ptr
          = std::make_shared<PtmSearchMng>(prsm_para_ptr, n_top, max_ptm_mass, min_ptm_mass,
                                           ptm_num, thread_num, "DIAG_FILTER", "PTM");
      PtmSearchProcessorPtr processor = std::make_shared<PtmSearchProcessor>(two_search_mng_ptr);
      processor->process();
      std::cout << "Two PTM search - finished." << std::endl;
      input_exts.push_back("PTM");
    }

    std::cout << "Combining PrSMs - started." << std::endl;
    int prsm_top_num = (ptm_num + 1) * 4;
    PrsmStrCombinePtr combine_ptr
        = std::make_shared<PrsmStrCombine>(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num);
    combine_ptr->process();
    combine_ptr = nullptr;
    std::cout << "Combining PrSMs - finished." << std::endl;

    std::cout << "E-value computation - started." << std::endl;
    bool variable_ptm = false;
    TdgfMngPtr tdgf_mng_ptr
        = std::make_shared<TdgfMng>(prsm_para_ptr, ptm_num,
                                    std::max(std::abs(max_ptm_mass), std::abs(min_ptm_mass)),
                                    use_gf, variable_ptm, thread_num, "RAW_RESULT", "EVALUE");
    EValueProcessorPtr processor = std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return EXIT_SUCCESS;
}

// proteoform clustering + FDR + HTML generation
int TopPIC_post(std::map<std::string, std::string> & arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];

    base_data::init(resource_dir);
    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    double min_ptm_mass = std::stod(arguments["minPtmMass"]);

    int n_top = std::stoi(arguments["numOfTopPrsms"]);

    bool localization = false;
    if (arguments["residueModFileName"] != "") {
      localization = true;
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

    std::cout << "Finding PrSM clusters - started." << std::endl;
    if (arguments["useFeatureFile"] == "true") {
      // TopFD msalign file with feature ID
      double prec_error_tole = 1.2;
      ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
      PrsmFeatureClusterPtr prsm_clusters
          = std::make_shared<PrsmFeatureCluster>(db_file_name,
                                                 sp_file_name,
                                                 "EVALUE",
                                                 "CLUSTERS",
                                                 fix_mod_list,
                                                 prec_error_tole,
                                                 prsm_para_ptr);
      prsm_clusters->process();
      prsm_clusters = nullptr;
    } else {
      double ppo;
      std::istringstream(arguments["errorTolerance"]) >> ppo;
      ppo = ppo / 1000000.0;
      PrsmClusterPtr prsm_clusters
          = std::make_shared<PrsmCluster>(db_file_name, sp_file_name,
                                          "EVALUE", prsm_para_ptr->getFixModPtrVec(),
                                          "CLUSTERS", ppo);
      prsm_clusters->process();
      prsm_clusters = nullptr;
    }
    std::cout << "Finding PrSM clusters - finished." << std::endl;

    if (arguments["searchType"] == "TARGET") {
      std::cout << "Top PrSM selecting - started" << std::endl;
      PrsmTopSelectorPtr selector
          = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name, "CLUSTERS", "TOP", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PrSM selecting - finished." << std::endl;
    } else {
      std::cout << "Top PrSM selecting - started " << std::endl;
      PrsmTopSelectorPtr selector
          = std::make_shared<PrsmTopSelector>(db_file_name, sp_file_name,
                                              "CLUSTERS", "TOP_PRE", n_top);
      selector->process();
      selector = nullptr;
      std::cout << "Top PrSM selecting - finished." << std::endl;

      std::cout << "FDR computation - started. " << std::endl;
      PrsmFdrPtr fdr = std::make_shared<PrsmFdr>(db_file_name, sp_file_name, "TOP_PRE", "TOP");
      fdr->process();
      fdr = nullptr;
      std::cout << "FDR computation - finished." << std::endl;
    }

    std::string cutoff_type = arguments["cutoffSpectralType"];
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    double cutoff_value;
    std::istringstream(arguments["cutoffSpectralValue"]) >> cutoff_value;
    PrsmCutoffSelectorPtr cutoff_selector
        = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP",
                                               "CUTOFF_RESULT_SPEC", cutoff_type, cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    std::string suffix = "CUTOFF_RESULT_SPEC";

    if (localization) {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng
          = std::make_shared<LocalMng>(prsm_para_ptr,
                                       std::stod(arguments["local_threshold"]),
                                       arguments["residueModFileName"],
                                       max_ptm_mass,
                                       min_ptm_mass,
                                       suffix, "LOCAL_RESULT_SPEC");
      LocalProcessorPtr local_ptr = std::make_shared<LocalProcessor>(local_mng);
      local_ptr->process();
      local_ptr = nullptr;
      std::cout << "PTM characterization - finished." << std::endl;
      suffix = "LOCAL_RESULT_SPEC";
    }

    std::time_t end = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["end_time"] = buf;

    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmTableWriterPtr table_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments, suffix, "OUTPUT_TABLE");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    std::cout << "Generating PrSM xml files - started." << std::endl;
    XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, suffix, "prsm_cutoff");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating PrSM xml files - finished." << std::endl;

    std::cout << "Converting PrSM xml files to html files - started." << std::endl;
    translate(arguments, "prsm_cutoff");
    std::cout << "Converting PrSM xml files to html files - finished." << std::endl;

    cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR": "EVALUE";
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    cutoff_selector = std::make_shared<PrsmCutoffSelector>(db_file_name, sp_file_name, "TOP",
                                                           "CUTOFF_RESULT_FORM", cutoff_type,
                                                           cutoff_value);
    cutoff_selector->process();
    cutoff_selector = nullptr;
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    suffix = "CUTOFF_RESULT_FORM";

    if (localization) {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng
          = std::make_shared<LocalMng>(prsm_para_ptr,
                                       std::stod(arguments["local_threshold"]),
                                       arguments["residueModFileName"],
                                       max_ptm_mass,
                                       min_ptm_mass,
                                       suffix, "LOCAL_RESULT_FORM");
      LocalProcessorPtr local_ptr = std::make_shared<LocalProcessor>(local_mng);
      local_ptr->process();
      local_ptr = nullptr;
      std::cout << "PTM characterization - finished." << std::endl;
      suffix = "LOCAL_RESULT_FORM";
    }

    std::cout << "Selecting top PrSMs for proteoforms - started." << std::endl;
    PrsmFormFilterPtr form_filter
        = std::make_shared<PrsmFormFilter>(db_file_name, sp_file_name, suffix,
                                           "FORM_FILTER_RESULT", "FORM_RESULT");
    form_filter->process();
    form_filter = nullptr;
    std::cout << "Selecting top PrSMs for proteoforms - finished." << std::endl;

    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmTableWriterPtr form_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, arguments,
                                            "FORM_RESULT", "FORM_OUTPUT_TABLE");
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

    std::cout << "Generating proteoform xml files - started." << std::endl;
    xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, "CUTOFF_RESULT_FORM", "proteoform_cutoff");
    xml_gene->process();
    xml_gene = nullptr;
    std::cout << "Generating proteoform xml files - finished." << std::endl;

    std::cout << "Converting proteoform xml files to html files - started." << std::endl;
    translate(arguments, "proteoform_cutoff");
    std::cout << "Converting proteoform xml files to html files - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPIC finished." << std::endl;
  return EXIT_SUCCESS;
}

int TopPICProgress(std::map<std::string, std::string> & arguments) {
  TopPIC_identify(arguments);

  TopPIC_post(arguments);

  return EXIT_SUCCESS;
}

int TopPICProgress_multi_file(std::map<std::string, std::string> & arguments,
                              const std::vector<std::string> & spec_file_lst) {
  std::string base_name = arguments["combinedOutputName"];

  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

  std::string start_time_bak = buf;

  std::cout << "TopPIC " << prot::version_number << std::endl;

  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    arguments["spectrumFileName"] = spec_file_lst[k];
    prot::TopPICProgress(arguments);
  }

  if (spec_file_lst.size() > 1) {
    arguments["start_time"] = start_time_bak;
    std::cout << "Merging files - started." << std::endl;
    int N = 100000;
    // merge msalign files
    prot::msalign_util::merge_msalign_files(spec_file_lst, N, base_name + "_ms2.msalign");
    // merge feature files
    std::vector<std::string> feature_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      std::string sp_file_name = spec_file_lst[i];
      feature_file_lst[i] = sp_file_name.substr(0, sp_file_name.length() - 12) + ".feature";
    }
    prot::feature_util::merge_feature_files(feature_file_lst, N, base_name + ".feature");
    // merge EVALUE files
    std::vector<std::string> prsm_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      prsm_file_lst[i] = prot::file_util::basename(spec_file_lst[i]) + ".EVALUE"; 
    }
    prot::prsm_util::merge_prsm_files(prsm_file_lst, N, base_name + "_ms2.EVALUE");
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = base_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;

    prot::TopPIC_post(arguments);
  }

  if (arguments["keepTempFiles"] != "true") {
    std::cout << "Deleting temporary files - started." << std::endl;
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      std::string sp_file_name = spec_file_lst[k];
      prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_proteoform_cutoff_xml");
      prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_prsm_cutoff_xml");
      prot::file_util::cleanDir(ori_db_file_name, sp_file_name);
    }

    std::string sp_file_name = base_name + "_ms2.msalign";
    prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_proteoform_cutoff_xml");
    prot::file_util::delDir(prot::file_util::basename(sp_file_name) + "_prsm_cutoff_xml");
    prot::file_util::cleanDir(ori_db_file_name, sp_file_name);

    std::cout << "Deleting temporary files - finished." << std::endl; 
  }

  return EXIT_SUCCESS;
}

}  // namespace prot

