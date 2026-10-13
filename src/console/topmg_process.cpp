// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include "console/topmg_process.hpp"

#include <algorithm>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include "common/base/base_data.hpp"
#include "common/base/mod_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/version.hpp"
#include "console/topmg_argument.hpp"
#include "filter/diag/diag_filter_processor.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "filter/mng/index_file_name.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "ms/env/env_base.hpp"
#include "ms/feature/feature_merge.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/msalign_util.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_decoy_filter.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_fdr_groups.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_match_num_recount.hpp"
#include "prsm/prsm_match_table_writer.hpp"
#include "prsm/prsm_post_mass_match.hpp"
#include "prsm/prsm_prot_cluster.hpp"
#include "prsm/prsm_prot_filter.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/prsm_sql_output.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_str_merge.hpp"
#include "prsm/simple_prsm_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"
#include "search/graph/graph.hpp"
#include "search/graph/proteo_graph.hpp"
#include "search/graphalign/graph_align_mng.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "search/graphalign/graph_post_processor.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"
#include "sql/sql_db_copy.hpp"
#include "sql/zstd_db.hpp"
#include "stat/mcmc/mcmc_dpr_processor.hpp"
#include "stat/mcmc/mcmc_mng.hpp"

namespace toppic {

namespace {

// The spectrum file whose databases a topmg search uses: the file itself, or
// <x>_ms2.msalign for toppic's post-matched <x>_post_ms2.msalign, so that
// TopFD's database is <x>.sqlite and topmg's <x>_topmg.sqlite
// (prsm_post_mass_match::sqlFileName and topmgSqlFileName).
std::string dbSpFileName(const std::string& sp_file_name) {
  std::string base = file_util::basename(sp_file_name);
  if (str_util::endsWith(base, "_post_ms2")) {
    base = base.substr(0, base.size() - 9) + "_ms2";
  }
  return base + ".msalign";
}

}  // namespace

void cleanTopmgDir(const std::string& fa_name, const std::string& sp_name,
                   bool keep_temp_files) {
  std::string abs_sp_name = file_util::absoluteName(sp_name);
  std::string sp_base = file_util::basename(abs_sp_name);
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');

  file_util::delFile(sp_base + "_topmg_proteoform.xml");
  std::filesystem::rename(sp_base + ".topmg_form_cutoff_form",
                          sp_base + "_topmg_proteoform.xml");

  file_util::delFile(sp_base + "_topmg_protein.xml");
  std::filesystem::rename(sp_base + ".topmg_prot_cutoff_prot",
                          sp_base + "_topmg_protein.xml");

  file_util::delFile(sp_base + "_topmg_prsm_cutoff.xml");
  std::filesystem::rename(sp_base + ".topmg_prsm_cutoff",
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
    file_util::delFile(sp_base + ".topmg_recount");
    file_util::delFile(sp_base + ".topmg_proteoform_cluster");
    file_util::delFile(sp_base + ".topmg_cluster");
    file_util::delFile(sp_base + ".topmg_cluster_fdr");
    file_util::delFile(sp_base + ".topmg_prsm_target");
    file_util::delFile(sp_base + ".topmg_form_cutoff");
    file_util::delFile(sp_base + ".topmg_prot_cutoff");
    file_util::delDir(sp_base + "_topmg_proteoform_cutoff_xml");
    file_util::delDir(sp_base + "_topmg_prsm_cutoff_xml");
  }
}

int TopMG_testModFile(std::map<std::string, std::string>& arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);
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

int TopMG_identify(std::map<std::string, std::string>& arguments) {
  try {
    std::time_t start = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

    arguments["startTime"] = buf;
    // the base data first: the printout takes the masses of C57 and C58 from
    // it
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);
    TopmgArgument::outputArguments(std::cout, " ", arguments);

    EnvBase::initBase(resource_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string var_mod_file_name = arguments["varModFileName"];

    std::string feature_file_name =
        file_util::basename(sp_file_name) + ".feature";

    if (arguments["useFeatureFile"] == "true") {
      if (!std::filesystem::exists(feature_file_name)) {
        LOG_ERROR(
            "TopFD feature file " << feature_file_name
                                  << " does not exist! Please use the -x "
                                     "option in the command line or check "
                                     "'No TopFD feature file' in the GUI.");
        return 1;
      }
    }

    // the ASF filtering keeps this many results per spectrum: TopMG's own
    // number, not the filteringResultNumber of toppic's stages
    int filter_result_num = std::stoi(arguments["asfFilteringResultNumber"]);
    double max_shift_mass = std::stod(arguments["maxShiftMass"]);

    int thread_num = std::stoi(arguments["threadNumber"]);

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
    // the unexpected shifts of the graph alignment: TopMG's own --num-shift,
    // which is not toppic's shiftNumber (PrsmPara::getTopmgShiftNum)
    int shift_num = prsm_para_ptr->getTopmgShiftNum();
    LOG_DEBUG("num of unknown shift " << shift_num);

    // index file name
    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);

    std::vector<std::string> input_exts;

    std::cout << "ASF-One PTM filtering - started." << std::endl;
    // The filtering steps require a large amount of memory, so each one asks
    // mem_check for its own thread number where it starts, from its own
    // per-thread budget and the memory available then.
    int filter_thread_num =
        mem_check::getStageThreadNum("ASF-One PTM filtering", thread_num);
    LOG_DEBUG("Filter thread number " << filter_thread_num);
    int var_ptm_num_in_filtering = 1;
    OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
        std::make_shared<OnePtmFilterMng>(
            prsm_para_ptr, index_file_para, "topmg_one_filter",
            filter_thread_num, -max_shift_mass, max_shift_mass,
            var_mod_file_name, var_ptm_num_in_filtering);
    one_ptm_filter_mng_ptr->inte_num_ = 4;
    one_ptm_filter_mng_ptr->pref_suff_num_ = 4;
    one_ptm_filter_mng_ptr->comp_num_ = 4;
    one_ptm_filter_processor::process(one_ptm_filter_mng_ptr);
    std::cout << "ASF-One PTM filtering - finished." << std::endl;

    input_exts.push_back("topmg_one_filter_COMPLETE");
    input_exts.push_back("topmg_one_filter_PREFIX");
    input_exts.push_back("topmg_one_filter_SUFFIX");
    input_exts.push_back("topmg_one_filter_INTERNAL");

    if (arguments["useAsfDiag"] == "true") {
      std::cout << "ASF-Diagonal PTM filtering - started." << std::endl;
      int diag_filter_thread_num =
          mem_check::getStageThreadNum("ASF-Diagonal PTM filtering",
                                       thread_num);
      LOG_DEBUG("Diagonal filter thread number " << diag_filter_thread_num);
      DiagFilterMngPtr diag_filter_mng_ptr = std::make_shared<DiagFilterMng>(
          prsm_para_ptr, index_file_para, filter_result_num,
          diag_filter_thread_num, "topmg_multi_filter", var_mod_file_name, 1);
      DiagFilterProcessorPtr diag_filter_processor =
          std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      diag_filter_processor = nullptr;
      std::cout << "ASF-Diagonal filtering - finished." << std::endl;

      input_exts.push_back("topmg_multi_filter");
    }

    std::cout << "Combining filtering results - started." << std::endl;
    SimplePrsmStrMergePtr asf_filter_merger =
        std::make_shared<SimplePrsmStrMerge>(sp_file_name, input_exts,
                                             "topmg_graph_filter",
                                             20 * input_exts.size());
    asf_filter_merger->process();
    asf_filter_merger = nullptr;
    std::cout << "Combining filtering results - finished." << std::endl;

    int max_mod_num = std::stoi(arguments["varPtmNumber"]);
    int gap = std::stoi(arguments["proteoGraphGap"]);
    int var_ptm_in_gap =
        std::min(std::stoi(arguments["varPtmNumInGap"]), max_mod_num);
    // whole proteins only is the proteoform type COMPLETE alone
    bool whole_protein_only = (arguments["allowProtType"] == "COMPLETE");
    GraphAlignMngPtr ga_mng_ptr = std::make_shared<GraphAlignMng>(
        prsm_para_ptr, var_mod_file_name, shift_num, max_mod_num, gap,
        var_ptm_in_gap, max_shift_mass,
        std::stod(arguments["topmgAlignErrorTolerance"]), thread_num,
        whole_protein_only, "topmg_graph_filter", "topmg_graph_align");
    std::cout << "Graph alignment - started." << std::endl;
    GraphAlignProcessorPtr ga_processor_ptr =
        std::make_shared<GraphAlignProcessor>(ga_mng_ptr);
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment - finished." << std::endl;

    std::cout << "Graph alignment post-processing - started." << std::endl;
    GraphPostProcessorPtr ga_post_processor_ptr =
        std::make_shared<GraphPostProcessor>(ga_mng_ptr, "topmg_graph_align",
                                             "topmg_graph_post");
    ga_post_processor_ptr->process();
    ga_post_processor_ptr = nullptr;
    std::cout << "Graph alignment post-processing - finished." << std::endl;

    std::cout << "E-value computation using MCMC - started." << std::endl;
    MCMCMngPtr mcmc_mng_ptr = std::make_shared<MCMCMng>(
        prsm_para_ptr, "topmg_graph_post", "topmg_evalue", var_mod_file_name,
        max_mod_num, thread_num);
    DprProcessorPtr processor = std::make_shared<DprProcessor>(mcmc_mng_ptr);
    processor->process();
    processor = nullptr;
    std::cout << "E-value computation using MCMC - finished." << std::endl;

    int n_top = std::stoi(arguments["numOfTopPrsms"]);

    std::cout << "Top PrSM selecting - started" << std::endl;
    prsm_top_selector::process(sp_file_name, "topmg_evalue", "topmg_raw_prsm",
                               n_top);
    std::cout << "Top PrSM selecting - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

int TopMG_post(std::map<std::string, std::string>& arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);
    LOG_DEBUG("Initialization completed");
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name =
        ori_db_file_name + "_idx" + file_util::getFileSeparator() +
        file_util::filenameFromEntirePath(arguments["databaseFileName"]);
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string base_name = file_util::basename(sp_file_name);
    std::string var_mod_file_name = arguments["varModFileName"];
    // the variable modifications, whose residues decide the variable PTM
    // sites of the proteoform classification
    ModPtrVec var_mod_ptrs;
    if (var_mod_file_name != "") {
      var_mod_ptrs = mod_util::readAnywhereModTxt(var_mod_file_name);
    }

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);

    std::cout << "Recounting matched masses and fragments - started."
              << std::endl;
    // the search counted only the masses above the EnvCNN cutoff; recount
    // against the full spectra
    prsm_match_num_recount::process(prsm_para_ptr, "topmg_raw_prsm",
                                    "topmg_recount");
    std::cout << "Recounting matched masses and fragments - finished."
              << std::endl;

    std::cout << "Finding PrSM proteoform clusters - started." << std::endl;
    bool is_proteoform_ppm_error = (arguments["proteoformPpmError"] == "true");
    double proteoform_error_tole =
        std::stod(arguments["proteoformErrorTolerance"]);
    if (is_proteoform_ppm_error) {
      proteoform_error_tole = proteoform_error_tole / 1000000;
    }

    if (arguments["useFeatureFile"] == "true") {
      // TopFD msalign file with feature ID
      ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
      prsm_feature_cluster::process(
          sp_file_name, "topmg_recount", "topmg_proteoform_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    } else {
      prsm_simple_cluster::process(
          db_file_name, sp_file_name, "topmg_recount",
          prsm_para_ptr->getFixModPtrVec(), "topmg_proteoform_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    }
    std::cout << "Finding PrSM proteoform clusters - finished." << std::endl;

    std::cout << "Finding PrSM protein clusters - started." << std::endl;
    prsm_prot_cluster::process(sp_file_name, "topmg_proteoform_cluster",
                               "topmg_cluster");
    std::cout << "Finding PrSM protein clusters - finished." << std::endl;

    std::string cur_suffix = "topmg_cluster";

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      // the FDR step writes the decoys too; they are removed below unless -K
      prsm_fdr::process(sp_file_name, "topmg_cluster", "topmg_cluster_fdr");
      // prsm_fdr_groups::process(sp_file_name, "topmg_cluster",
      // "topmg_cluster_fdr");
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "topmg_cluster_fdr";
      if (arguments["keepDecoyResults"] != "true") {
        std::cout << "Removing decoy identifications - started." << std::endl;
        prsm_decoy_filter::process(base_name, cur_suffix, "topmg_prsm_target");
        cur_suffix = "topmg_prsm_target";
        std::cout << "Removing decoy identifications - finished."
                  << std::endl;
      }
    }

    std::string cutoff_type = arguments["cutoffSpectralType"];
    std::cout << "PrSM filtering by " << cutoff_type << " - started."
              << std::endl;
    double cutoff_value = std::stod(arguments["cutoffSpectralValue"]);
    prsm_cutoff_selector::process(db_file_name, base_name, cur_suffix,
                                  "topmg_prsm_cutoff", cutoff_type,
                                  cutoff_value, var_mod_ptrs);
    std::cout << "PrSM filtering by " << cutoff_type << " - finished."
              << std::endl;
    cur_suffix = "topmg_prsm_cutoff";

    std::time_t end = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["endTime"] = buf;

    std::string argu_str = TopmgArgument::outputTsvArguments(arguments);

    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmMatchTableWriterPtr table_out = std::make_shared<PrsmMatchTableWriter>(
        prsm_para_ptr, argu_str, "topmg_prsm_cutoff");
    table_out->write("_topmg_prsm_single.tsv", false);
    table_out->write("_topmg_prsm.tsv", true);
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    cutoff_type =
        (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR" : "EVALUE";
    std::cout << "Proteoform filtering by " << cutoff_type << " - started."
              << std::endl;
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    prsm_cutoff_selector::process(db_file_name, base_name, cur_suffix,
                                  "topmg_form_cutoff", cutoff_type,
                                  cutoff_value, var_mod_ptrs);
    std::cout << "Proteoform filtering by " << cutoff_type << " - finished."
              << std::endl;

    std::cout << "Selecting top PrSMs for proteoforms - started." << std::endl;
    prsm_form_filter::process(db_file_name, base_name, "topmg_form_cutoff",
                              "topmg_form_cutoff_form");
    std::cout << "Selecting top PrSMs for proteoforms - finished." << std::endl;
    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmMatchTableWriterPtr form_out = std::make_shared<PrsmMatchTableWriter>(
        prsm_para_ptr, argu_str, "topmg_form_cutoff_form");
    form_out->write("_topmg_proteoform_single.tsv", false);
    form_out->write("_topmg_proteoform.tsv", true);
    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

    cutoff_type =
        (arguments["cutoffProteinType"] == "FDR") ? "PROTFDR" : "EVALUE";
    std::cout << "Protein filtering by " << cutoff_type << " - started."
              << std::endl;
    std::istringstream(arguments["cutoffProteinValue"]) >> cutoff_value;
    prsm_cutoff_selector::process(db_file_name, base_name,
                                  "topmg_form_cutoff", "topmg_prot_cutoff",
                                  cutoff_type, cutoff_value, var_mod_ptrs);
    std::cout << "Protein filtering by " << cutoff_type << " - finished."
              << std::endl;

    std::cout << "Selecting top PrSMs for proteins - started." << std::endl;
    prsm_prot_filter::process(db_file_name, base_name, "topmg_prot_cutoff",
                              "topmg_prot_cutoff_prot");
    std::cout << "Selecting top PrSMs for proteins - finished." << std::endl;
    std::cout << "Outputting protein table - started." << std::endl;
    PrsmMatchTableWriterPtr prot_out = std::make_shared<PrsmMatchTableWriter>(
        prsm_para_ptr, argu_str, "topmg_prot_cutoff_prot");
    prot_out->write("_topmg_protein_single.tsv", false);
    prot_out->write("_topmg_protein.tsv", true);
    prot_out = nullptr;
    std::cout << "Outputting protein table - finished." << std::endl;

    // topmg's identifications go into its database of the spectrum file,
    // <x>_topmg.sqlite (TopMGProgress_multi_file; toppic_test experiment
    // 170: into TopFD's database before it); none for the combined result
    // of -c, nor with --no-sql
    if (arguments["noSql"] == "true") {
      std::cout << "The PrSM tables are not written to an SQLite database "
                   "(--no-sql)."
                << std::endl;
    } else {
      prsm_sql_output::write(prsm_para_ptr, arguments["topmgDbFileName"],
                             base_name, ori_db_file_name, "topmg_prsm_cutoff",
                             "topmg_form_cutoff_form",
                             "topmg_prot_cutoff_prot");
    }

  } catch (const char* e) {
    LOG_ERROR("[Exception]" << e);
  }
  return 0;
}

int TopMGProcess(std::map<std::string, std::string>& arguments) {
  if (TopMG_identify(arguments) != 0) {
    return 1;
  }
  return TopMG_post(arguments);
}

int TopMGProgress_multi_file(std::map<std::string, std::string>& arguments,
                             const std::vector<std::string>& spec_file_lst) {
  std::string base_path = file_util::absoluteDir(spec_file_lst[0]);

  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
  std::string combined_start_time = buf;

  std::cout << "TopMG " << Version::getVersion() << std::endl;
  arguments["version"] = Version::getVersion();

  TopMG_testModFile(arguments);

  // check if a combined file name given in -c parameter
  // is the same as one of the input spectrum file. If so, throw error.
  if (arguments["combinedOutputName"] != "") {
    std::string merged_file_name =
        arguments["combinedOutputName"] + "_ms2.msalign";
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      if (merged_file_name == spec_file_lst[k]) {
        std::string raw_file_name =
            spec_file_lst[k].substr(0, spec_file_lst[k].find("_ms2.msalign"));
        LOG_ERROR(
            "A combined file name cannot be the same as one of the input file "
            "names '"
            << raw_file_name
            << "'. Please choose a different name for a combined file and "
               "retry.");
        return 1;
      }
    }
  }

  bool keep_temp_files = (arguments["keepTempFiles"] == "true");
  std::string ori_db_file_name = arguments["oriDatabaseFileName"];

  // each file's database (<x>_topmg.sqlite), compressed at the end
  std::vector<std::string> topmg_db_lst;
  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
    std::string start_time = buf;
    arguments["startTime"] = start_time;
    arguments["spectrumFileName"] = spec_file_lst[k];
    // topmg's database of the file, <x>_topmg.sqlite (toppic_test experiment
    // 170, as toppic's <x>_toppic.sqlite): an earlier one is deleted and,
    // unless --no-sql, a new one written with TopFD's tables ms2_env and
    // ms2_env_peak, into which TopMG_post writes the PrSM tables; TopFD's
    // database is only read, and a reader takes its other tables from it
    std::string db_sp_name = dbSpFileName(spec_file_lst[k]);
    std::string topmg_db_name =
        prsm_post_mass_match::topmgSqlFileName(db_sp_name);
    if (sql_db_copy::removeDb(topmg_db_name)) {
      std::cout << "Deleted the earlier database " << topmg_db_name << "."
                << std::endl;
    }
    if (arguments["noSql"] != "true") {
      sql_db_copy::writeEnvelopeDb(
          prsm_post_mass_match::sqlFileName(db_sp_name), topmg_db_name);
    }
    arguments["topmgDbFileName"] = topmg_db_name;
    topmg_db_lst.push_back(topmg_db_name);
    if (TopMGProcess(arguments) != 0) {
      return 1;
    }
    cleanTopmgDir(ori_db_file_name, spec_file_lst[k], keep_temp_files);
  }

  if (arguments["combinedOutputName"] != "") {
    std::string merged_file_name = arguments["combinedOutputName"];
    std::vector<std::string> raw_file_list;
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      std::string raw_file_name =
          spec_file_lst[k].substr(0, spec_file_lst[k].find("_ms2.msalign"));
      raw_file_list.push_back(raw_file_name);
    }
    std::filesystem::path merge_file_path(merged_file_name);
    if (!merge_file_path.is_absolute()) {
      merged_file_name =
          base_path + file_util::getFileSeparator() + merged_file_name;
    }
    std::string para_str = "";
    std::cout << "Merging files started." << std::endl;
    std::cout << "Merging msalign files started." << std::endl;
    msalign_frac_merge::mergeFractions(raw_file_list, merged_file_name,
                                       para_str);
    std::cout << "Merging msalign files finished." << std::endl;

    if (arguments["useFeatureFile"] ==
        "true") {  // only when feature files are being used
      std::cout << "Merging feature files started." << std::endl;
      feature_merge::process(raw_file_list, merged_file_name, para_str);
      std::cout << "Merging feature files finished." << std::endl;
    }
    // merge TOP files
    std::cout << "Merging identification files started." << std::endl;
    std::vector<std::string> prsm_file_lst(spec_file_lst.size());
    for (size_t i = 0; i < spec_file_lst.size(); i++) {
      prsm_file_lst[i] =
          file_util::basename(spec_file_lst[i]) + ".topmg_raw_prsm";
    }
    prsm_util::mergePrsmFiles(prsm_file_lst, SpPara::getMaxSpecNumPerFile(),
                              SpPara::getMaxFeatureNumPerFile(),
                              merged_file_name + "_ms2.topmg_raw_prsm");
    std::cout << "Merging identification files finished." << std::endl;
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = merged_file_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;
    arguments["startTime"] = combined_start_time;
    // the combined result has no database
    arguments["topmgDbFileName"] = "";
    TopMG_post(arguments);
    cleanTopmgDir(ori_db_file_name, sp_file_name, keep_temp_files);
  }

  // each file's database compressed in place when complete (toppic_test
  // experiment 172), unless --no-sql or the hidden --disable-zstd
  if (arguments["noSql"] != "true" && arguments["disableZstd"] != "true") {
    for (const std::string& name : topmg_db_lst) {
      zstd_db::compressInPlace(name);
    }
  }

  base_data::release();

  std::cout << "TopMG finished." << std::endl << std::flush;
  return 0;
}

}  // namespace toppic
