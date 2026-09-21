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

#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "common/base/base_data.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/ptm_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/version.hpp"
#include "console/toppic_argument.hpp"
#include "filter/diag/diag_filter_processor.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "filter/mng/index_file_name.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "filter/mng/zero_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "filter/varptm/var_ptm_filter_processor.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "ms/env/env_base.hpp"
#include "ms/feature/feature_merge.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/msalign_util.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm_coverage.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_fdr_groups.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_match_num_recount.hpp"
#include "prsm/prsm_match_table_writer.hpp"
#include "prsm/prsm_post_mass_match.hpp"
#include "prsm/prsm_prot_cluster.hpp"
#include "prsm/prsm_prot_filter.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/prsm_sql_output.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "search/oneptmsearch/one_ptm_search_processor.hpp"
#include "search/oneptmsearch/ptm_search_mng.hpp"
#include "search/ptmsearch/ptm_search_processor.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"
#include "search/varptmsearch/var_ptm_search_processor.hpp"
#include "search/zeroptmsearch/zero_ptm_search_mng.hpp"
#include "search/zeroptmsearch/zero_ptm_search_processor.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"
#include "stat/local/local_mng.hpp"
#include "stat/local/local_processor.hpp"
#include "stat/tdgf/evalue_processor.hpp"
#include "stat/tdgf/tdgf_mng.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

// post_sp_name is the spectrum file written by the post mass matching (on by
// default; --disable-post-match turns it off), whose name the result files
// carry; "" if it did not run. toppic_all_raw_prsm (the best PrSM of every
// spectrum before the post-processing) is always kept: a later combined run
// (-c) merges it.
void cleanToppicDir(const std::string& fa_name, const std::string& sp_name,
                    bool keep_temp_files, const std::string& post_sp_name) {
  std::string abs_sp_name = file_util::absoluteName(sp_name);
  std::string sp_base = file_util::basename(abs_sp_name);
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');
  std::string out_base = sp_base;
  if (post_sp_name != "") {
    out_base = file_util::basename(file_util::absoluteName(post_sp_name));
    std::replace(out_base.begin(), out_base.end(), '\\', '/');
  }
  file_util::delFile(out_base + "_toppic_proteoform.xml");
  bool overwrite = true;
  file_util::copyFile(out_base + ".toppic_form_cutoff_form",
                      out_base + "_toppic_proteoform.xml", overwrite);
  file_util::delFile(out_base + "_toppic_prsm.xml");
  file_util::copyFile(out_base + "." + "toppic_prsm_cutoff",
                      out_base + "_toppic_prsm.xml", overwrite);
  file_util::delFile(out_base + "_toppic_protein.xml");
  file_util::copyFile(out_base + ".toppic_prot_cutoff_prot",
                      out_base + "_toppic_protein.xml", overwrite);
  if (!keep_temp_files && post_sp_name != "") {
    file_util::delFile(file_util::absoluteName(post_sp_name) + "_index");
    file_util::delFile(out_base + ".feature");
    file_util::delFile(out_base + ".toppic_proteoform_cluster");
    file_util::delFile(out_base + ".toppic_cluster");
    file_util::delFile(out_base + ".toppic_cluster_fdr");
    file_util::delFile(out_base + ".toppic_cluster_local");
    file_util::delFile(out_base + ".toppic_prsm_cutoff");
    file_util::delFile(out_base + ".toppic_form_cutoff");
    file_util::delFile(out_base + ".toppic_prot_cutoff");
    file_util::delFile(out_base + ".toppic_prot_cutoff_prot");
    file_util::delFile(out_base + ".toppic_form_cutoff_form");
    file_util::delDir(out_base + "_toppic_proteoform_cutoff_xml");
    file_util::delDir(out_base + "_toppic_prsm_cutoff_xml");
  }
  if (!keep_temp_files) {
    file_util::cleanPrefix(sp_name, sp_base + ".msalign_");
    file_util::delFile(abs_sp_name + "_index");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_zero_filter_");
    file_util::delFile(sp_base + ".toppic_zero_shift");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_zero_shift_");

    file_util::cleanPrefix(sp_name, sp_base + ".toppic_var_filter_");
    file_util::delFile(sp_base + ".toppic_var_ptm");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_var_ptm_");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_one_filter_");
    file_util::delFile(sp_base + ".toppic_one_shift");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_one_shift_");
    file_util::delFile(sp_base + ".toppic_multi_filter");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_multi_filter_");
    file_util::delFile(sp_base + ".toppic_multi_combined");
    file_util::delFile(sp_base + ".toppic_multi_evalue");
    file_util::delFile(sp_base + ".toppic_multi_first");
    file_util::delFile(sp_base + ".toppic_multi_skipped");
    file_util::delFile(sp_base + ".toppic_multi_evalue_all");
    file_util::delFile(sp_base + ".toppic_multi_raw_prsm");
    file_util::delFile(sp_base + ".toppic_raw_prsm");
    file_util::delFile(sp_base + ".toppic_multi_shift");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_multi_shift_");
    file_util::delFile(sp_base + ".toppic_combined");
    file_util::delFile(sp_base + ".toppic_evalue");
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_evalue_");
    file_util::delFile(sp_base + ".toppic_recount");
    file_util::delFile(sp_base + ".toppic_proteoform_cluster");
    file_util::delFile(sp_base + ".toppic_cluster");
    file_util::delFile(sp_base + ".toppic_cluster_fdr");
    file_util::delFile(sp_base + ".toppic_cluster_local");
    file_util::delFile(sp_base + ".toppic_prsm_cutoff");
    file_util::delFile(sp_base + ".toppic_form_cutoff");
    file_util::delFile(sp_base + ".toppic_prot_cutoff");
    file_util::delFile(sp_base + ".toppic_prot_cutoff_prot");
    file_util::delFile(sp_base + ".toppic_form_cutoff_form");
    file_util::delDir(sp_base + "_toppic_proteoform_cutoff_xml");
    file_util::delDir(sp_base + "_toppic_prsm_cutoff_xml");
  }
}

// Test modification files.
int TopPIC_testModFile(std::map<std::string, std::string>& arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);
    LOG_DEBUG("Init base data completed");

    // Test arguments
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

    if (arguments["variablePtmFileName"] != "") {
      mod_util::readModTxt(arguments["variablePtmFileName"]);
      ptm_util::readPtmTxt(arguments["variablePtmFileName"]);
    }

    if (arguments["localPtmFileName"] != "") {
      mod_util::readModTxt(arguments["localPtmFileName"]);
      ptm_util::readPtmTxt(arguments["localPtmFileName"]);
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    exit(EXIT_FAILURE);
  }
  return 0;
}

// protein filtering + database searching + E-value computation
int TopPIC_identify(std::map<std::string, std::string>& arguments) {
  try {
    std::time_t start = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

    arguments["startTime"] = buf;
    ToppicArgument::outputArguments(std::cout, " ", arguments);

    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);
    EnvBase::initBase(resource_dir);

    LOG_DEBUG("Init base data completed");

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string feature_file_name =
        file_util::basename(sp_file_name) + ".feature";

    if (arguments["useFeatureFile"] == "true") {
      if (!std::filesystem::exists(feature_file_name)) {
        LOG_ERROR("TopFD feature file does not exist!.");
        LOG_ERROR(
            "Please use -x option in command line or select 'Missing MS1 "
            "feature file in GUI'.");
        return 1;
      }
    }

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    LOG_DEBUG("Number of top PrSMs to retain for each spectrum: " << n_top);
    int var_ptm_num = std::stoi(arguments["variablePtmNum"]);
    int shift_num = std::stoi(arguments["shiftNumber"]);
    std::string var_ptm_file_name = arguments["variablePtmFileName"];

    double max_shift_mass = std::stod(arguments["maxShiftMass"]);
    double min_shift_mass = std::stod(arguments["minShiftMass"]);

    int thread_num = std::stoi(arguments["threadNumber"]);

    // Filter steps requires a large amount of memory.
    // We use only one thread to reduce the memory requirement.
    int filter_thread_num = mem_check::getMaxThreads("zero_one_shift_filter");
    if (filter_thread_num > thread_num) {
      filter_thread_num = thread_num;
    }
    LOG_DEBUG("Filter thread number " << filter_thread_num);

    // var_ptm_type number is used for E-value computation
    int var_ptm_type_num = 0;

    // index file name
    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    LOG_DEBUG("prsm para inited");

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int max_frag_len = std::stoi(arguments["maxFragmentLength"]);
    int min_block_num = std::stoi(arguments["minBlockNum"]);
    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy,
                             db_block_size, max_frag_len, min_block_num);
    msalign_util::geneSpIndex(sp_file_name);

    std::vector<std::string> input_exts;

    std::cout << "Zero unexpected shift filtering - started." << std::endl;
    ZeroPtmFilterMngPtr zero_filter_mng_ptr =
        std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, index_file_para,
                                           filter_thread_num,
                                           "toppic_zero_filter");
    zero_ptm_filter_processor::process(zero_filter_mng_ptr);
    std::cout << "Zero unexpected shift filtering - finished." << std::endl;

    std::cout << "Zero unexpected shift search - started." << std::endl;
    ZeroPtmSearchMngPtr zero_search_mng_ptr =
        std::make_shared<ZeroPtmSearchMng>(prsm_para_ptr, thread_num,
                                           "toppic_zero_filter",
                                           "toppic_zero_shift");
    ZeroPtmSearchProcessorPtr zero_search_processor =
        std::make_shared<ZeroPtmSearchProcessor>(zero_search_mng_ptr);
    zero_search_processor->process();
    zero_search_processor = nullptr;
    std::cout << "Zero unexpected shift search - finished." << std::endl;

    input_exts.push_back("toppic_zero_shift_COMPLETE");
    input_exts.push_back("toppic_zero_shift_PREFIX");
    input_exts.push_back("toppic_zero_shift_SUFFIX");
    input_exts.push_back("toppic_zero_shift_INTERNAL");

    if (var_ptm_num >= 1 && var_ptm_file_name != "") {
      std::cout << "Variable PTM filtering - started." << std::endl;
      bool use_approx_spec = false;
      if (arguments["useApproxSpectra"] == "true") {
        use_approx_spec = true;
      }
      VarPtmFilterMngPtr var_filter_mng_ptr = std::make_shared<VarPtmFilterMng>(
          prsm_para_ptr, index_file_para, var_ptm_file_name, var_ptm_num,
          filter_thread_num, use_approx_spec, "toppic_var_filter");
      var_ptm_type_num = var_filter_mng_ptr->getSingleShiftNum();
      var_ptm_filter_processor::process(var_filter_mng_ptr);
      std::cout << "Variable PTM filtering - finished." << std::endl;

      std::cout << "Variable PTM search - started." << std::endl;
      VarPtmSearchMngPtr var_ptm_search_mng_ptr =
          std::make_shared<VarPtmSearchMng>(
              prsm_para_ptr, n_top, var_ptm_file_name, var_ptm_num, thread_num,
              "toppic_var_filter", "toppic_var_ptm");
      VarPtmSearchProcessorPtr var_ptm_search_processor =
          std::make_shared<VarPtmSearchProcessor>(var_ptm_search_mng_ptr);
      var_ptm_search_processor->process();
      var_ptm_search_processor = nullptr;
      std::cout << "Variable PTM search - finished." << std::endl;
      input_exts.push_back("toppic_var_ptm_COMPLETE");
      input_exts.push_back("toppic_var_ptm_PREFIX");
      input_exts.push_back("toppic_var_ptm_SUFFIX");
      input_exts.push_back("toppic_var_ptm_INTERNAL");
    }

    if (shift_num >= 1) {
      std::cout << "One unexpected shift filtering - started." << std::endl;
      OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
          std::make_shared<OnePtmFilterMng>(
              prsm_para_ptr, index_file_para, "toppic_one_filter",
              filter_thread_num, min_shift_mass, max_shift_mass);
      one_ptm_filter_processor::process(one_ptm_filter_mng_ptr);
      std::cout << "One unexpected shift filtering - finished." << std::endl;

      std::cout << "One unexpected shift search - started." << std::endl;
      int search_shift_num = 1;
      PtmSearchMngPtr one_search_mng_ptr = std::make_shared<PtmSearchMng>(
          prsm_para_ptr, n_top, min_shift_mass, max_shift_mass,
          search_shift_num, thread_num, "toppic_one_filter",
          "toppic_one_shift");
      OnePtmSearchProcessorPtr one_search_processor =
          std::make_shared<OnePtmSearchProcessor>(one_search_mng_ptr);
      one_search_processor->process();
      one_search_processor = nullptr;
      std::cout << "One unexpected shift search - finished." << std::endl;
      input_exts.push_back("toppic_one_shift_COMPLETE");
      input_exts.push_back("toppic_one_shift_PREFIX");
      input_exts.push_back("toppic_one_shift_SUFFIX");
      input_exts.push_back("toppic_one_shift_INTERNAL");
    }

    // the multiple unexpected shifts stage runs after the post-processing
    // of these stages, on the spectra they left unidentified
    // (TopPIC_multi_ptm_search); it needs the variable-PTM type count
    arguments["varPtmTypeNum"] = std::to_string(var_ptm_type_num);

    std::cout << "Merging PrSMs - started." << std::endl;
    int prsm_top_num = std::numeric_limits<int>::max();
    bool norm = false;
    bool remove_dup = false;
    PrsmStrMergePtr merge_ptr = std::make_shared<PrsmStrMerge>(
        sp_file_name, input_exts, "toppic_combined", prsm_top_num, norm,
        remove_dup);
    merge_ptr->process();
    merge_ptr = nullptr;
    std::cout << "Merging PrSMs - finished." << std::endl;

    std::cout << "E-value computation - started." << std::endl;
    TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(
        prsm_para_ptr, shift_num,
        std::max(std::abs(max_shift_mass), std::abs(min_shift_mass)),
        var_ptm_type_num, thread_num, "toppic_combined", "toppic_evalue");
    EValueProcessorPtr processor =
        std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run
    processor->process(false);
    processor = nullptr;
    std::cout << "E-value computation - finished." << std::endl;

    std::cout << "Top PrSM selecting - started" << std::endl;
    prsm_top_selector::process(sp_file_name, "toppic_evalue", "toppic_raw_prsm",
                               n_top);
    std::cout << "Top PrSM selecting - finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

// The multiple unexpected shifts stage (-s 2): the diagonal filter and the
// PTM search on the spectra that TopPIC_identify left without a PrSM passing
// the spectrum-level cutoff, then the merge, E-value computation and top
// PrSM selection of the new PrSMs alone.
int TopPIC_multi_ptm_search(std::map<std::string, std::string>& arguments) {
  try {
    int shift_num = std::stoi(arguments["shiftNumber"]);
    if (shift_num < 2) {
      return 0;
    }
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);

    std::string sp_file_name = arguments["spectrumFileName"];
    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    double max_shift_mass = std::stod(arguments["maxShiftMass"]);
    double min_shift_mass = std::stod(arguments["minShiftMass"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
    int thread_num = std::stoi(arguments["threadNumber"]);
    int diag_filter_thread_num = mem_check::getMaxThreads("diag_filter");
    if (diag_filter_thread_num > thread_num) {
      diag_filter_thread_num = thread_num;
    }
    int var_ptm_type_num = std::stoi(arguments["varPtmTypeNum"]);

    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);

    // The spectra the earlier stages identified are not searched again: a
    // spectrum is skipped when its best first-pass PrSM (toppic_raw_prsm,
    // targets and decoys alike) has an E-value strictly below the skip
    // E-value (hidden --multi-skip-evalue, default 0.02, below the smallest
    // decoy E-value seen so far, 0.033 on 480_F3_02): every spectrum a
    // decoy could compete with is searched again, so the decoys of the
    // final set are those of a full search, while the spectra no decoy can
    // reach keep their first-pass PrSM (report 0024).
    std::string base_name = file_util::basename(sp_file_name);
    std::set<int> identified_spec_ids;
    double skip_evalue = std::stod(arguments["multiSkipEvalue"]);
    int skipped_target_num = 0;
    int skipped_decoy_num = 0;
    PrsmStrPtrVec first_ptrs =
        prsm_reader_util::readAllPrsmStrs(base_name + ".toppic_raw_prsm");
    for (size_t i = 0; i < first_ptrs.size(); i++) {
      if (first_ptrs[i]->getEValue() < skip_evalue) {
        int spec_id = first_ptrs[i]->getSpectrumId();
        if (identified_spec_ids.insert(spec_id).second) {
          if (first_ptrs[i]->getSeqName().find("DECOY_") == 0) {
            skipped_decoy_num++;
          } else {
            skipped_target_num++;
          }
        }
      }
    }
    std::cout << "Multiple unexpected shifts stage - "
              << identified_spec_ids.size()
              << " spectra identified by the earlier stages are skipped ("
              << skipped_target_num << " target, " << skipped_decoy_num
              << " decoy; first-pass E-value < " << skip_evalue << ")."
              << std::endl;

    std::vector<std::string> input_exts;
    {
      std::cout << "Multiple unexpected shifts filtering - started."
                << std::endl;
      // thread number is used because diagonal filter uses only one index
      DiagFilterMngPtr diag_filter_mng_ptr = std::make_shared<DiagFilterMng>(
          prsm_para_ptr, index_file_para, filter_result_num,
          diag_filter_thread_num, "toppic_multi_filter");
      diag_filter_mng_ptr->skip_spec_ids_ = identified_spec_ids;
      DiagFilterProcessorPtr diag_filter_processor =
          std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      diag_filter_processor = nullptr;
      std::cout << "Multiple unexpected shifts filtering - finished."
                << std::endl;

      std::cout << "Multiple unexpected shifts search - started." << std::endl;
      PtmSearchMngPtr multi_search_mng_ptr = std::make_shared<PtmSearchMng>(
          prsm_para_ptr, n_top, min_shift_mass, max_shift_mass, shift_num,
          thread_num, "toppic_multi_filter", "toppic_multi_shift");
      PtmSearchProcessorPtr processor =
          std::make_shared<PtmSearchProcessor>(multi_search_mng_ptr);
      processor->process();
      processor = nullptr;
      std::cout << "Multiple unexpected shifts search - finished." << std::endl;

      input_exts.push_back("toppic_multi_shift_COMPLETE_2");
      input_exts.push_back("toppic_multi_shift_PREFIX_2");
      input_exts.push_back("toppic_multi_shift_SUFFIX_2");
      input_exts.push_back("toppic_multi_shift_INTERNAL_2");
    }


    std::cout << "Merging multiple-shift PrSMs - started." << std::endl;
    int prsm_top_num = std::numeric_limits<int>::max();
    bool norm = false;
    bool remove_dup = false;
    PrsmStrMergePtr merge_ptr = std::make_shared<PrsmStrMerge>(
        sp_file_name, input_exts, "toppic_multi_combined", prsm_top_num, norm,
        remove_dup);
    merge_ptr->process();
    merge_ptr = nullptr;
    std::cout << "Merging multiple-shift PrSMs - finished." << std::endl;

    std::cout << "Multiple-shift E-value computation - started." << std::endl;
    TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(
        prsm_para_ptr, shift_num,
        std::max(std::abs(max_shift_mass), std::abs(min_shift_mass)),
        var_ptm_type_num, thread_num, "toppic_multi_combined",
        "toppic_multi_evalue");
    EValueProcessorPtr processor =
        std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    processor->process(false);
    processor = nullptr;
    std::cout << "Multiple-shift E-value computation - finished." << std::endl;

    // The first-pass PrSMs of the searched spectra (their E-values come
    // from the same tables) compete with the multiple-shift PrSMs for the
    // best PrSM of each spectrum, so that a spectrum a first-pass PrSM
    // matched better keeps that PrSM, target or decoy, in the FDR
    // estimation of the multiple-shift pass.
    std::cout << "Multiple-shift top PrSM selecting - started" << std::endl;
    {
      // the first-pass PrSMs of the searched spectra join the multiple-shift
      // PrSMs (toppic_multi_first); those of the skipped spectra go straight
      // to the final set (toppic_multi_skipped, read by TopPIC_multi_ptm_post)
      PrsmXmlWriter writer(base_name + ".toppic_multi_first");
      PrsmXmlWriter skipped_writer(base_name + ".toppic_multi_skipped");
      for (size_t i = 0; i < first_ptrs.size(); i++) {
        if (identified_spec_ids.count(first_ptrs[i]->getSpectrumId()) == 0) {
          writer.write(first_ptrs[i]);
        } else {
          skipped_writer.write(first_ptrs[i]);
        }
      }
      writer.close();
      skipped_writer.close();
      std::vector<std::string> all_exts = {"toppic_multi_evalue",
                                           "toppic_multi_first"};
      PrsmStrMergePtr all_merge_ptr = std::make_shared<PrsmStrMerge>(
          sp_file_name, all_exts, "toppic_multi_evalue_all",
          std::numeric_limits<int>::max(), false, false);
      all_merge_ptr->process();
      all_merge_ptr = nullptr;
    }
    prsm_top_selector::process(sp_file_name, "toppic_multi_evalue_all",
                               "toppic_multi_raw_prsm", n_top);
    std::cout << "Multiple-shift top PrSM selecting - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

// Post-processing of the identifications. With -s 2 the set is the best PrSM
// of every spectrum: the first-pass PrSM of the skipped spectra and, for the
// searched spectra, the best of their first-pass and multiple-shift PrSMs
// (toppic_all_raw_prsm); without -s 2 (and for the merged file of a -c run)
// it is toppic_raw_prsm. The set is recounted, clustered, characterized and
// given one FDR and the spectrum-level cutoff; then the post mass matching,
// the proteoform and protein cutoffs, the tables and the SQLite file.
int TopPIC_multi_ptm_post(std::map<std::string, std::string>& arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);

    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name =
        ori_db_file_name + "_idx" + file_util::getFileSeparator() +
        file_util::filenameFromEntirePath(arguments["databaseFileName"]);
    double max_shift_mass = std::stod(arguments["maxShiftMass"]);
    double min_shift_mass = std::stod(arguments["minShiftMass"]);
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);
    std::string base_name = file_util::basename(sp_file_name);

    // toppic_all_raw_prsm, the best PrSM of every spectrum before the
    // post-processing, is the file a later combined run (-c) merges; it is
    // the only PrSM file kept without -k (toppic_raw_prsm is deleted)
    if (std::filesystem::exists(base_name + ".toppic_multi_raw_prsm")) {
      std::cout << "Merging the first-pass and multiple-shift PrSMs - started."
                << std::endl;
      std::vector<std::string> raw_exts = {"toppic_multi_skipped",
                                           "toppic_multi_raw_prsm"};
      PrsmStrMergePtr merge_ptr = std::make_shared<PrsmStrMerge>(
          sp_file_name, raw_exts, "toppic_all_raw_prsm",
          std::numeric_limits<int>::max(), false, false);
      merge_ptr->process();
      merge_ptr = nullptr;
      std::cout << "Merging the first-pass and multiple-shift PrSMs - finished."
                << std::endl;
    } else {
      file_util::copyFile(base_name + ".toppic_raw_prsm",
                          base_name + ".toppic_all_raw_prsm", true);
    }

    std::cout << "Recounting matched masses and fragments - started."
              << std::endl;
    // the search counted only the masses above the EnvCNN cutoff; recount
    // against the full spectra
    prsm_match_num_recount::process(prsm_para_ptr, "toppic_all_raw_prsm",
                                    "toppic_recount");
    std::cout << "Recounting matched masses and fragments - finished."
              << std::endl;

    bool is_proteoform_ppm_error = (arguments["proteoformPpmError"] == "true");
    double proteoform_error_tole =
        std::stod(arguments["proteoformErrorTolerance"]);
    if (is_proteoform_ppm_error) {
      proteoform_error_tole = proteoform_error_tole / 1000000;
    }
    std::cout << "Finding PrSM proteoform clusters - started." << std::endl;
    if (arguments["useFeatureFile"] == "true") {
      // TopFD msalign file with feature ID
      prsm_feature_cluster::process(
          sp_file_name, "toppic_recount", "toppic_proteoform_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    } else {
      prsm_simple_cluster::process(
          db_file_name, sp_file_name, "toppic_recount",
          prsm_para_ptr->getFixModPtrVec(), "toppic_proteoform_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    }
    std::cout << "Finding PrSM proteoform clusters - finished." << std::endl;
    std::cout << "Finding PrSM protein clusters - started." << std::endl;
    prsm_prot_cluster::process(sp_file_name, "toppic_proteoform_cluster",
                               "toppic_cluster");
    std::cout << "Finding PrSM protein clusters - finished." << std::endl;
    std::string cur_suffix = "toppic_cluster";

    if (arguments["localPtmFileName"] != "") {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng = std::make_shared<LocalMng>(
          prsm_para_ptr, std::stod(arguments["localThreshold"]),
          arguments["localPtmFileName"], min_shift_mass, max_shift_mass,
          cur_suffix, "toppic_cluster_local");
      LocalProcessorPtr local_ptr = std::make_shared<LocalProcessor>(local_mng);
      local_ptr->process();
      local_ptr = nullptr;
      std::cout << "PTM characterization - finished." << std::endl;
      cur_suffix = "toppic_cluster_local";
    }

    if (arguments["outputRawPrsmTable"] == "true") {
      std::cout << "Outputting Raw PrSM table - started." << std::endl;
      std::string raw_argu_str = ToppicArgument::outputTsvArguments(arguments);
      PrsmMatchTableWriterPtr raw_table_out =
          std::make_shared<PrsmMatchTableWriter>(prsm_para_ptr, raw_argu_str,
                                                 cur_suffix);
      raw_table_out->write("_toppic_raw_prsm_single.tsv", false);
      raw_table_out->write("_toppic_raw_prsm.tsv", true);
      raw_table_out = nullptr;
      std::cout << "Outputting Raw PrSM table - finished." << std::endl;
    }

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      prsm_fdr::process(sp_file_name, cur_suffix, "toppic_cluster_fdr",
                        arguments["keepDecoyResults"]);
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "toppic_cluster_fdr";
    }

    std::string spectra_cutoff_type = arguments["cutoffSpectralType"];
    double spectra_cutoff_value;
    std::istringstream(arguments["cutoffSpectralValue"]) >> spectra_cutoff_value;
    std::cout << "PrSM filtering by " << spectra_cutoff_type << " - started."
              << std::endl;
    prsm_cutoff_selector::process(db_file_name, sp_file_name, cur_suffix,
                                  "toppic_prsm_cutoff", spectra_cutoff_type,
                                  spectra_cutoff_value);
    std::cout << "PrSM filtering by " << spectra_cutoff_type << " - finished."
              << std::endl;

    if (arguments["postMassMatch"] == "true") {
      // On the reported PrSMs only: the FDRs and the spectrum-level cutoff
      // use the E-values, which the added masses do not change.
      std::cout << "Post mass matching - started." << std::endl;
      int min_peak_num = std::stoi(arguments["postMinPeakNum"]);
      // the EnvCNN model scores the added masses
      onnx_env_cnn::initModel(resource_dir,
                              std::stoi(arguments["threadNumber"]));
      sp_file_name = prsm_post_mass_match::process(
          prsm_para_ptr, "toppic_prsm_cutoff", "toppic_prsm_cutoff",
          min_peak_num);
      // The steps below read the spectra with the added masses and the PrSMs
      // with the recounted matches, and name their files after the new
      // spectrum file.
      arguments["spectrumFileName"] = sp_file_name;
      prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
      msalign_util::geneSpIndex(sp_file_name);
      std::cout << "Post mass matching - finished." << std::endl;
    }

    {
      double cutoff_value;
      std::string form_cutoff_type =
          (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR" : "EVALUE";
      std::cout << "Proteoform filtering by " << form_cutoff_type
                << " - started." << std::endl;
      std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
      prsm_cutoff_selector::process(db_file_name, sp_file_name,
                                    "toppic_prsm_cutoff", "toppic_form_cutoff",
                                    form_cutoff_type, cutoff_value);
      std::cout << "Proteoform filtering by " << form_cutoff_type
                << " - finished." << std::endl;
      std::cout << "Selecting top PrSMs for proteoforms - started."
                << std::endl;
      prsm_form_filter::process(db_file_name, sp_file_name,
                                "toppic_form_cutoff",
                                "toppic_form_cutoff_form");
      std::cout << "Selecting top PrSMs for proteoforms - finished."
                << std::endl;

      std::string prot_cutoff_type =
          (arguments["cutoffProteinType"] == "FDR") ? "PROTFDR" : "EVALUE";
      std::cout << "Protein filtering by " << prot_cutoff_type << " - started."
                << std::endl;
      std::istringstream(arguments["cutoffProteinValue"]) >> cutoff_value;
      prsm_cutoff_selector::process(db_file_name, sp_file_name,
                                    "toppic_form_cutoff", "toppic_prot_cutoff",
                                    prot_cutoff_type, cutoff_value);
      std::cout << "Protein filtering by " << prot_cutoff_type << " - finished."
                << std::endl;
      std::cout << "Selecting top PrSMs for proteins - started." << std::endl;
      prsm_prot_filter::process(db_file_name, sp_file_name,
                                "toppic_prot_cutoff",
                                "toppic_prot_cutoff_prot");
      std::cout << "Selecting top PrSMs for proteins - finished." << std::endl;
    }

    std::time_t end = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["endTime"] = buf;
    std::string argu_str = ToppicArgument::outputTsvArguments(arguments);

    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmMatchTableWriterPtr table_out = std::make_shared<PrsmMatchTableWriter>(
        prsm_para_ptr, argu_str, "toppic_prsm_cutoff");
    table_out->write("_toppic_prsm_single.tsv", false);
    table_out->write("_toppic_prsm.tsv", true);
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    if (arguments["outputPrsmCoverage"] == "true") {
      std::cout << "Outputting PrSM coverage - started." << std::endl;
      PrsmCoveragePtr cov_ptr = std::make_shared<PrsmCoverage>(
          prsm_para_ptr, "toppic_prsm_cutoff", "toppic_prsm_coverage");
      cov_ptr->processSingleCoverage();
      cov_ptr = nullptr;
      std::cout << "Outputting PrSM coverage - finished." << std::endl;
    }

    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmMatchTableWriterPtr form_out = std::make_shared<PrsmMatchTableWriter>(
        prsm_para_ptr, argu_str, "toppic_form_cutoff_form");
    form_out->write("_toppic_proteoform_single.tsv", false);
    form_out->write("_toppic_proteoform.tsv", true);
    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

    std::cout << "Outputting protein table - started." << std::endl;
    PrsmMatchTableWriterPtr prot_out = std::make_shared<PrsmMatchTableWriter>(
        prsm_para_ptr, argu_str, "toppic_prot_cutoff_prot");
    prot_out->write("_toppic_protein_single.tsv", false);
    prot_out->write("_toppic_protein.tsv", true);
    prot_out = nullptr;
    std::cout << "Outputting protein table - finished." << std::endl;

    prsm_sql_output::write(prsm_para_ptr, sp_file_name, ori_db_file_name,
                           "toppic_prsm_cutoff", "toppic_form_cutoff_form",
                           "toppic_prot_cutoff_prot");

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

int TopPICProgress(std::map<std::string, std::string>& arguments) {
  TopPIC_identify(arguments);
  TopPIC_multi_ptm_search(arguments);
  TopPIC_multi_ptm_post(arguments);
  return 0;
}

int TopPICProgress_multi_file(std::map<std::string, std::string>& arguments,
                              const std::vector<std::string>& spec_file_lst) {
  std::string base_path = file_util::absoluteDir(spec_file_lst[0]);

  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
  std::string combined_start_time = buf;

  std::cout << "TopPIC " << toppic::Version::getVersion() << std::endl;
  arguments["version"] = toppic::Version::getVersion();

  TopPIC_testModFile(arguments);

  // check if a combined file name given in -c parameter is the same as
  // one of the input spectrum file. If so, throw error.
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

  bool combine_result_only = (arguments["combineResultOnly"] == "true");

  if (!combine_result_only) {
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
      std::string start_time = buf;
      arguments["startTime"] = start_time;
      arguments["spectrumFileName"] = spec_file_lst[k];
      if (toppic::TopPICProgress(arguments) != 0) {
        return 1;
      }
      // the post mass matching replaces the spectrum file name
      std::string post_sp_name = "";
      if (arguments["spectrumFileName"] != spec_file_lst[k]) {
        post_sp_name = arguments["spectrumFileName"];
      }
      cleanToppicDir(ori_db_file_name, spec_file_lst[k], keep_temp_files,
                     post_sp_name);
    }
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
    // Post mass matching runs only on the reported PrSMs of a fraction, so
    // the merge takes TopFD's spectra and the PrSMs before the
    // post-processing, <frac>_ms2.toppic_all_raw_prsm (always kept; the
    // best PrSM of every spectrum over both passes with -s 2). The merged
    // file has no SQLite database, so post mass matching is not applied to
    // the merged results.
    std::string ms2_suffix = "_ms2.msalign";
    std::string prsm_ext = "_ms2.toppic_all_raw_prsm";
    std::string para_str = "";
    std::cout << "Merging files started." << std::endl;
    std::cout << "Merging msalign files started." << std::endl;
    msalign_frac_merge::mergeFractions(raw_file_list, merged_file_name,
                                       para_str, ms2_suffix);
    std::cout << "Merging msalign files finished." << std::endl;

    if (arguments["useFeatureFile"] ==
        "true") {  // only when feature files are being used
      std::cout << "Merging feature files started." << std::endl;
      feature_merge::process(raw_file_list, merged_file_name, para_str);
      std::cout << "Merging feature files finished." << std::endl;
    }
    // merge TOP files
    std::cout << "Merging identification files started." << std::endl;
    std::vector<std::string> prsm_file_lst(raw_file_list.size());
    for (size_t i = 0; i < raw_file_list.size(); i++) {
      prsm_file_lst[i] = raw_file_list[i] + prsm_ext;
    }
    prsm_util::mergePrsmFiles(prsm_file_lst, SpPara::getMaxSpecNumPerFile(),
                              SpPara::getMaxFeatureNumPerFile(),
                              merged_file_name + "_ms2.toppic_raw_prsm");
    std::cout << "Merging identification files finished." << std::endl;
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = merged_file_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;
    arguments["startTime"] = combined_start_time;

    // no SQLite database for the merged spectra: no post mass matching
    std::string post_mass_match = arguments["postMassMatch"];
    arguments["postMassMatch"] = "false";
    TopPIC_multi_ptm_post(arguments);
    arguments["postMassMatch"] = post_mass_match;
    sp_file_name = merged_file_name + "_ms2.msalign";
    std::string post_sp_name = "";
    if (arguments["spectrumFileName"] != sp_file_name) {
      post_sp_name = arguments["spectrumFileName"];
    }
    cleanToppicDir(ori_db_file_name, sp_file_name, keep_temp_files,
                   post_sp_name);
  }

  base_data::release();

  std::cout << "TopPIC finished." << std::endl << std::flush;

  return 0;
}

}  // namespace toppic
