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

#include "console/toppic_process.hpp"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "common/base/base_data.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/ptm_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"
#include "console/toppic_argument.hpp"
#include "filter/diag/diag_filter_processor.hpp"
#include "filter/mng/diag_filter_mng.hpp"
#include "filter/mng/index_file_name.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "filter/mng/zero_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_processor.hpp"
#include "filter/tag/tag_filter_merge.hpp"
#include "filter/tag/tag_filter_mng.hpp"
#include "filter/tag/tag_filter_processor.hpp"
#include "filter/varptm/var_ptm_filter_processor.hpp"
#include "filter/varptm/zero_var_ptm_filter_processor.hpp"
#include "filter/zeroptm/zero_ptm_filter_processor.hpp"
#include "ms/calib/mass_calib_fit.hpp"
#include "ms/calib/mass_calib_model.hpp"
#include "ms/calib/msalign_mass_calib.hpp"
#include "ms/calib/sqlite_mass_calib.hpp"
#include "ms/env/avg_mass_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/msalign_util.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm_coverage.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_deamid_selector.hpp"
#include "prsm/prsm_decoy_filter.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_fdr_groups.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_mass_calib.hpp"
#include "prsm/prsm_match_num_recount.hpp"
#include "prsm/prsm_match_table_writer.hpp"
#include "prsm/prsm_post_mass_match.hpp"
#include "prsm/prsm_prot_cluster.hpp"
#include "prsm/prsm_prot_filter.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_second_prec.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/simple_prsm_str_merge.hpp"
#include "prsm/prsm_sql_output.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "search/graphalign/graph_align_mng.hpp"
#include "search/graphalign/graph_align_processor.hpp"
#include "search/graphalign/graph_post_processor.hpp"
#include "search/oneptmsearch/one_ptm_search_processor.hpp"
#include "search/oneptmsearch/ptm_search_mng.hpp"
#include "search/ptmsearch/ptm_search_processor.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"
#include "search/varptmsearch/var_ptm_search_processor.hpp"
#include "search/zeroptmsearch/zero_ptm_search_mng.hpp"
#include "search/zeroptmsearch/zero_ptm_search_processor.hpp"
#include "seq/fasta_reader.hpp"
#include "seq/fasta_util.hpp"
#include "sql/sql_db_copy.hpp"
#include "sql/zstd_db.hpp"
#include "stat/local/local_mng.hpp"
#include "stat/local/local_processor.hpp"
#include "stat/local_align/local_align_processor.hpp"
#include "stat/mcmc/mcmc_dpr_processor.hpp"
#include "stat/mcmc/mcmc_mng.hpp"
#include "stat/rescore/prsm_rescore.hpp"
#include "stat/tdgf/cluster_prot_corrector.hpp"
#include "stat/tdgf/evalue_processor.hpp"
#include "stat/tdgf/tdgf_mng.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

// post_sp_name is the spectrum file written by the post mass matching (on by
// default; --disable-post-match turns it off), whose name the result files
// carry; "" if it did not run. toppic_record_prsm is never deleted.
// copy_outputs: the cutoff files become the XML outputs (false for the
// residual spectra of --second-precursor, which have no report)
void cleanToppicDir(const std::string& fa_name, const std::string& sp_name,
                    bool keep_temp_files, const std::string& post_sp_name,
                    bool copy_outputs = true) {
  std::string abs_sp_name = file_util::absoluteName(sp_name);
  std::string sp_base = file_util::basename(abs_sp_name);
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');
  std::string out_base = sp_base;
  if (post_sp_name != "") {
    out_base = file_util::basename(file_util::absoluteName(post_sp_name));
    std::replace(out_base.begin(), out_base.end(), '\\', '/');
  }
  if (copy_outputs) {
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
  }
  if (!keep_temp_files && post_sp_name != "") {
    file_util::delFile(file_util::absoluteName(post_sp_name) + "_index");
    file_util::delFile(out_base + ".feature");
    file_util::delFile(out_base + ".toppic_recount");
    file_util::delFile(out_base + ".toppic_proteoform_cluster");
    file_util::delFile(out_base + ".toppic_cluster");
    file_util::delFile(out_base + ".toppic_cluster_fdr");
    file_util::delFile(out_base + ".toppic_cluster_local");
    file_util::delFile(out_base + ".toppic_rescore");
    file_util::delFile(out_base + ".toppic_cluster_prot");
    // the second precursor search's merge (--second-precursor)
    file_util::delFile(out_base + ".toppic_first_record_prsm");
    file_util::delFile(out_base + ".toppic_second_merged");
    file_util::delFile(out_base + ".toppic_second_cluster");
    file_util::delFile(out_base + ".toppic_second_prot");
    file_util::delFile(out_base + ".toppic_second_rescore");
    file_util::delFile(out_base + ".toppic_second_fdr");
    file_util::delFile(out_base + ".toppic_prsm_target");
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
    // the tag filter's files and the joined lists (--tag-filter)
    file_util::cleanPrefix(sp_name, sp_base + ".toppic_tag_filter_");
    file_util::delFile(sp_base + ".toppic_topmg_tag_filter");
    file_util::delFile(sp_base + ".toppic_multi_combined");
    file_util::delFile(sp_base + ".toppic_multi_evalue");
    file_util::delFile(sp_base + ".toppic_multi_evalue_all");
    file_util::delFile(sp_base + ".toppic_pre_raw_prsm");
    file_util::delFile(sp_base + ".toppic_multi_ptm_raw_prsm");
    file_util::delFile(sp_base + ".toppic_topmg_evalue_all");
    file_util::delFile(sp_base + ".toppic_topmg_raw_prsm");
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
    file_util::delFile(sp_base + ".toppic_rescore");
    file_util::delFile(sp_base + ".toppic_cluster_prot");
    file_util::delFile(sp_base + ".toppic_first_record_prsm");
    file_util::delFile(sp_base + ".toppic_second_merged");
    file_util::delFile(sp_base + ".toppic_second_cluster");
    file_util::delFile(sp_base + ".toppic_second_prot");
    file_util::delFile(sp_base + ".toppic_second_rescore");
    file_util::delFile(sp_base + ".toppic_second_fdr");
    file_util::delFile(sp_base + ".toppic_prsm_target");
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

namespace {
// The start of the search of one spectrum file, whichever stage searches
// first: the start time, the echo of the arguments, the base data and the
// check that TopFD's feature file exists. Returns false if it does not.
bool TopPIC_init(std::map<std::string, std::string>& arguments) {
  std::time_t start = time(nullptr);
  char buf[50];
  std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));

  arguments["startTime"] = buf;
  // the base data first: the printout takes the masses of C57 and C58 from it
  std::string resource_dir = arguments["resourceDir"];
  base_data::init(resource_dir);
  ToppicArgument::outputArguments(std::cout, " ", arguments);

  EnvBase::initBase(resource_dir);

  LOG_DEBUG("Init base data completed");

  std::string feature_file_name =
      file_util::basename(arguments["spectrumFileName"]) + ".feature";
  if (arguments["useFeatureFile"] == "true") {
    if (!std::filesystem::exists(feature_file_name)) {
      LOG_ERROR("TopFD feature file " << feature_file_name
                                      << " does not exist!");
      LOG_ERROR(
          "Please use the -x option in the command line or check 'No TopFD "
          "feature file' in the GUI.");
      return false;
    }
  }
  return true;
}

// The research switches of --mass-calib, from the environment variable
// TOPPIC_MASS_CALIB_PARA ("name=value,..."): max_evalue (1e-5, the largest
// E-value of a calibrant PrSM), tol_ppm (20, the tolerance of the fragment
// calibrants) and min_tole (10, the smallest precursor and fragment tolerance
// in ppm of the calibration search).
struct MassCalibPara {
  double max_evalue = 1e-5;
  double tol_ppm = 20.0;
  double min_tole = 10.0;
};

MassCalibPara getMassCalibPara() {
  MassCalibPara para;
  const char* env = std::getenv("TOPPIC_MASS_CALIB_PARA");
  if (env == nullptr) {
    return para;
  }
  for (const std::string& item : str_util::split(env, ",")) {
    size_t eq = item.find('=');
    if (eq == std::string::npos) {
      continue;
    }
    std::string name = item.substr(0, eq);
    double value = std::stod(item.substr(eq + 1));
    if (name == "max_evalue") {
      para.max_evalue = value;
    } else if (name == "tol_ppm") {
      para.tol_ppm = value;
    } else if (name == "min_tole") {
      para.min_tole = value;
    } else {
      LOG_ERROR("Unknown TOPPIC_MASS_CALIB_PARA switch " << name << "!");
      exit(EXIT_FAILURE);
    }
  }
  return para;
}

// The rescoring's treatment of the PrSMs with the maximum E-value
// (prsm_rescore::EValueFloor): on unless --disable-rescore-evalue-floor;
// research switches from the environment variable TOPPIC_RESCORE_PARA
// ("name=value,..."): floor_train (0; 1: the model is trained on every PrSM
// and only their scores are replaced) and floor_value (their rescore; by
// default the lowest double, -DBL_MAX).
prsm_rescore::EValueFloor getRescoreEValueFloor(
    std::map<std::string, std::string>& arguments) {
  prsm_rescore::EValueFloor floor;
  floor.enabled = arguments["rescoreEvalueFloor"] != "false";
  const char* env = std::getenv("TOPPIC_RESCORE_PARA");
  if (env == nullptr) {
    return floor;
  }
  for (const std::string& item : str_util::split(env, ",")) {
    size_t eq = item.find('=');
    if (eq == std::string::npos) {
      continue;
    }
    std::string name = item.substr(0, eq);
    double value = std::stod(item.substr(eq + 1));
    if (name == "floor_train") {
      floor.train = value != 0;
    } else if (name == "floor_value") {
      floor.value = value;
    } else {
      LOG_ERROR("Unknown TOPPIC_RESCORE_PARA switch " << name << "!");
      exit(EXIT_FAILURE);
    }
  }
  return floor;
}

// The cutoffs at which a rescored FDR computation's two rounds (E-value and
// rescore) count the targets that pass: each level's cutoff value when its
// cutoff type is FDR, else 1% (the level is then cut by the E-value and its
// FDR only reported).
prsm_fdr::FdrCutoffs getFdrCutoffs(
    std::map<std::string, std::string>& arguments) {
  auto value = [&](const std::string& type, const std::string& val) {
    return arguments[type] == "FDR" ? std::stod(arguments[val]) : 0.01;
  };
  prsm_fdr::FdrCutoffs cutoffs;
  cutoffs.prsm = value("cutoffSpectralType", "cutoffSpectralValue");
  cutoffs.proteoform = value("cutoffProteoformType", "cutoffProteoformValue");
  cutoffs.protein = value("cutoffProteinType", "cutoffProteinValue");
  return cutoffs;
}

// a tolerance as the arguments write it ("10", "7.5"): the index file names
// carry the fragment tolerance's string
std::string toleranceStr(double value) {
  std::ostringstream out;
  out << value;
  return out.str();
}

// By default the generating function of the variable PTM PrSMs takes the
// residues and their forms with the variable PTMs of the file (toppic_test
// experiment 129), unless the hidden --disable-var-ptm-evalue-mod-residues
// is given; toppic's E-value passes call this for their TdgfMng.
void setVarPtmModResidues(std::map<std::string, std::string>& arguments,
                          const TdgfMngPtr& tdgf_mng_ptr) {
  std::string var_ptm_file_name = arguments["variablePtmFileName"];
  if (arguments["varPtmEvalueModResidues"] == "true" &&
      var_ptm_file_name != "" && std::stoi(arguments["variablePtmNum"]) >= 1) {
    tdgf_mng_ptr->var_ptm_mod_ptr_vec_ =
        mod_util::readAnywhereModTxt(var_ptm_file_name);
  }
}
}  // namespace

// protein filtering + database searching + E-value computation
// The tag filter (filter/tag, hidden --tag-filter; toppic_test experiment
// 150): the stages whose candidate lists it joins.
std::set<TagFilterMng::Stage> getTagFilterStages(
    std::map<std::string, std::string>& arguments) {
  std::set<TagFilterMng::Stage> stages;
  TagFilterMng::parseStages(arguments["tagFilterStages"], &stages);
  return stages;
}

// The tag filter of one stage: its settings from the arguments, the run
// and, for the var, one and multi stages, the joining of its candidates to
// the stage's filter file <filter_ext> into <filter_ext>_tag, the name the
// search then reads. The topmg stage adds its file to the graph filter
// merge itself. Returns the file extension the search reads.
std::string runTagFilter(std::map<std::string, std::string>& arguments,
                         const PrsmParaPtr& prsm_para_ptr,
                         TagFilterMng::Stage stage, int thread_num,
                         const std::string& filter_ext,
                         const std::set<int>& search_spec_ids,
                         bool search_all, const std::set<int>& skip_spec_ids) {
  std::string stage_name = TagFilterMng::getStageName(stage);
  std::string tag_ext = stage == TagFilterMng::Stage::kTopmg
                            ? "toppic_topmg_tag_filter"
                            : "toppic_tag_filter_" + stage_name;
  TagFilterMngPtr mng_ptr = std::make_shared<TagFilterMng>(
      prsm_para_ptr, stage, thread_num, tag_ext);
  mng_ptr->tag_len_ = std::stoi(arguments["tagFilterLength"]);
  mng_ptr->cand_num_ = std::stoi(arguments["tagFilterCandNum"]);
  mng_ptr->var_ptm_file_name_ = arguments["variablePtmFileName"];
  mng_ptr->var_ptm_num_ = std::stoi(arguments["variablePtmNum"]);
  mng_ptr->search_all_ = search_all;
  mng_ptr->search_spec_ids_ = search_spec_ids;
  mng_ptr->skip_spec_ids_ = skip_spec_ids;
  std::cout << "Tag filtering (" << stage_name << ") - started." << std::endl;
  tag_filter_processor::process(mng_ptr);
  std::cout << "Tag filtering (" << stage_name << ") - finished." << std::endl;
  if (stage == TagFilterMng::Stage::kTopmg) {
    return tag_ext;
  }
  std::string out_ext = filter_ext + "_tag";
  std::cout << "Tag filtering (" << stage_name
            << ") - joining the candidate lists - started." << std::endl;
  tag_filter_merge::process(prsm_para_ptr->getSpectrumFileName(), filter_ext,
                            tag_ext, out_ext, stage);
  std::cout << "Tag filtering (" << stage_name
            << ") - joining the candidate lists - finished." << std::endl;
  return out_ext;
}

std::string TopPIC_identify(std::map<std::string, std::string>& arguments) {
  try {
    if (!TopPIC_init(arguments)) {
      return "";
    }
    std::set<TagFilterMng::Stage> tag_stages = getTagFilterStages(arguments);

    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    LOG_DEBUG("Number of top PrSMs to retain for each spectrum: " << n_top);
    int var_ptm_num = std::stoi(arguments["variablePtmNum"]);
    int shift_num = std::stoi(arguments["shiftNumber"]);
    std::string var_ptm_file_name = arguments["variablePtmFileName"];

    double max_shift_mass = std::stod(arguments["maxShiftMass"]);
    double min_shift_mass = std::stod(arguments["minShiftMass"]);

    int thread_num = std::stoi(arguments["threadNumber"]);

    // var_ptm_type number is used for E-value computation
    int var_ptm_type_num = 0;
    // and the number of proteoform N-terminal variable PTM masses
    int n_term_var_ptm_type_num = 0;

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

    bool use_var_ptm = var_ptm_num >= 1 && var_ptm_file_name != "";
    bool use_approx_spec = arguments["useApproxSpectra"] == "true";
    // The E-value stage's CountTestNum (a 9 s database scan on one thread,
    // cached for the later E-value stages) is built in a background thread
    // during the filter stages; its parameters do not depend on the
    // variable PTM type number (EValueProcessor::prepareCountTestNum,
    // toppic_test experiment 150).
    std::thread count_test_num_thread([&]() {
      TdgfMngPtr early_tdgf_mng_ptr = std::make_shared<TdgfMng>(
          prsm_para_ptr, shift_num,
          std::max(std::abs(max_shift_mass), std::abs(min_shift_mass)), 0,
          thread_num, "toppic_combined", "toppic_evalue");
      EValueProcessor::prepareCountTestNum(early_tdgf_mng_ptr);
    });
    // The zero-shift and variable PTM filters score the same four index
    // files of every block, so by default they run as one pass over the
    // blocks (zero_var_ptm_filter_processor, toppic_test experiment 150);
    // the hidden --separate-var-ptm-filter runs them one after the other,
    // as -A (approximate spectra) always does. The candidates are the same.
    bool combined_filter = use_var_ptm && !use_approx_spec &&
                           arguments["separateVarPtmFilter"] != "true";
    // the variable PTM filter's manager, built where its stage starts
    VarPtmFilterMngPtr var_filter_mng_ptr;
    auto makeVarFilterMng = [&](int var_filter_thread_num) {
      var_filter_mng_ptr = std::make_shared<VarPtmFilterMng>(
          prsm_para_ptr, index_file_para, var_ptm_file_name, var_ptm_num,
          var_filter_thread_num, use_approx_spec, "toppic_var_filter");
      var_ptm_type_num = var_filter_mng_ptr->getSingleShiftNum();
      n_term_var_ptm_type_num = var_filter_mng_ptr->getNTermShiftNum();
      // the hidden --var-ptm-filter-cand-num <complete>,<prefix/suffix>,
      // <internal> replaces the candidates kept per spectrum (5,20,30);
      // ToppicArgument has checked the three numbers
      if (arguments["varPtmFilterCandNum"] != "") {
        std::vector<std::string> nums =
            str_util::split(arguments["varPtmFilterCandNum"], ",");
        var_filter_mng_ptr->comp_num_ = std::stoi(nums[0]);
        var_filter_mng_ptr->pref_suff_num_ = std::stoi(nums[1]);
        var_filter_mng_ptr->internal_num_ = std::stoi(nums[2]);
      }
    };

    // with one unexpected shift allowed, the one-shift filter joins the
    // combined pass too (its index files are the zero-shift ones cut at its
    // column limit)
    bool combined_one_filter = combined_filter && shift_num >= 1;
    const std::string combined_stage_name =
        combined_one_filter
            ? "Zero unexpected shift, variable PTM and one unexpected shift "
              "filtering"
            : "Zero unexpected shift and variable PTM filtering";
    if (combined_filter) {
      std::cout << combined_stage_name << " - started." << std::endl;
    } else {
      std::cout << "Zero unexpected shift filtering - started." << std::endl;
    }
    // The filtering steps require a large amount of memory, so each one asks
    // mem_check for its own thread number where it starts, from its own
    // per-thread budget and the memory available then; mem_check reports the
    // number when it differs from the one it announced at the start of the
    // run. The thread count cannot change the output of a filtering step -
    // each (spectrum, block) score is a pure function of the index and the
    // spectrum, and the in-memory block combiner is order-independent
    // (toppic_test reports 0022/0066).
    int zero_filter_thread_num = mem_check::getStageThreadNum(
        "Zero unexpected shift filtering", thread_num);
    LOG_DEBUG("Zero-shift filter thread number " << zero_filter_thread_num);
    ZeroPtmFilterMngPtr zero_filter_mng_ptr =
        std::make_shared<ZeroPtmFilterMng>(prsm_para_ptr, index_file_para,
                                           zero_filter_thread_num,
                                           "toppic_zero_filter");
    if (combined_filter) {
      makeVarFilterMng(zero_filter_thread_num);
      OnePtmFilterMngPtr combined_one_filter_mng_ptr;
      if (combined_one_filter) {
        combined_one_filter_mng_ptr = std::make_shared<OnePtmFilterMng>(
            prsm_para_ptr, index_file_para, "toppic_one_filter",
            zero_filter_thread_num, min_shift_mass, max_shift_mass);
      }
      zero_var_ptm_filter_processor::process(zero_filter_mng_ptr,
                                             var_filter_mng_ptr,
                                             combined_one_filter_mng_ptr);
      std::cout << combined_stage_name << " - finished." << std::endl;
    } else {
      zero_ptm_filter_processor::process(zero_filter_mng_ptr);
      std::cout << "Zero unexpected shift filtering - finished." << std::endl;
    }

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

    if (use_var_ptm) {
      if (!combined_filter) {
        std::cout << "Variable PTM filtering - started." << std::endl;
        int var_filter_thread_num =
            mem_check::getStageThreadNum("Variable PTM filtering", thread_num);
        LOG_DEBUG("Variable PTM filter thread number "
                  << var_filter_thread_num);
        makeVarFilterMng(var_filter_thread_num);
        var_ptm_filter_processor::process(var_filter_mng_ptr);
        std::cout << "Variable PTM filtering - finished." << std::endl;
      }

      // the tag filter's candidates join the list (--tag-filter var)
      std::string var_filter_ext = "toppic_var_filter";
      if (tag_stages.count(TagFilterMng::Stage::kVar) > 0) {
        var_filter_ext = runTagFilter(
            arguments, prsm_para_ptr, TagFilterMng::Stage::kVar, thread_num,
            var_filter_ext, std::set<int>(), true, std::set<int>());
      }

      std::cout << "Variable PTM search - started." << std::endl;
      VarPtmSearchMngPtr var_ptm_search_mng_ptr =
          std::make_shared<VarPtmSearchMng>(
              prsm_para_ptr, n_top, var_ptm_file_name, var_ptm_num, thread_num,
              var_filter_ext, "toppic_var_ptm");
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
      if (!combined_one_filter) {
        std::cout << "One unexpected shift filtering - started." << std::endl;
        int one_filter_thread_num = mem_check::getStageThreadNum(
            "One unexpected shift filtering", thread_num);
        LOG_DEBUG("One-shift filter thread number " << one_filter_thread_num);
        OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
            std::make_shared<OnePtmFilterMng>(
                prsm_para_ptr, index_file_para, "toppic_one_filter",
                one_filter_thread_num, min_shift_mass, max_shift_mass);
        one_ptm_filter_processor::process(one_ptm_filter_mng_ptr);
        std::cout << "One unexpected shift filtering - finished."
                  << std::endl;
      }

      // the tag filter's candidates join the list (--tag-filter one)
      std::string one_filter_ext = "toppic_one_filter";
      if (tag_stages.count(TagFilterMng::Stage::kOne) > 0) {
        one_filter_ext = runTagFilter(
            arguments, prsm_para_ptr, TagFilterMng::Stage::kOne, thread_num,
            one_filter_ext, std::set<int>(), true, std::set<int>());
      }

      std::cout << "One unexpected shift search - started." << std::endl;
      int search_shift_num = 1;
      PtmSearchMngPtr one_search_mng_ptr = std::make_shared<PtmSearchMng>(
          prsm_para_ptr, n_top, min_shift_mass, max_shift_mass,
          search_shift_num, thread_num, one_filter_ext,
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
    arguments["nTermVarPtmTypeNum"] = std::to_string(n_term_var_ptm_type_num);

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
    count_test_num_thread.join();
    TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(
        prsm_para_ptr, shift_num,
        std::max(std::abs(max_shift_mass), std::abs(min_shift_mass)),
        var_ptm_type_num, thread_num, "toppic_combined", "toppic_evalue");
    tdgf_mng_ptr->n_term_var_ptm_type_num_ = n_term_var_ptm_type_num;
    setVarPtmModResidues(arguments, tdgf_mng_ptr);
    EValueProcessorPtr processor =
        std::make_shared<EValueProcessor>(tdgf_mng_ptr);
    processor->init();
    // compute E-value for a set of prsms each run (or one PrSM at a time
    // at its own precursor with --evalue-complement-tolerance separate)
    processor->process(prsm_para_ptr->useEvalueSeparate());
    processor = nullptr;
    std::cout << "E-value computation - finished." << std::endl;

    // With the hidden --deamid-selector and deamidation as a variable PTM, a
    // spectrum's zero-PTM PrSM and deamidation PrSM of the same protein and
    // positions (the unmodified proteoform one isotope below the precursor
    // and the deamidated one at it) are compared by matched fragments and
    // precursor mass errors, and only one kind goes on to the top PrSM
    // selection (toppic_test experiment 123). Off by default.
    std::string evalue_ext = "toppic_evalue";
    if (arguments["deamidSelector"] == "true" && var_ptm_num >= 1 &&
        var_ptm_file_name != "") {
      std::set<std::string> deamid_ptm_names =
          prsm_deamid_selector::getDeamidPtmNames(var_ptm_file_name);
      if (!deamid_ptm_names.empty()) {
        std::cout << "Deamidation PrSM selection - started." << std::endl;
        prsm_deamid_selector::process(sp_file_name, evalue_ext,
                                      "toppic_evalue_deamid",
                                      deamid_ptm_names);
        evalue_ext = "toppic_evalue_deamid";
        std::cout << "Deamidation PrSM selection - finished." << std::endl;
      }
    }

    std::cout << "Top PrSM selecting - started" << std::endl;
    prsm_top_selector::process(sp_file_name, evalue_ext, "toppic_pre_raw_prsm",
                               n_top);
    std::cout << "Top PrSM selecting - finished." << std::endl;

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    return "";
  }
  return "toppic_pre_raw_prsm";
}

// TopMG's identification hosted in toppic, behind the hidden --topmg option
// (pipeline or topmg-only). Every stage of TopMG_identify, in its order: the
// database preprocessing and the spectrum index, the ASF-One PTM filtering
// with the variable modification file, the ASF-Diagonal filtering when it is
// asked for, the merge of their results, the proteoform-graph alignment, its
// post-processing, the E-value computation by the generating function (by
// MCMC with the hidden --topmg-use-mcmc) and the top PrSM selection. Only
// the initialisation of the base data and the echo of the arguments are left
// out, because toppic has done them before this runs.
//
// The files are named toppic_topmg_*, so they collide neither with TopMG's
// own nor with toppic's. The graph alignment searches the spectra that the
// set raw_prsm_ext names (the stages before this one) left unidentified, and
// the top PrSM selection chooses between its PrSMs and that set's, writing
// toppic_topmg_raw_prsm, which the post-processing reads. With the hidden
// --topmg topmg-only there is no set before this stage (raw_prsm_ext is
// empty): every spectrum is searched and only this stage's PrSMs are
// selected.
std::string TopPIC_topmg_identify(
    std::map<std::string, std::string>& arguments,
    const std::string& raw_prsm_ext) {
  try {
    std::string db_file_name = arguments["databaseFileName"];
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string var_mod_file_name = arguments["variablePtmFileName"];
    if (var_mod_file_name == "") {
      LOG_ERROR("TopMG's identification needs a variable modification file. "
                "Please use the -b option.");
      return raw_prsm_ext;
    }

    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    // the ASF filtering of these stages keeps this many results per
    // spectrum: TopMG's own number, not toppic's filteringResultNumber
    int filter_result_num = std::stoi(arguments["asfFilteringResultNumber"]);
    double max_shift_mass = std::stod(arguments["maxShiftMass"]);
    int thread_num = std::stoi(arguments["threadNumber"]);

    bool decoy = false;
    if (arguments["searchType"] == "TARGET+DECOY") {
      decoy = true;
    }
    int db_block_size = std::stoi(arguments["databaseBlockSize"]);
    int max_frag_len = std::stoi(arguments["maxFragmentLength"]);
    int min_block_num = std::stoi(arguments["minBlockNum"]);
    // every step of dbPreprocess is guarded by an existence check, so calling
    // it again after TopPIC_identify rebuilds nothing and does not reshuffle
    // the decoy database
    fasta_util::dbPreprocess(ori_db_file_name, db_file_name, decoy,
                             db_block_size, max_frag_len, min_block_num);
    msalign_util::geneSpIndex(sp_file_name);

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    // the unexpected shifts of the graph alignment: TopMG's own number, not
    // toppic's --num-shift
    int shift_num = prsm_para_ptr->getTopmgShiftNum();
    // The +-1 Da precursor error of the one-shift filtering is off in SpPara
    // and in toppic's arguments alike, so these stages filter as TopMG
    // filters without setting it here. It used to be on in SpPara, which is
    // what toppic_test experiment 71 had to work around.

    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);

    // The spectra to search: the ones no earlier stage has identified, as the
    // multiple-shift stage chooses them. A spectrum counts as identified when
    // its best PrSM in the set this stage was given has an E-value strictly
    // below the skip E-value, which with a decoy database is the smallest
    // E-value of a decoy PrSM of that set. Every spectrum a decoy could
    // compete with is searched, so TopMG's decoys are those of a full search
    // over the spectra it searches, and the spectra no decoy can reach keep
    // the PrSM they have. The hidden --multi-skip-evalue overrides it. A
    // spectrum with no PrSM at all is searched, which is why the list is built
    // from every spectrum of the file and not from the PrSMs.
    std::string base_name = file_util::basename(sp_file_name);
    std::set<int> search_spec_ids;
    // the ids of the file, which need not be 0 .. spectrum number - 1
    // (issue #598)
    std::vector<int> spec_ids = msalign_util::readSpIds(sp_file_name);
    if (raw_prsm_ext.empty()) {
      // --topmg topmg-only: no earlier stage has searched, every spectrum
      // is searched here
      search_spec_ids.insert(spec_ids.begin(), spec_ids.end());
      std::cout << "TopMG's identification - " << search_spec_ids.size()
                << " spectra to search, the only search (--topmg topmg-only)."
                << std::endl;
    } else {
      std::set<int> identified_spec_ids;
      PrsmStrPtrVec prev_ptrs =
          prsm_reader_util::readAllPrsmStrs(base_name + "." + raw_prsm_ext);
      double skip_evalue;
      std::string skip_source;
      if (arguments["multiSkipEvalue"] != "") {
        skip_evalue = std::stod(arguments["multiSkipEvalue"]);
        skip_source = "--multi-skip-evalue";
      } else if (arguments["searchType"] == "TARGET+DECOY") {
        skip_evalue = std::numeric_limits<double>::max();
        for (size_t i = 0; i < prev_ptrs.size(); i++) {
          if (prev_ptrs[i]->getSeqName().find("DECOY_") == 0 &&
              prev_ptrs[i]->getEValue() < skip_evalue) {
            skip_evalue = prev_ptrs[i]->getEValue();
          }
        }
        skip_source = "the smallest E-value of a decoy PrSM";
      } else {
        skip_evalue = 0.01;
        skip_source = "the default without a decoy database";
      }
      int skipped_target_num = 0;
      int skipped_decoy_num = 0;
      for (size_t i = 0; i < prev_ptrs.size(); i++) {
        if (prev_ptrs[i]->getEValue() < skip_evalue) {
          int spec_id = prev_ptrs[i]->getSpectrumId();
          if (identified_spec_ids.insert(spec_id).second) {
            if (prev_ptrs[i]->getSeqName().find("DECOY_") == 0) {
              skipped_decoy_num++;
            } else {
              skipped_target_num++;
            }
          }
        }
      }
      for (size_t i = 0; i < spec_ids.size(); i++) {
        if (identified_spec_ids.count(spec_ids[i]) == 0) {
          search_spec_ids.insert(spec_ids[i]);
        }
      }
      std::cout << "TopMG's identification - " << search_spec_ids.size()
                << " spectra to search, " << identified_spec_ids.size()
                << " identified by the earlier stages (" << skipped_target_num
                << " target, " << skipped_decoy_num << " decoy; E-value < "
                << skip_evalue << ", " << skip_source << ")." << std::endl;
    }

    std::vector<std::string> input_exts;

    std::cout << "TopMG ASF-One PTM filtering - started." << std::endl;
    int filter_thread_num =
        mem_check::getStageThreadNum("ASF-One PTM filtering", thread_num);
    LOG_DEBUG("Filter thread number " << filter_thread_num);
    int var_ptm_num_in_filtering = 1;
    OnePtmFilterMngPtr one_ptm_filter_mng_ptr =
        std::make_shared<OnePtmFilterMng>(
            prsm_para_ptr, index_file_para, "toppic_topmg_one_filter",
            filter_thread_num, -max_shift_mass, max_shift_mass,
            var_mod_file_name, var_ptm_num_in_filtering);
    one_ptm_filter_mng_ptr->inte_num_ = 4;
    one_ptm_filter_mng_ptr->pref_suff_num_ = 4;
    one_ptm_filter_mng_ptr->comp_num_ = 4;
    // The filter scores only the spectra the graph alignment searches. The
    // candidates of a spectrum do not depend on the other spectra, and the
    // graph alignment drops those of every other spectrum, so from the graph
    // alignment on nothing changes; the filtering files hold the searched
    // spectra only.
    one_ptm_filter_mng_ptr->search_spec_ids_ = search_spec_ids;
    one_ptm_filter_mng_ptr->search_all_ = false;
    one_ptm_filter_processor::process(one_ptm_filter_mng_ptr);
    std::cout << "TopMG ASF-One PTM filtering - finished." << std::endl;

    input_exts.push_back("toppic_topmg_one_filter_COMPLETE");
    input_exts.push_back("toppic_topmg_one_filter_PREFIX");
    input_exts.push_back("toppic_topmg_one_filter_SUFFIX");
    input_exts.push_back("toppic_topmg_one_filter_INTERNAL");

    if (arguments["useAsfDiag"] == "true") {
      std::cout << "TopMG ASF-Diagonal PTM filtering - started." << std::endl;
      int diag_filter_thread_num = mem_check::getStageThreadNum(
          "ASF-Diagonal PTM filtering", thread_num);
      LOG_DEBUG("Diagonal filter thread number " << diag_filter_thread_num);
      DiagFilterMngPtr diag_filter_mng_ptr = std::make_shared<DiagFilterMng>(
          prsm_para_ptr, index_file_para, filter_result_num,
          diag_filter_thread_num, "toppic_topmg_multi_filter",
          var_mod_file_name, 1);
      DiagFilterProcessorPtr diag_filter_processor =
          std::make_shared<DiagFilterProcessor>(diag_filter_mng_ptr);
      diag_filter_processor->process();
      diag_filter_processor = nullptr;
      std::cout << "TopMG ASF-Diagonal filtering - finished." << std::endl;

      input_exts.push_back("toppic_topmg_multi_filter");
    }

    // the tag filter's candidates join the graph filter merge
    // (--tag-filter topmg)
    if (getTagFilterStages(arguments).count(TagFilterMng::Stage::kTopmg) >
        0) {
      input_exts.push_back(runTagFilter(
          arguments, prsm_para_ptr, TagFilterMng::Stage::kTopmg, thread_num,
          "", search_spec_ids, false, std::set<int>()));
    }

    std::cout << "TopMG combining filtering results - started." << std::endl;
    SimplePrsmStrMergePtr asf_filter_merger =
        std::make_shared<SimplePrsmStrMerge>(sp_file_name, input_exts,
                                             "toppic_topmg_graph_filter",
                                             20 * input_exts.size());
    asf_filter_merger->process();
    asf_filter_merger = nullptr;
    std::cout << "TopMG combining filtering results - finished." << std::endl;

    int max_mod_num = std::stoi(arguments["variablePtmNum"]);
    int gap = std::stoi(arguments["proteoGraphGap"]);
    int var_ptm_in_gap =
        std::min(std::stoi(arguments["varPtmNumInGap"]), max_mod_num);
    // whole proteins only is the proteoform type COMPLETE alone, which is one
    // of the two values the argument check allows with --topmg
    bool whole_protein_only = (arguments["allowProtType"] == "COMPLETE");
    GraphAlignMngPtr ga_mng_ptr = std::make_shared<GraphAlignMng>(
        prsm_para_ptr, var_mod_file_name, shift_num, max_mod_num, gap,
        var_ptm_in_gap, max_shift_mass,
        std::stod(arguments["topmgAlignErrorTolerance"]), thread_num,
        whole_protein_only, "toppic_topmg_graph_filter", "toppic_topmg_graph_align");
    // the graph alignment searches these spectra and no others
    ga_mng_ptr->search_spec_ids_ = search_spec_ids;
    ga_mng_ptr->search_all_ = false;
    std::cout << "Graph alignment - started." << std::endl;
    GraphAlignProcessorPtr ga_processor_ptr =
        std::make_shared<GraphAlignProcessor>(ga_mng_ptr);
    ga_processor_ptr->process();
    ga_processor_ptr = nullptr;
    std::cout << "Graph alignment - finished." << std::endl;

    std::cout << "Graph alignment post-processing - started." << std::endl;
    GraphPostProcessorPtr ga_post_processor_ptr =
        std::make_shared<GraphPostProcessor>(ga_mng_ptr,
                                             "toppic_topmg_graph_align",
                                             "toppic_topmg_graph_post");
    ga_post_processor_ptr->process();
    ga_post_processor_ptr = nullptr;
    std::cout << "Graph alignment post-processing - finished." << std::endl;

    if (arguments["topmgUseMcmc"] != "true") {
      // The generating function of toppic's own E-value stages. Its table
      // has one layer per unexpected shift, so it takes the graph
      // alignment's shift number; the maximum shift mass the graph alignment
      // allowed bounds the candidate count of a shifted PrSM, and the
      // variable PTM type number that TopPIC_identify counted from the same
      // modification file multiplies it for every variable PTM, as it does
      // for toppic's own PrSMs.
      std::cout << "E-value computation using TDGF - started." << std::endl;
      int var_ptm_type_num = arguments["varPtmTypeNum"] == ""
                                 ? 0
                                 : std::stoi(arguments["varPtmTypeNum"]);
      TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(
          prsm_para_ptr, shift_num, max_shift_mass, var_ptm_type_num,
          thread_num, "toppic_topmg_graph_post", "toppic_topmg_evalue");
      // the truncation N-terminal modifications' own type number, as
      // toppic's own stages count it
      tdgf_mng_ptr->n_term_var_ptm_type_num_ =
          arguments["nTermVarPtmTypeNum"] == ""
              ? 0
              : std::stoi(arguments["nTermVarPtmTypeNum"]);
      EValueProcessorPtr processor =
          std::make_shared<EValueProcessor>(tdgf_mng_ptr);
      processor->init();
      processor->process(prsm_para_ptr->useEvalueSeparate());
      processor = nullptr;
      std::cout << "E-value computation using TDGF - finished." << std::endl;
    } else {
      std::cout << "E-value computation using MCMC - started." << std::endl;
      MCMCMngPtr mcmc_mng_ptr = std::make_shared<MCMCMng>(
          prsm_para_ptr, "toppic_topmg_graph_post", "toppic_topmg_evalue",
          var_mod_file_name, max_mod_num, thread_num);
      DprProcessorPtr processor = std::make_shared<DprProcessor>(mcmc_mng_ptr);
      processor->process();
      processor = nullptr;
      std::cout << "E-value computation using MCMC - finished." << std::endl;
    }

    // the PrSMs of this stage, whose E-values the generating function (or
    // the MCMC) gave, and those of the stage before it, so that the best
    // PrSM of every spectrum is chosen between the two searches
    std::cout << "Merging TopMG and TopPIC PrSMs - started." << std::endl;
    {
      // with --topmg topmg-only there is no stage before this one
      std::vector<std::string> all_exts = {"toppic_topmg_evalue"};
      if (!raw_prsm_ext.empty()) {
        all_exts.push_back(raw_prsm_ext);
      }
      PrsmStrMergePtr all_merge_ptr = std::make_shared<PrsmStrMerge>(
          sp_file_name, all_exts, "toppic_topmg_evalue_all",
          std::numeric_limits<int>::max(), false, false);
      all_merge_ptr->process();
      all_merge_ptr = nullptr;
    }
    std::cout << "Merging TopMG and TopPIC PrSMs - finished." << std::endl;

    std::cout << "TopMG top PrSM selecting - started" << std::endl;
    prsm_top_selector::process(sp_file_name, "toppic_topmg_evalue_all",
                               "toppic_topmg_raw_prsm", n_top);
    std::cout << "TopMG top PrSM selecting - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  // its PrSMs, merged with those of the stage before it, are the set now
  return "toppic_topmg_raw_prsm";
}

// The multiple unexpected shifts stage (-s 2): the diagonal filter and the
// PTM search on the spectra that TopPIC_identify left without a PrSM passing
// the spectrum-level cutoff, then the merge, E-value computation and top
// PrSM selection of the new PrSMs alone.
std::string TopPIC_multi_ptm_search(
    std::map<std::string, std::string>& arguments,
    const std::string& raw_prsm_ext) {
  try {
    // the search and the E-value computation of this stage take the number of
    // unexpected shifts; TopPICProgress calls the stage only when it is 2
    int shift_num = std::stoi(arguments["shiftNumber"]);
    std::string resource_dir = arguments["resourceDir"];
    base_data::init(resource_dir);

    std::string sp_file_name = arguments["spectrumFileName"];
    int n_top = std::stoi(arguments["numOfTopPrsms"]);
    double max_shift_mass = std::stod(arguments["maxShiftMass"]);
    double min_shift_mass = std::stod(arguments["minShiftMass"]);
    int filter_result_num = std::stoi(arguments["filteringResultNumber"]);
    int thread_num = std::stoi(arguments["threadNumber"]);
    int var_ptm_type_num = std::stoi(arguments["varPtmTypeNum"]);

    IndexFileNamePtr file_name_ptr = std::make_shared<IndexFileName>();
    std::string index_file_para = file_name_ptr->geneFileName(arguments);
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);

    // The spectra the earlier stages identified are not searched again: a
    // spectrum is skipped when its best first-pass PrSM (toppic_pre_raw_prsm,
    // targets and decoys alike) has an E-value strictly below the skip
    // E-value. With a decoy database that is the smallest E-value of a
    // first-pass decoy PrSM: every spectrum a decoy could compete with is
    // searched again, so the decoys of the final set are those of a full
    // search, while the spectra no decoy can reach keep their first-pass
    // PrSM (report 0024); without one it is 0.01. The hidden
    // --multi-skip-evalue <x> overrides both.
    std::string base_name = file_util::basename(sp_file_name);
    std::set<int> identified_spec_ids;
    PrsmStrPtrVec first_ptrs =
        prsm_reader_util::readAllPrsmStrs(base_name + ".toppic_pre_raw_prsm");
    double skip_evalue;
    if (arguments["multiSkipEvalue"] != "") {
      skip_evalue = std::stod(arguments["multiSkipEvalue"]);
    } else if (arguments["searchType"] == "TARGET+DECOY") {
      // the smallest E-value of a decoy PrSM
      skip_evalue = std::numeric_limits<double>::max();
      for (size_t i = 0; i < first_ptrs.size(); i++) {
        if (first_ptrs[i]->getSeqName().find("DECOY_") == 0 &&
            first_ptrs[i]->getEValue() < skip_evalue) {
          skip_evalue = first_ptrs[i]->getEValue();
        }
      }
    } else {
      // the default without a decoy database
      skip_evalue = 0.01;
    }
    for (size_t i = 0; i < first_ptrs.size(); i++) {
      if (first_ptrs[i]->getEValue() < skip_evalue) {
        identified_spec_ids.insert(first_ptrs[i]->getSpectrumId());
      }
    }
    std::cout << "Multiple unexpected shifts stage - "
              << identified_spec_ids.size()
              << " spectra identified by the earlier stages are skipped."
              << std::endl;

    std::vector<std::string> input_exts;
    {
      std::cout << "Multiple unexpected shifts filtering - started."
                << std::endl;
      // thread number is used because diagonal filter uses only one index
      int diag_filter_thread_num = mem_check::getStageThreadNum(
          "Multiple unexpected shifts filtering", thread_num);
      LOG_DEBUG("Diagonal filter thread number " << diag_filter_thread_num);
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

      // the tag filter's candidates join the list (--tag-filter multi), on
      // the spectra this stage searches
      std::string multi_filter_ext = "toppic_multi_filter";
      if (getTagFilterStages(arguments).count(TagFilterMng::Stage::kMulti) >
          0) {
        multi_filter_ext = runTagFilter(
            arguments, prsm_para_ptr, TagFilterMng::Stage::kMulti, thread_num,
            multi_filter_ext, std::set<int>(), true, identified_spec_ids);
      }

      std::cout << "Multiple unexpected shifts search - started." << std::endl;
      PtmSearchMngPtr multi_search_mng_ptr = std::make_shared<PtmSearchMng>(
          prsm_para_ptr, n_top, min_shift_mass, max_shift_mass, shift_num,
          thread_num, multi_filter_ext, "toppic_multi_shift");
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
    setVarPtmModResidues(arguments, tdgf_mng_ptr);
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
      // the PrSMs of this stage and those of the stage before it
      std::vector<std::string> all_exts = {"toppic_multi_evalue",
                                           raw_prsm_ext};
      PrsmStrMergePtr all_merge_ptr = std::make_shared<PrsmStrMerge>(
          sp_file_name, all_exts, "toppic_multi_evalue_all",
          std::numeric_limits<int>::max(), false, false);
      all_merge_ptr->process();
      all_merge_ptr = nullptr;
    }
    prsm_top_selector::process(sp_file_name, "toppic_multi_evalue_all",
                               "toppic_multi_ptm_raw_prsm", n_top);
    std::cout << "Multiple-shift top PrSM selecting - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    // the set is still the one this stage was given
    return raw_prsm_ext;
  }
  return "toppic_multi_ptm_raw_prsm";
}

// Post-processing of the identifications. With -s 2 the set is the best PrSM
// of every spectrum: the first-pass PrSM of the skipped spectra and, for the
// searched spectra, the best of their first-pass and multiple-shift PrSMs
// (toppic_multi_ptm_raw_prsm); without -s 2 it is toppic_pre_raw_prsm.
// Either becomes toppic_raw_prsm, the name the steps below read. The set is
// clustered and protein-corrected, recounted against the full spectra, then
// post-matched (the added masses are in every later step), characterized
// and given one FDR; the best PrSM of every spectrum is kept as
// toppic_record_prsm, from which TopPIC_report writes the results.
int TopPIC_postprocess(std::map<std::string, std::string>& arguments,
                       const std::string& raw_prsm_ext) {
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
    // the E-value model of the cluster protein correction below. The PrSMs it
    // scores carry the proteoform of their cluster's best PrSM: one of
    // toppic's own, or with --topmg one of the graph alignment's,
    // which carries up to TopMG's own number of unexpected shifts. The
    // generating function needs a table layer for every shift, so it takes
    // the larger of the two numbers.
    int shift_num = std::stoi(arguments["shiftNumber"]);
    if (arguments["topmgMode"] != "disabled" &&
        arguments["topmgShiftNumber"] != "") {
      shift_num = std::max(shift_num, std::stoi(arguments["topmgShiftNumber"]));
    }
    // The other proteins' PrSMs of a cluster whose best PrSM comes from the
    // graph alignment are dropped, not corrected, only when the graph
    // alignment's E-values come from MCMC (--topmg-use-mcmc): a corrected
    // PrSM gets a generating-function E-value, which is not comparable with
    // the MCMC one of the PrSM it would be counted with.
    bool remove_topmg_others = arguments["topmgMode"] != "disabled" &&
                               arguments["topmgUseMcmc"] == "true";

    // raw_prsm_ext is the file the last stage that produced a PrSM set
    // reported: toppic_pre_raw_prsm, or toppic_multi_ptm_raw_prsm after the
    // multiple-shift stage. It is copied to toppic_raw_prsm, the name every
    // step below reads, unconditionally, so that a toppic_raw_prsm left in
    // the directory by an earlier run is never picked up.
    file_util::copyFile(base_name + "." + raw_prsm_ext,
                        base_name + ".toppic_raw_prsm", true);

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
          sp_file_name, "toppic_raw_prsm", "toppic_proteoform_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    } else {
      prsm_simple_cluster::process(
          db_file_name, sp_file_name, "toppic_raw_prsm",
          prsm_para_ptr->getFixModPtrVec(), "toppic_proteoform_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    }
    std::cout << "Finding PrSM proteoform clusters - finished." << std::endl;
    std::cout << "Finding PrSM protein clusters - started." << std::endl;
    prsm_prot_cluster::process(sp_file_name, "toppic_proteoform_cluster",
                               "toppic_cluster");
    std::cout << "Finding PrSM protein clusters - finished." << std::endl;
    std::string cur_suffix = "toppic_cluster";

    // In every proteoform cluster, the PrSMs whose protein differs from the
    // protein of the cluster's best-E-value PrSM (the same proteoform mass
    // matched to another protein) get that PrSM's proteoform and its tool;
    // their matched fragment numbers and E-values are recomputed on the
    // spectra as the search read them (EnvCNN cutoff, deduplicated masses).
    // The recount,
    // the post mass matching, the PTM characterization and the FDR run on
    // the corrected PrSMs.
    std::cout << "Cluster protein correction - started." << std::endl;
    {
      int var_ptm_type_num = arguments["varPtmTypeNum"] == ""
                                 ? 0
                                 : std::stoi(arguments["varPtmTypeNum"]);
      TdgfMngPtr tdgf_mng_ptr = std::make_shared<TdgfMng>(
          prsm_para_ptr, shift_num,
          std::max(std::abs(max_shift_mass), std::abs(min_shift_mass)),
          var_ptm_type_num, std::stoi(arguments["threadNumber"]), "", "");
      tdgf_mng_ptr->n_term_var_ptm_type_num_ =
          arguments["nTermVarPtmTypeNum"] == ""
              ? 0
              : std::stoi(arguments["nTermVarPtmTypeNum"]);
      setVarPtmModResidues(arguments, tdgf_mng_ptr);
      int removed_num = 0;
      int corrected_num = cluster_prot_corrector::process(
          prsm_para_ptr, tdgf_mng_ptr, cur_suffix, "toppic_cluster_prot",
          remove_topmg_others, &removed_num);
      if (removed_num > 0) {
        std::cout << "Cluster protein correction - removed " << removed_num
                  << " PrSMs of other proteins in clusters whose best PrSM "
                     "comes from the graph alignment." << std::endl;
      }
      std::cout << "Cluster protein correction - corrected " << corrected_num
                << " PrSMs." << std::endl;
    }
    std::cout << "Cluster protein correction - finished." << std::endl;
    cur_suffix = "toppic_cluster_prot";

    std::cout << "Recounting matched masses and fragments - started."
              << std::endl;
    // the search and the cluster protein correction counted only the masses
    // above the EnvCNN cutoff; recount every PrSM against the full spectra
    prsm_match_num_recount::process(prsm_para_ptr, cur_suffix,
                                    "toppic_recount");
    std::cout << "Recounting matched masses and fragments - finished."
              << std::endl;
    cur_suffix = "toppic_recount";

    if (arguments["postMassMatch"] == "true") {
      // On the best PrSM of every spectrum, targets and decoys alike, after
      // the clustering, the cluster protein correction and the recount: it
      // adds matched masses and fragments but changes no E-value or cluster,
      // so the FDR and the cutoffs give the identifications of a run without
      // it.
      std::cout << "Post mass matching - started." << std::endl;
      int min_peak_num = std::stoi(arguments["postMinPeakNum"]);
      // the EnvCNN model scores the added masses
      onnx_env_cnn::initModel(resource_dir,
                              std::stoi(arguments["threadNumber"]));
      // the peaks of TopFD's database, read only, with the m/z values
      // corrected in memory after --mass-calib (TopPIC_mass_calib); the
      // envelopes go to toppic's database of the file, <base>_toppic.sqlite
      // (TopPICProgress_multi_file), none with --no-sql
      sp_file_name = prsm_post_mass_match::process(
          prsm_para_ptr, cur_suffix, cur_suffix, min_peak_num,
          arguments["topfdDbFileName"], arguments["massCalibModelFileName"],
          arguments["noSql"] == "true" ? "" : arguments["toppicDbFileName"]);
      // The steps below read the spectra with the added masses and the PrSMs
      // with the recounted matches, and name their files after the new
      // spectrum file.
      arguments["spectrumFileName"] = sp_file_name;
      prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
      msalign_util::geneSpIndex(sp_file_name);
      std::cout << "Post mass matching - finished." << std::endl;
    }

    if (arguments["localPtmFileName"] != "") {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng = std::make_shared<LocalMng>(
          prsm_para_ptr, std::stod(arguments["localThreshold"]),
          arguments["localPtmFileName"], min_shift_mass, max_shift_mass,
          cur_suffix, "toppic_cluster_local");
      if (arguments["localAlign"] == "true") {
        // the default: the variable PTM search of stat/local_align; hidden
        // --local-original: the MIScore candidate search of stat/local
        LocalAlignProcessorPtr local_ptr =
            std::make_shared<LocalAlignProcessor>(local_mng);
        local_ptr->process();
      } else {
        LocalProcessorPtr local_ptr =
            std::make_shared<LocalProcessor>(local_mng);
        local_ptr->process();
      }
      std::cout << "PTM characterization - finished." << std::endl;
      cur_suffix = "toppic_cluster_local";
    }

    bool use_rescore = arguments["rescore"] == "true";
    if (arguments["searchType"] == "TARGET+DECOY" && use_rescore) {
      // the rescoring (--rescore, with -d; off by default): the
      // best PrSM of every spectrum, targets and decoys, after the
      // recount, the post mass matching and the PTM characterization, gets
      // a score of a linear SVM trained with cross-validation, which the
      // FDR ranks by instead of the E-value
      std::cout << "Rescoring - started." << std::endl;
      prsm_rescore::process(prsm_para_ptr, cur_suffix, "toppic_rescore", "",
                            getRescoreEValueFloor(arguments));
      std::cout << "Rescoring - finished." << std::endl;
      cur_suffix = "toppic_rescore";
    }

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      // the FDR step keeps the decoys so that toppic_record_prsm carries
      // them; without -K TopPIC_report removes them
      if (arguments["useFdrGroups"] == "true") {
        // hidden option --fdr-groups: separate FDRs for the PrSM groups
        // (0 shifts without/with variable PTMs, 1 shift, 2 shifts)
        prsm_fdr_groups::process(sp_file_name, cur_suffix,
                                 "toppic_cluster_fdr");
      } else {
        prsm_fdr::process(sp_file_name, cur_suffix, "toppic_cluster_fdr",
                          use_rescore, getFdrCutoffs(arguments),
                          arguments["rescoreFallback"] != "false");
      }
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "toppic_cluster_fdr";
    }
    // Copy the PrSM file the spectrum-level cutoff reads -- the best PrSM of
    // every spectrum after the clustering, the cluster protein correction,
    // the recount, the post mass matching, the PTM characterization and the
    // FDR computation,
    // <base>.<cur_suffix> (toppic_cluster_fdr with -d; toppic_recount, or
    // toppic_cluster_local with -l, without), <base> being the spectrum
    // file's name without its extension (the _post_ms2 name when the post
    // mass matching ran) -- to <base>.toppic_record_prsm, which is kept
    // after the run (cleanToppicDir does not delete it).
    std::string record_base_name = file_util::basename(sp_file_name);
    file_util::copyFile(record_base_name + "." + cur_suffix,
                        record_base_name + ".toppic_record_prsm", true);
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

// Post-processing of the merged PrSMs of a combined run (-c): the merged
// record files of the fractions (<name>_ms2.toppic_raw_prsm) are clustered
// again into proteoforms and proteins across the fractions (the fractions'
// cluster ids overlap) and given one FDR, and the result is written to
// <name>_ms2.toppic_record_prsm for TopPIC_report. No recount, post mass
// matching or PTM characterization: the fractions' values are kept.
int TopPIC_multi_postprocess(std::map<std::string, std::string>& arguments) {
  try {
    std::string sp_file_name = arguments["spectrumFileName"];
    std::string base_name = file_util::basename(sp_file_name);
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name =
        ori_db_file_name + "_idx" + file_util::getFileSeparator() +
        file_util::filenameFromEntirePath(arguments["databaseFileName"]);
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);

    bool is_proteoform_ppm_error = (arguments["proteoformPpmError"] == "true");
    double proteoform_error_tole =
        std::stod(arguments["proteoformErrorTolerance"]);
    if (is_proteoform_ppm_error) {
      proteoform_error_tole = proteoform_error_tole / 1000000;
    }
    std::cout << "Finding PrSM proteoform clusters - started." << std::endl;
    // by protein, precursor mass and proteoform sequence overlap, whether or
    // not the fractions had feature files
    prsm_simple_cluster::process(
        db_file_name, sp_file_name, "toppic_raw_prsm",
        prsm_para_ptr->getFixModPtrVec(), "toppic_proteoform_cluster",
        is_proteoform_ppm_error, proteoform_error_tole);
    std::cout << "Finding PrSM proteoform clusters - finished." << std::endl;
    std::cout << "Finding PrSM protein clusters - started." << std::endl;
    prsm_prot_cluster::process(sp_file_name, "toppic_proteoform_cluster",
                               "toppic_cluster");
    std::cout << "Finding PrSM protein clusters - finished." << std::endl;
    std::string cur_suffix = "toppic_cluster";

    // the rescoring: the merged PrSMs are rescored together on the merged
    // spectra, as in TopPIC_postprocess
    bool use_rescore = arguments["rescore"] == "true";
    if (arguments["searchType"] == "TARGET+DECOY" && use_rescore) {
      std::cout << "Rescoring - started." << std::endl;
      prsm_rescore::process(prsm_para_ptr, cur_suffix, "toppic_rescore", "",
                            getRescoreEValueFloor(arguments));
      std::cout << "Rescoring - finished." << std::endl;
      cur_suffix = "toppic_rescore";
    }

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      // decoys kept, as in TopPIC_postprocess; TopPIC_report removes them
      if (arguments["useFdrGroups"] == "true") {
        prsm_fdr_groups::process(sp_file_name, cur_suffix,
                                 "toppic_cluster_fdr");
      } else {
        prsm_fdr::process(sp_file_name, cur_suffix, "toppic_cluster_fdr",
                          use_rescore, getFdrCutoffs(arguments),
                          arguments["rescoreFallback"] != "false");
      }
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "toppic_cluster_fdr";
    }
    file_util::copyFile(base_name + "." + cur_suffix,
                        base_name + ".toppic_record_prsm", true);
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

// Reporting: the cutoffs, the tables and the SQLite file, from the
// toppic_record_prsm file that TopPIC_postprocess (or, for a combined run,
// TopPIC_multi_postprocess) left. The decoys are removed first unless -K
// keeps them.
int TopPIC_report(std::map<std::string, std::string>& arguments) {
  try {
    std::string base_name = file_util::basename(arguments["spectrumFileName"]);
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name =
        ori_db_file_name + "_idx" + file_util::getFileSeparator() +
        file_util::filenameFromEntirePath(arguments["databaseFileName"]);
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    std::string cur_suffix = "toppic_record_prsm";
    // the variable modifications (-b, and the localization's -B), whose
    // residues decide the variable PTM sites of the proteoform
    // classification
    ModPtrVec var_mod_ptrs;
    for (const char* key : {"variablePtmFileName", "localPtmFileName"}) {
      if (arguments[key] != "") {
        ModPtrVec mod_ptrs = mod_util::readAnywhereModTxt(arguments[key]);
        var_mod_ptrs.insert(var_mod_ptrs.end(), mod_ptrs.begin(),
                            mod_ptrs.end());
      }
    }

    if (arguments["searchType"] == "TARGET+DECOY" &&
        arguments["keepDecoyResults"] != "true") {
      std::cout << "Removing decoy identifications - started." << std::endl;
      prsm_decoy_filter::process(base_name, cur_suffix, "toppic_prsm_target");
      cur_suffix = "toppic_prsm_target";
      std::cout << "Removing decoy identifications - finished." << std::endl;
    }

    std::string spectra_cutoff_type = arguments["cutoffSpectralType"];
    double spectra_cutoff_value;
    std::istringstream(arguments["cutoffSpectralValue"]) >> spectra_cutoff_value;
    std::cout << "PrSM filtering by " << spectra_cutoff_type << " - started."
              << std::endl;
    prsm_cutoff_selector::process(db_file_name, base_name, cur_suffix,
                                  "toppic_prsm_cutoff", spectra_cutoff_type,
                                  spectra_cutoff_value, var_mod_ptrs);
    std::cout << "PrSM filtering by " << spectra_cutoff_type << " - finished."
              << std::endl;

    {
      double cutoff_value;
      std::string form_cutoff_type =
          (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR" : "EVALUE";
      std::cout << "Proteoform filtering by " << form_cutoff_type
                << " - started." << std::endl;
      std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
      prsm_cutoff_selector::process(db_file_name, base_name,
                                    "toppic_prsm_cutoff", "toppic_form_cutoff",
                                    form_cutoff_type, cutoff_value,
                                    var_mod_ptrs);
      std::cout << "Proteoform filtering by " << form_cutoff_type
                << " - finished." << std::endl;
      std::cout << "Selecting top PrSMs for proteoforms - started."
                << std::endl;
      prsm_form_filter::process(db_file_name, base_name, "toppic_form_cutoff",
                                "toppic_form_cutoff_form");
      std::cout << "Selecting top PrSMs for proteoforms - finished."
                << std::endl;

      std::string prot_cutoff_type =
          (arguments["cutoffProteinType"] == "FDR") ? "PROTFDR" : "EVALUE";
      std::cout << "Protein filtering by " << prot_cutoff_type << " - started."
                << std::endl;
      std::istringstream(arguments["cutoffProteinValue"]) >> cutoff_value;
      prsm_cutoff_selector::process(db_file_name, base_name,
                                    "toppic_form_cutoff", "toppic_prot_cutoff",
                                    prot_cutoff_type, cutoff_value,
                                    var_mod_ptrs);
      std::cout << "Protein filtering by " << prot_cutoff_type << " - finished."
                << std::endl;
      std::cout << "Selecting top PrSMs for proteins - started." << std::endl;
      prsm_prot_filter::process(db_file_name, base_name, "toppic_prot_cutoff",
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

    // toppic's database of the file, <base>_toppic.sqlite; none for the
    // combined result of -c, nor with --no-sql
    if (arguments["noSql"] == "true") {
      std::cout << "The PrSM tables are not written to an SQLite database "
                   "(--no-sql)."
                << std::endl;
    } else {
      prsm_sql_output::write(prsm_para_ptr, arguments["toppicDbFileName"],
                             base_name, ori_db_file_name, "toppic_prsm_cutoff",
                             "toppic_form_cutoff_form",
                             "toppic_prot_cutoff_prot");
    }

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

// The identification stages of one spectrum file: toppic's own
// (TopPIC_identify, the multiple-shift stage) or the hosted TopMG's, as
// --topmg says. Returns the extension of the raw PrSM file they leave,
// empty on an error.
std::string TopPIC_identifyAll(std::map<std::string, std::string>& arguments) {
  std::string raw_prsm_ext;
  if (arguments["topmgMode"] == "topmg-only") {
    // --topmg topmg-only: TopMG's identification is the only database
    // search, and toppic's own identification stages (TopPIC_identify and
    // TopPIC_multi_ptm_search) are skipped
    if (!TopPIC_init(arguments)) {
      return "";
    }
    // The variable PTM type number by which the generating function (the
    // hosted stage's E-values and the cluster protein correction) multiplies
    // the candidate number of a PrSM for each variable PTM. TopPIC_identify
    // counts it in its variable-PTM stage, from the same file.
    int var_ptm_type_num = 0;
    // and the number of proteoform N-terminal PTM masses (Position N-term),
    // by which it multiplies the candidate number of a PrSM with a
    // truncation N-terminal modification (VarPtmFilterMng::getNTermShiftNum)
    int n_term_var_ptm_type_num = 0;
    if (std::stoi(arguments["variablePtmNum"]) >= 1 &&
        arguments["variablePtmFileName"] != "") {
      var_ptm_type_num = static_cast<int>(
          mod_util::readModTxtToShiftList(arguments["variablePtmFileName"])
              .size());
      ModPtrVec n_term_mod_ptrs =
          mod_util::readNTermModTxt(arguments["variablePtmFileName"]);
      std::vector<double> n_term_shifts;
      for (size_t i = 0; i < n_term_mod_ptrs.size(); i++) {
        double shift = n_term_mod_ptrs[i]->getShift();
        if (std::find(n_term_shifts.begin(), n_term_shifts.end(), shift) ==
            n_term_shifts.end()) {
          n_term_shifts.push_back(shift);
        }
      }
      n_term_var_ptm_type_num = static_cast<int>(n_term_shifts.size());
    }
    arguments["varPtmTypeNum"] = std::to_string(var_ptm_type_num);
    arguments["nTermVarPtmTypeNum"] = std::to_string(n_term_var_ptm_type_num);
    raw_prsm_ext = TopPIC_topmg_identify(arguments, "");
  } else {
    raw_prsm_ext = TopPIC_identify(arguments);
    if (raw_prsm_ext.empty()) {
      return "";
    }
    // the multiple unexpected shifts stage is the -s 2 stage
    if (std::stoi(arguments["shiftNumber"]) >= 2) {
      raw_prsm_ext = TopPIC_multi_ptm_search(arguments, raw_prsm_ext);
    }
    if (arguments["topmgMode"] == "pipeline") {
      raw_prsm_ext = TopPIC_topmg_identify(arguments, raw_prsm_ext);
    }
  }
  return raw_prsm_ext;
}

// The second precursor search (hidden --second-precursor; prsm_second_prec).
// After the first pass's post-processing (arguments["spectrumFileName"] is
// its post-matched spectrum file; ori_sp_file_name the file it searched):
// the residual spectra <base>_resid_ms2.msalign are written (the second
// precursor first, the masses of the accepted first-pass PrSM removed) and
// searched and post-processed like the first pass (without the rescoring,
// which runs on the merged set); then the two passes' record PrSMs are
// merged into one set of (spectrum id, precursor id) units, clustered again
// (no cluster protein correction, recount or post mass matching: each pass
// did its own), rescored (with --rescore; the second units with
// their residual spectra), given one FDR and written as the file's
// toppic_record_prsm (the first pass's is kept as toppic_first_record_prsm),
// which TopPIC_report reads; the masses the second pass's post mass
// matching added join the first pass's post-matched spectrum file, one
// spectrum per scan. The residual pass's files are cleaned up as the first
// pass's.
int TopPIC_second_precursor(std::map<std::string, std::string>& arguments,
                            const std::string& ori_sp_file_name) {
  try {
    std::string post_sp_name = arguments["spectrumFileName"];
    std::string post_base = file_util::basename(post_sp_name);
    std::string ori_base = file_util::basename(ori_sp_file_name);
    if (str_util::endsWith(ori_base, "_ms2")) {
      ori_base = ori_base.substr(0, ori_base.size() - 4);
    }
    std::string resid_sp_name = ori_base + "_resid_ms2.msalign";
    std::string resid_base = file_util::basename(resid_sp_name);
    bool keep_temp_files = (arguments["keepTempFiles"] == "true");
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];

    std::cout << "Second precursor search - started." << std::endl;
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    int spec_num = prsm_second_prec::writeResidualSpectra(
        prsm_para_ptr, post_base + ".toppic_record_prsm", ori_sp_file_name,
        arguments["cutoffSpectralType"],
        std::stod(arguments["cutoffSpectralValue"]), resid_sp_name);
    if (arguments["useFeatureFile"] == "true") {
      file_util::copyFile(file_util::basename(ori_sp_file_name) + ".feature",
                          resid_base + ".feature", true);
    }
    if (spec_num == 0) {
      std::cout << "Second precursor search - finished (no spectrum with a "
                   "second precursor)."
                << std::endl;
      return 0;
    }

    // the second pass, on its own copy of the arguments
    std::map<std::string, std::string> second_arguments = arguments;
    second_arguments["spectrumFileName"] = resid_sp_name;
    second_arguments["rescore"] = "false";
    std::string raw_prsm_ext = TopPIC_identifyAll(second_arguments);
    if (raw_prsm_ext.empty()) {
      return 1;
    }
    TopPIC_postprocess(second_arguments, raw_prsm_ext);
    std::string resid_post_sp_name = second_arguments["spectrumFileName"];
    std::string resid_post_base = file_util::basename(resid_post_sp_name);

    // the merge of the passes' PrSMs into (spectrum, precursor) units
    std::cout << "Second precursor search - merging the passes." << std::endl;
    file_util::copyFile(post_base + ".toppic_record_prsm",
                        post_base + ".toppic_first_record_prsm", true);
    prsm_second_prec::mergePrsms(
        post_base + ".toppic_first_record_prsm",
        resid_post_base + ".toppic_record_prsm",
        post_base + ".toppic_second_merged",
        file_util::filenameFromEntirePath(post_sp_name));

    // the clusters of the merged set, as TopPIC_postprocess builds them
    std::string db_file_name =
        ori_db_file_name + "_idx" + file_util::getFileSeparator() +
        file_util::filenameFromEntirePath(arguments["databaseFileName"]);
    bool is_proteoform_ppm_error = (arguments["proteoformPpmError"] == "true");
    double proteoform_error_tole =
        std::stod(arguments["proteoformErrorTolerance"]);
    if (is_proteoform_ppm_error) {
      proteoform_error_tole = proteoform_error_tole / 1000000;
    }
    if (arguments["useFeatureFile"] == "true") {
      prsm_feature_cluster::process(
          post_sp_name, "toppic_second_merged", "toppic_second_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    } else {
      prsm_simple_cluster::process(
          db_file_name, post_sp_name, "toppic_second_merged",
          prsm_para_ptr->getFixModPtrVec(), "toppic_second_cluster",
          is_proteoform_ppm_error, proteoform_error_tole);
    }
    prsm_prot_cluster::process(post_sp_name, "toppic_second_cluster",
                               "toppic_second_prot");
    std::string cur_suffix = "toppic_second_prot";

    bool use_rescore = arguments["rescore"] == "true";
    if (arguments["searchType"] == "TARGET+DECOY" && use_rescore) {
      std::cout << "Rescoring - started." << std::endl;
      // the second units with their residual spectra, the masses their
      // search used
      prsm_rescore::process(prsm_para_ptr, cur_suffix,
                            "toppic_second_rescore", resid_post_sp_name,
                            getRescoreEValueFloor(arguments));
      std::cout << "Rescoring - finished." << std::endl;
      cur_suffix = "toppic_second_rescore";
    }
    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      if (arguments["useFdrGroups"] == "true") {
        prsm_fdr_groups::process(post_sp_name, cur_suffix,
                                 "toppic_second_fdr");
      } else {
        prsm_fdr::process(post_sp_name, cur_suffix, "toppic_second_fdr",
                          use_rescore, getFdrCutoffs(arguments),
                          arguments["rescoreFallback"] != "false");
      }
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = "toppic_second_fdr";
    }
    file_util::copyFile(post_base + "." + cur_suffix,
                        post_base + ".toppic_record_prsm", true);

    // the output spectra: the first pass's post-matched file with the
    // masses the second pass's post mass matching added (not with -D, when
    // the spectrum file is the input file and nothing was added)
    if (resid_post_sp_name != resid_sp_name) {
      std::string merged_tmp_name = post_sp_name + "_second_tmp";
      prsm_second_prec::mergePostMatchedMasses(post_sp_name, resid_post_sp_name,
                                               merged_tmp_name,
                                               prsm_para_ptr->getSpParaPtr());
      std::filesystem::rename(merged_tmp_name, post_sp_name);
      msalign_util::geneSpIndex(post_sp_name);
    }

    // the residual pass's files: the intermediates as the first pass's,
    // and without -k the residual spectra too
    std::string resid_post_arg =
        resid_post_sp_name == resid_sp_name ? "" : resid_post_sp_name;
    cleanToppicDir(ori_db_file_name, resid_sp_name, keep_temp_files,
                   resid_post_arg, false);
    if (!keep_temp_files) {
      file_util::delFile(resid_sp_name);
      file_util::delFile(resid_base + ".feature");
      file_util::delFile(resid_sp_name + "_index");
      if (resid_post_arg != "") {
        file_util::delFile(resid_post_sp_name);
        file_util::delFile(resid_post_base + ".feature");
        file_util::delFile(resid_post_sp_name + "_index");
      }
    }
    std::cout << "Second precursor search - finished." << std::endl;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    return 1;
  }
  return 0;
}

// The two-pass m/z recalibration of one spectrum file (hidden --mass-calib;
// ms/calib, prsm/prsm_mass_calib; toppic_test experiment 152, amendment 3).
// A zero-shift calibration search of arguments["spectrumFileName"] on its own
// copy of the arguments: no unexpected shift, no variable PTM, no tag filter
// or TopMG stage, and precursor and fragment tolerances of at least 10 ppm
// (min_tole), so that an instrument's offset of a few ppm does not cost the
// calibrants their fragments and the index of the default tolerance serves.
// Its target PrSMs with E-value <= 1e-5, no unexpected shift and no variable
// PTM are the calibrants: the MS1 model is fitted on their precursor errors,
// the MS/MS model on their fragment errors (mass_calib_fit::fitModel), and
// written to <base>_mass_calib.json, whose name the post mass matching gets
// (arguments["massCalibModelFileName"]) to correct TopFD's peaks; the
// corrected spectra are written to <base>_calib_ms2.msalign, with a copy of
// the feature file, and TopFD's database arguments["topfdDbFileName"], with
// every m/z and deconvoluted mass corrected, to toppic's database
// arguments["toppicDbFileName"] (<base>_toppic.sqlite; sqlite_mass_calib;
// toppic_test experiments 154, 156 and 168: the one case in which toppic's
// database holds every table), to which the search of the corrected spectra
// writes (the post-matched envelopes, the PrSM tables). The calibration
// search's files are deleted (kept with -k). Returns the corrected spectrum
// file name; with fewer than 3 precursor or 4 fragment calibrants the file
// is searched uncorrected: its own name, and toppic's database has TopFD's
// envelope tables only, as without --mass-calib. "" on an error.
std::string TopPIC_mass_calib(std::map<std::string, std::string>& arguments) {
  std::string sp_file_name = arguments["spectrumFileName"];
  try {
    MassCalibPara para = getMassCalibPara();
    std::cout << "Mass calibration of " << sp_file_name << " - started."
              << std::endl;
    std::map<std::string, std::string> calib_arguments = arguments;
    calib_arguments["shiftNumber"] = "0";
    calib_arguments["variablePtmFileName"] = "";
    calib_arguments["tagFilterStages"] = "";
    calib_arguments["topmgMode"] = "disabled";
    calib_arguments["deamidSelector"] = "false";
    calib_arguments["massCalib"] = "false";
    if (std::stod(arguments["fragErrorTolerance"]) < para.min_tole) {
      calib_arguments["fragErrorTolerance"] = toleranceStr(para.min_tole);
    }
    if (arguments["precDaltonError"] == "" &&
        std::stod(arguments["precErrorTolerance"]) < para.min_tole) {
      calib_arguments["precErrorTolerance"] = toleranceStr(para.min_tole);
    }
    std::cout << "Mass calibration - zero unexpected shift search (fragment "
                 "tolerance "
              << calib_arguments["fragErrorTolerance"] << " ppm) - started."
              << std::endl;
    std::string raw_prsm_ext = TopPIC_identify(calib_arguments);
    if (raw_prsm_ext.empty()) {
      return "";
    }
    std::cout << "Mass calibration - zero unexpected shift search - finished."
              << std::endl;
    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(calib_arguments);
    MassCalibrantSet set = prsm_mass_calib::readCalibrants(
        prsm_para_ptr, raw_prsm_ext, para.max_evalue, para.tol_ppm);
    std::cout << "Mass calibration - calibrants: " << set.used_prsm_num
              << " of " << set.prsm_num
              << " PrSMs (targets with E-value <= " << para.max_evalue
              << ", no unexpected shift and no variable PTM), "
              << set.prec.size() << " precursors and " << set.frag.size()
              << " fragments." << std::endl;
    cleanToppicDir(arguments["oriDatabaseFileName"], sp_file_name,
                   arguments["keepTempFiles"] == "true", "", false);
    std::string topfd_db_name = arguments["topfdDbFileName"];
    std::string toppic_db_name = arguments["toppicDbFileName"];
    bool no_sql = arguments["noSql"] == "true";
    // the post mass matching corrects TopFD's peaks with the model file;
    // none without a model
    arguments["massCalibModelFileName"] = "";
    if (set.prec.size() < 3 || set.frag.size() < 4) {
      std::cout << "WARNING: too few calibrants for the mass calibration; "
                << sp_file_name << " is searched without it." << std::endl;
      if (!no_sql) {
        sql_db_copy::writeEnvelopeDb(topfd_db_name, toppic_db_name);
      }
      return sp_file_name;
    }
    MassCalibModel model = mass_calib_fit::fitModel(set);
    std::ostringstream info;
    info << std::setprecision(3) << "Mass calibration - model: MS1 "
         << model.ms1.c0 << " ppm, " << model.ms1.c_rt_per_hour
         << " ppm/h; MS/MS " << model.ms2.c0 << " ppm, " << model.ms2.c_mz
         << " ppm per 1000 m/z, " << model.ms2.c_rt_per_hour
         << " ppm/h; robust SD of the residuals: precursors "
         << model.prec_robust_sd << " ppm, fragments " << model.frag_robust_sd
         << " ppm.";
    std::cout << info.str() << std::endl;
    std::string model_file_name = prsm_mass_calib::modelFileName(sp_file_name);
    model.write(model_file_name);
    std::string calib_sp_file_name =
        prsm_mass_calib::calibSpFileName(sp_file_name);
    int sp_num =
        msalign_mass_calib::process(model, sp_file_name, calib_sp_file_name);
    if (arguments["useFeatureFile"] == "true") {
      file_util::copyFile(file_util::basename(sp_file_name) + ".feature",
                          file_util::basename(calib_sp_file_name) + ".feature",
                          true);
    }
    // TopFD's database, read only, and toppic's corrected copy of it for the
    // search of the corrected spectra (TopFD's masses are monoisotopic, also
    // after -A); with --no-sql no copy
    arguments["massCalibModelFileName"] = model_file_name;
    if (!no_sql) {
      if (std::filesystem::exists(topfd_db_name)) {
        sqlite_mass_calib::process(model, topfd_db_name, toppic_db_name);
      } else {
        std::cout << "Mass calibration - the database " << topfd_db_name
                  << " was not found: no " << toppic_db_name << " is written."
                  << std::endl;
      }
    }
    std::cout << "Mass calibration - finished: " << sp_num
              << " spectra corrected in " << calib_sp_file_name
              << ", the model in " << model_file_name << "." << std::endl;
    return calib_sp_file_name;
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
    return "";
  }
}

int TopPICProgress(std::map<std::string, std::string>& arguments) {
  std::string ori_sp_file_name = arguments["spectrumFileName"];
  std::string raw_prsm_ext = TopPIC_identifyAll(arguments);
  if (raw_prsm_ext.empty()) {
    return 1;
  }
  TopPIC_postprocess(arguments, raw_prsm_ext);
  if (arguments["secondPrecursor"] == "true") {
    if (TopPIC_second_precursor(arguments, ori_sp_file_name) != 0) {
      return 1;
    }
  }
  // In a combined run (-c) the file is reported after the combined result,
  // with its proteoform ids (TopPICProgress_multi_file).
  if (arguments["combinedOutputName"] == "") {
    TopPIC_report(arguments);
  }
  return 0;
}

int TopPICProgress_multi_file(std::map<std::string, std::string>& arguments,
                              const std::vector<std::string>& ori_spec_file_lst) {
  std::vector<std::string> spec_file_lst = ori_spec_file_lst;
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

  // Average-mass spectrum files (topfd --average-mass): every step below,
  // including the merge of a combined run, works on the converted
  // <name>_mono_ms2.msalign files with monoisotopic masses.
  if (arguments["averageMass"] == "true") {
    std::cout << "Converting average masses to monoisotopic masses - started."
              << std::endl;
    EnvBase::initBase(arguments["resourceDir"]);
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      spec_file_lst[k] = avg_mass_util::convertToMonoMass(spec_file_lst[k]);
      std::cout << "Converted " << ori_spec_file_lst[k] << " to "
                << spec_file_lst[k] << "." << std::endl;
    }
    std::cout << "Converting average masses to monoisotopic masses - finished."
              << std::endl;
  }

  // toppic's database of each spectrum file (toppic_test experiments 156 and
  // 168): <name>_toppic.sqlite for <name>_ms2.msalign and TopFD's
  // <name>.sqlite, whatever the options. Before the file is searched an
  // earlier one is deleted and a new one written with TopFD's tables ms2_env
  // and ms2_env_peak (with --mass-calib, when it fits a model, the corrected
  // copy of every table, by TopPIC_mass_calib); the post mass matching, which
  // reads its peaks from TopFD's database, appends its envelopes to it and
  // the report writes toppic's tables into it, so that TopFD's database is
  // only read. A reader takes the tables toppic's database lacks from
  // TopFD's. Every step gets the names from arguments["topfdDbFileName"] and
  // arguments["toppicDbFileName"]. With -C (no search) the earlier run's
  // database is kept: it holds that run's post-matched envelopes.
  std::vector<std::string> topfd_db_lst(spec_file_lst.size());
  std::vector<std::string> toppic_db_lst(spec_file_lst.size());
  for (size_t k = 0; k < spec_file_lst.size(); k++) {
    topfd_db_lst[k] = prsm_post_mass_match::sqlFileName(ori_spec_file_lst[k]);
    toppic_db_lst[k] =
        prsm_post_mass_match::toppicSqlFileName(ori_spec_file_lst[k]);
    if (combine_result_only) {
      continue;
    }
    if (sql_db_copy::removeDb(toppic_db_lst[k])) {
      std::cout << "Deleted the earlier database " << toppic_db_lst[k] << "."
                << std::endl;
    }
    // --no-sql: no database is written (an earlier one is still deleted, so
    // that none of another run sits next to the new results)
    if (arguments["massCalib"] != "true" && arguments["noSql"] != "true") {
      // TopFD's tables ms2_env and ms2_env_peak (toppic_test experiment 168);
      // a reader takes the other tables from TopFD's database
      sql_db_copy::writeEnvelopeDb(topfd_db_lst[k], toppic_db_lst[k]);
    }
  }
  // --mass-calib: the calibration model of each file, which the post mass
  // matching applies to TopFD's peaks
  std::vector<std::string> calib_model_lst(spec_file_lst.size());

  // --mass-calib: every file is recalibrated before the searches, and every
  // step below, including the merge of a combined run, works on the
  // corrected <name>_calib_ms2.msalign files (after -A
  // <name>_mono_calib_ms2.msalign; with -C, on the earlier run's) and the
  // corrected <name>_toppic.sqlite
  if (arguments["massCalib"] == "true") {
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      if (combine_result_only) {
        std::string calib_name =
            prsm_mass_calib::calibSpFileName(spec_file_lst[k]);
        if (std::filesystem::exists(calib_name)) {
          spec_file_lst[k] = calib_name;
        }
        continue;
      }
      arguments["spectrumFileName"] = spec_file_lst[k];
      arguments["topfdDbFileName"] = topfd_db_lst[k];
      arguments["toppicDbFileName"] = toppic_db_lst[k];
      std::string calib_name = TopPIC_mass_calib(arguments);
      if (calib_name.empty()) {
        return 1;
      }
      spec_file_lst[k] = calib_name;
      calib_model_lst[k] = arguments["massCalibModelFileName"];
    }
  }

  if (!combine_result_only) {
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
      std::string start_time = buf;
      arguments["startTime"] = start_time;
      arguments["spectrumFileName"] = spec_file_lst[k];
      arguments["topfdDbFileName"] = topfd_db_lst[k];
      arguments["toppicDbFileName"] = toppic_db_lst[k];
      arguments["massCalibModelFileName"] = calib_model_lst[k];
      if (toppic::TopPICProgress(arguments) != 0) {
        return 1;
      }
      // with -c the file is reported, and then cleaned up, after the
      // combined result (below): the report needs the spectra and their
      // index
      if (arguments["combinedOutputName"] == "") {
        // the post mass matching replaces the spectrum file name
        std::string post_sp_name = "";
        if (arguments["spectrumFileName"] != spec_file_lst[k]) {
          post_sp_name = arguments["spectrumFileName"];
        }
        cleanToppicDir(ori_db_file_name, spec_file_lst[k], keep_temp_files,
                       post_sp_name);
      }
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
    // The record file of each fraction and the spectra it was matched
    // against are named after the post-matched spectrum file,
    // <frac>_post_ms2 (prsm_post_mass_match::postBaseName: <frac> the base
    // of the spectrum file, with the _mono suffix of a -g run and the _calib
    // suffix of --mass-calib), or after the spectrum file with -D,
    // <frac>_ms2.
    bool post_matched = arguments["postMassMatch"] == "true";
    std::string ms2_suffix =
        post_matched ? "_post_ms2.msalign" : "_ms2.msalign";
    std::vector<std::string> record_base_list(raw_file_list.size());
    std::vector<std::string> prsm_file_lst(raw_file_list.size());
    // the MS1 spectra are TopFD's, <frac>_ms1.msalign without toppic's
    // suffixes
    std::vector<std::string> ms1_base_list(raw_file_list.size());
    for (size_t i = 0; i < raw_file_list.size(); i++) {
      ms1_base_list[i] = prsm_post_mass_match::topfdBaseName(spec_file_lst[i]);
      std::string record_base = raw_file_list[i];
      if (post_matched) {
        std::string post_base =
            prsm_post_mass_match::postBaseName(spec_file_lst[i]);
        record_base = post_base.substr(0, post_base.size() - 9);
      }
      record_base_list[i] = record_base;
      record_base += post_matched ? "_post_ms2" : "_ms2";
      prsm_file_lst[i] = record_base + ".toppic_record_prsm";
    }

    std::cout << "Merging files - started." << std::endl;
    // The merged spectra (<name>_ms2.msalign, with the masses the post mass
    // matching added) are what the report attaches to the merged PrSMs.
    std::cout << "Merging msalign files started." << std::endl;
    std::string para_str = "";
    msalign_frac_merge::mergeFractions(record_base_list, merged_file_name,
                                       para_str, ms2_suffix, ms1_base_list);
    std::cout << "Merging msalign files finished." << std::endl;

    std::cout << "Merging identification files started." << std::endl;
    prsm_util::mergePrsmFiles(prsm_file_lst, SpPara::getMaxSpecNumPerFile(),
                              SpPara::getMaxFeatureNumPerFile(),
                              merged_file_name + "_ms2.toppic_raw_prsm");
    std::cout << "Merging identification files finished." << std::endl;
    std::cout << "Merging files - finished." << std::endl;

    std::string sp_file_name = merged_file_name + "_ms2.msalign";
    arguments["spectrumFileName"] = sp_file_name;
    arguments["startTime"] = combined_start_time;
    // the combined result has no database
    arguments["toppicDbFileName"] = "";

    TopPIC_multi_postprocess(arguments);
    TopPIC_report(arguments);
    cleanToppicDir(ori_db_file_name, sp_file_name, keep_temp_files, "");

    // The files' reports. The PrSMs of each file's record file get the
    // proteoform ids of the combined result, so that a proteoform has one
    // id in the combined tables and in every file's tables; then each file
    // is reported and cleaned up as in a run without -c (with -C, from the
    // record files of the earlier run). The files' FDRs, protein cluster
    // ids and everything else in the record files are kept.
    std::cout << "Updating the proteoform ids of the files - started."
              << std::endl;
    prsm_util::updateProteoClusterIds(prsm_file_lst,
                                      SpPara::getMaxSpecNumPerFile(),
                                      merged_file_name +
                                          "_ms2.toppic_record_prsm");
    std::cout << "Updating the proteoform ids of the files - finished."
              << std::endl;
    for (size_t k = 0; k < spec_file_lst.size(); k++) {
      // the spectra the file's PrSMs were matched against: the post-matched
      // ones (<frac>_post_ms2.msalign), or the file itself with -D
      std::string file_sp_name =
          post_matched ? record_base_list[k] + ms2_suffix : spec_file_lst[k];
      std::cout << "Reporting " << file_sp_name << " - started." << std::endl;
      arguments["spectrumFileName"] = file_sp_name;
      arguments["startTime"] = combined_start_time;
      arguments["toppicDbFileName"] = toppic_db_lst[k];
      msalign_util::geneSpIndex(file_sp_name);
      TopPIC_report(arguments);
      cleanToppicDir(ori_db_file_name, spec_file_lst[k], keep_temp_files,
                     post_matched ? file_sp_name : "");
      std::cout << "Reporting " << file_sp_name << " - finished." << std::endl;
    }
  }

  // each file's database compressed in place when complete, after the
  // files' reports of a combined run (toppic_test experiment 172), unless
  // --no-sql or the hidden --disable-zstd; with -C the earlier run's
  // database may be compressed already, and stays so (the PrSM tables the
  // report wrote into it again are plain tables)
  if (arguments["noSql"] != "true" && arguments["disableZstd"] != "true") {
    for (const std::string& name : toppic_db_lst) {
      zstd_db::compressInPlace(name);
    }
  }

  base_data::release();

  std::cout << "TopPIC finished." << std::endl << std::flush;

  return 0;
}

}  // namespace toppic
