//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include "console/topmerge_argument.hpp"
#include "console/topmerge_process.hpp"

#include "common/base/base_data.hpp"
#include "common/util/file_util.hpp"
#include "common/util/version.hpp"

#include "ms/spec/msalign_util.hpp"
#include "ms/spec/msalign_frac_merge.hpp"
#include "ms/spec/deconv_json_merge.hpp"
#include "ms/feature/feature_merge.hpp"

#include "para/prsm_para.hpp"

#include "para/prsm_para.hpp"
#include "prsm/prsm_str_merge.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_simple_cluster.hpp"
#include "prsm/prsm_feature_cluster.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"
#include "prsm/prsm_form_filter.hpp"
#include "prsm/prsm_util.hpp"

#include "stat/local/local_mng.hpp"
#include "stat/local/local_processor.hpp"

#include "visual/xml_generator.hpp"
#include "visual/json_transformer.hpp"

namespace toppic{
void copyTopMSV(std::map<std::string, std::string> &arguments) {
  std::string spectrum_file_name = arguments["spectrumFileName"];
  std::string base_name = file_util::basename(spectrum_file_name);
  std::string base_name_short = base_name.substr(0, base_name.length() - 4);
  std::string topmsv_dir = base_name_short + "_html" +  file_util::getFileSeparator() + "topmsv";
  if (file_util::exists(topmsv_dir)) {
    LOG_WARN("The TopMSV directory " << topmsv_dir << " exists!");
    //file_util::delDir(topmsv_dir);
  }
  else{
    if (!file_util::exists(base_name_short + "_html")){//if _html folder was not created with topfd
      file_util::createFolder(base_name_short + "_html");
    }
    std::string resource_dir = arguments["resourceDir"];
    // copy resources 
    std::string from_path(resource_dir + file_util::getFileSeparator() + "topmsv");
    file_util::copyDir(from_path, topmsv_dir);
  }
}
void cleanTopMergeDir(const std::string &fa_name, 
                    const std::string & sp_name,
                    bool keep_temp_files,
                    std::string tool_name) {
  //std::string fa_base = file_util::absoluteName(fa_name);
  //std::replace(fa_base.begin(), fa_base.end(), '\\', '/');
  std::string abs_sp_name = file_util::absoluteName(sp_name);
  std::string sp_base = file_util::basename(abs_sp_name);
  std::replace(sp_base.begin(), sp_base.end(), '\\', '/');
  file_util::rename(sp_base + "." + tool_name + "_form_cutoff_form",
                    sp_base + "_" + tool_name + "_proteoform.xml");
  file_util::rename(sp_base + "." + tool_name + "_prsm",
                    sp_base + "_" + tool_name + "_prsm.xml");
  if (!keep_temp_files) {
    //file_util::cleanPrefix(fa_name, fa_base + "_");
    file_util::cleanPrefix(sp_name, sp_base + ".msalign_");
    file_util::delFile(abs_sp_name + "_index");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_zero_filter_");
    file_util::delFile(sp_base + "." + tool_name + "_zero_ptm");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_zero_ptm_");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_one_filter_");
    file_util::delFile(sp_base + "." + tool_name + "_one_ptm");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_one_ptm_");
    file_util::delFile(sp_base + "." + tool_name + "_multi_filter");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_multi_filter_");
    file_util::delFile(sp_base + "." + tool_name + "_multi_ptm");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_multi_ptm_");
    file_util::delFile(sp_base + "." + tool_name + "_combined");
    file_util::delFile(sp_base + "." + tool_name + "_evalue");
    file_util::cleanPrefix(sp_name, sp_base + "." + tool_name + "_evalue_");
    //file_util::delFile(sp_base + ".toppic_top");
    file_util::delFile(sp_base + "." + tool_name + "_cluster");
    file_util::delFile(sp_base + "." + tool_name + "_cluster_fdr");
    file_util::delFile(sp_base + "." + tool_name + "_prsm_cutoff");
    file_util::delFile(sp_base + "." + tool_name + "_prsm_cutoff_local");
    file_util::delFile(sp_base + "." + tool_name + "_form_cutoff");
    file_util::delDir(sp_base + "_" + tool_name + "_proteoform_cutoff_xml");
    file_util::delDir(sp_base + "_" + tool_name + "_prsm_cutoff_xml");
  }
}

void getToolNames(std::vector<std::string> &tool_vec, std::string base_path, const std::vector<std::string> & spec_file_lst) {
  bool is_toppic_data = true;
  bool is_topmg_data = true;
  std::cout << file_util::basename(spec_file_lst[0]) << std::endl;
  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    if (!file_util::exists(file_util::basename(spec_file_lst[i]) + "_toppic_prsm.tsv")) {
      is_toppic_data = false;
    }
  }
  for (size_t i = 0; i < spec_file_lst.size(); i++) {
    if (!file_util::exists(file_util::basename(spec_file_lst[i]) + "_topmg_prsm.tsv")) {
      is_topmg_data = false;
    }
  }
  if (is_toppic_data) {
    tool_vec.push_back("toppic");
  }
  if (is_topmg_data) {
    tool_vec.push_back("topmg");
  }
}

// proteoform clustering + FDR + HTML generation
int TopMerge_post(std::map<std::string, std::string> & arguments) {
  try {
    std::string resource_dir = arguments["resourceDir"];

    base_data::init();
    LOG_DEBUG("Init base data completed");

    std::string sp_file_name = arguments["spectrumFileName"];
    std::string ori_db_file_name = arguments["oriDatabaseFileName"];
    std::string db_file_name = ori_db_file_name + "_idx" + file_util::getFileSeparator() 
      + file_util::filenameFromEntirePath(arguments["databaseFileName"]);
    double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
    double min_ptm_mass = std::stod(arguments["minPtmMass"]);
    bool localization = false;
    if (arguments["residueModFileName"] != "" && arguments["toolName"] == "toppic") {
      localization = true;
    }

    PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);
    msalign_util::geneSpIndex(sp_file_name);
    LOG_DEBUG("prsm para inited");

    std::cout << "Finding PrSM clusters - started." << std::endl;
    double form_error_tole = std::stod(arguments["proteoformErrorTolerance"]);
    LOG_DEBUG("form error tole " << form_error_tole);
    std::string id_file_name = "";

    if (arguments["toolName"] == "toppic") {
      id_file_name = "toppic_prsm";
    }
    else if (arguments["toolName"] == "topmg") {
      id_file_name = "topmg_prsm";
    }

    if (arguments["useFeatureFile"] == "true") {
      // TopFD msalign file with feature ID
      ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
      prsm_feature_cluster::process(sp_file_name,
                                    id_file_name,
                                    arguments["toolName"] + "_cluster",
                                    form_error_tole);
    } 
    else {
      prsm_simple_cluster::process(db_file_name, sp_file_name,
                                   id_file_name, prsm_para_ptr->getFixModPtrVec(),
                                   arguments["toolName"] + "_cluster", form_error_tole);
    }
    std::cout << "Finding PrSM clusters - finished." << std::endl;
    std::string cur_suffix = arguments["toolName"] + "_cluster";

    if (arguments["searchType"] == "TARGET+DECOY") {
      std::cout << "FDR computation - started. " << std::endl;
      prsm_fdr::process(sp_file_name, arguments["toolName"] + "_cluster", arguments["toolName"] + "_cluster_fdr", arguments["keepDecoyResults"]);
      std::cout << "FDR computation - finished." << std::endl;
      cur_suffix = arguments["toolName"] + "_cluster_fdr";
    }

    std::string cutoff_type = arguments["cutoffSpectralType"];
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    double cutoff_value;
    std::istringstream(arguments["cutoffSpectralValue"]) >> cutoff_value;
    prsm_cutoff_selector::process(db_file_name, sp_file_name, cur_suffix,
                                  arguments["toolName"] + "_prsm_cutoff", cutoff_type, cutoff_value);
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;
    cur_suffix = arguments["toolName"] + "_prsm_cutoff";

    if (localization) {
      std::cout << "PTM characterization - started." << std::endl;
      LocalMngPtr local_mng
          = std::make_shared<LocalMng>(prsm_para_ptr,
                                       std::stod(arguments["localThreshold"]),
                                       arguments["residueModFileName"],
                                       min_ptm_mass,
                                       max_ptm_mass,
                                       arguments["toolName"] + "_prsm_cutoff", 
                                       arguments["toolName"] + "_prsm_cutoff_local");
      LocalProcessorPtr local_ptr = std::make_shared<LocalProcessor>(local_mng);
      local_ptr->process();
      local_ptr = nullptr;
      std::cout << "PTM characterization - finished." << std::endl;
      cur_suffix = arguments["toolName"] + "_prsm_cutoff_local";
    }


    std::time_t end = time(nullptr);
    char buf[50];
    std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&end));
    arguments["endTime"] = buf;

    std::string argu_str = Argument::outputTsvArguments(arguments);

    std::cout << "Outputting PrSM table - started." << std::endl;
    PrsmTableWriterPtr table_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str, cur_suffix, "_" + arguments["toolName"] + "_prsm.tsv");
    table_out->write();
    table_out = nullptr;
    std::cout << "Outputting PrSM table - finished." << std::endl;

    XmlGeneratorPtr xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, 
                                                                cur_suffix, arguments["toolName"] + "_prsm_cutoff");
    if (arguments["geneHTMLFolder"] == "true"){
      std::cout << "Generating PrSM xml files - started." << std::endl;
    
      xml_gene->process();
      xml_gene = nullptr;
      std::cout << "Generating PrSM xml files - finished." << std::endl;

      copyTopMSV(arguments);
  
      std::cout << "Converting PrSM xml files to json files - started." << std::endl;
      jsonTranslate(arguments, arguments["toolName"] + "_prsm_cutoff");
      std::cout << "Converting PrSM xml files to json files - finished." << std::endl;
    }

    cutoff_type = (arguments["cutoffProteoformType"] == "FDR") ? "FORMFDR": "EVALUE";
    std::cout << "PrSM filtering by " << cutoff_type << " - started." << std::endl;
    std::istringstream(arguments["cutoffProteoformValue"]) >> cutoff_value;
    prsm_cutoff_selector::process(db_file_name, sp_file_name, cur_suffix,
                                  arguments["toolName"]+ "_form_cutoff", cutoff_type,
                                  cutoff_value);
    std::cout << "PrSM filtering by " << cutoff_type << " - finished." << std::endl;

    std::cout << "Selecting top PrSMs for proteoforms - started." << std::endl;
    prsm_form_filter::process(db_file_name, sp_file_name, arguments["toolName"] + "_form_cutoff",
                              arguments["toolName"] + "_form_cutoff_form");
    std::cout << "Selecting top PrSMs for proteoforms - finished." << std::endl;

    std::cout << "Outputting proteoform table - started." << std::endl;
    PrsmTableWriterPtr form_out
        = std::make_shared<PrsmTableWriter>(prsm_para_ptr, argu_str,
                                            arguments["toolName"] + "_form_cutoff_form", "_" + arguments["toolName"] + "_proteoform.tsv");
    form_out->write();
    form_out = nullptr;
    std::cout << "Outputting proteoform table - finished." << std::endl;

    if (arguments["geneHTMLFolder"] == "true"){

      std::cout << "Generating proteoform xml files - started." << std::endl;
      xml_gene = std::make_shared<XmlGenerator>(prsm_para_ptr, resource_dir, 
                                              arguments["toolName"] + "_form_cutoff", 
                                              arguments["toolName"] + "_proteoform_cutoff");
    
      xml_gene->process();
      xml_gene = nullptr;
      std::cout << "Generating proteoform xml files - finished." << std::endl;

    
      std::cout << "Converting proteoform xml files to html files - started." << std::endl;
      jsonTranslate(arguments, arguments["toolName"] + "_proteoform_cutoff");
      std::cout << "Converting proteoform xml files to html files - finished." << std::endl;
    }
  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  return 0;
}

void TopMergeProcess(std::map<std::string, std::string> & arguments,
                     const std::vector<std::string> & spec_file_lst){
  std::vector<std::string> tool_vec;
  std::string base_path = file_util::absoluteDir(spec_file_lst[0]);

  getToolNames(tool_vec, base_path, spec_file_lst);

  for (auto it = std::begin(tool_vec); it != std::end(tool_vec); it++) {
    std::cout << *it << std::endl;
  }

  for (auto it = std::begin(tool_vec); it != std::end(tool_vec); it++) {
    std::string tool_name = *it;
    std::string input_prsm_name = "";
    std::string output_prsm_name = "";
    try {
      std::cout << "TopMerge " << Version::getVersion() << std::endl;
      arguments["version"] = toppic::Version::getVersion();

      xercesc::XMLPlatformUtils::Initialize(); 

      std::string full_combined_name = base_path + file_util::getFileSeparator() 
        +  arguments["combinedOutputName"];
    
      std::time_t start = time(nullptr);
      char buf[50];
      std::strftime(buf, 50, "%a %b %d %H:%M:%S %Y", std::localtime(&start));
      arguments["startTime"] = buf;
      arguments["spectrumFileName"] = ""; //reset after the processing for the first app output finishes
      for (size_t i = 0; i < spec_file_lst.size(); i++) {
        arguments["spectrumFileName"] = arguments["spectrumFileName"] + spec_file_lst[i] + " ";
      }    
      arguments["toolName"] = tool_name;

      if (arguments["toolName"] == "toppic") {
        input_prsm_name = "toppic_prsm.xml";
        output_prsm_name = "toppic_prsm";
      }
      else if (arguments["toolName"] == "topmg") {
        input_prsm_name = "topmg_prsm.xml";
        output_prsm_name = "topmg_prsm";
      }

      Argument::outputArguments(std::cout, arguments);

      base_data::init();

      PrsmParaPtr prsm_para_ptr = std::make_shared<PrsmPara>(arguments);

      if (arguments["combinedOutputName"] != "") {
        std::string merged_file_name = arguments["combinedOutputName"]; 
        std::string para_str = "";
        std::cout << "Merging files started." << std::endl;
        std::cout << "Merging msalign files started." << std::endl;
        MsAlignFracMerge::mergeFiles(spec_file_lst, full_combined_name + "_ms2.msalign", para_str);
        std::cout << "Merging msalign files finished." << std::endl;
        if (arguments["geneHTMLFolder"] == "true"){
          std::cout << "Merging json files started." << std::endl;
          DeconvJsonMergePtr json_merger 
              = std::make_shared<DeconvJsonMerge>(spec_file_lst, full_combined_name);
          json_merger->process();
          json_merger = nullptr;
          std::cout << "Merging json files finished." << std::endl;
        }
        if (arguments["useFeatureFile"] == "true") {//only when feature files are being used
          std::cout << "Merging feature files started." << std::endl;
          feature_merge::process(spec_file_lst, full_combined_name, para_str);
          std::cout << "Merging feature files finished." << std::endl;
        }
        // merge TOP files
        std::cout << "Merging identification files started." << std::endl;
        std::vector<std::string> prsm_file_lst(spec_file_lst.size());
        for (size_t i = 0; i < spec_file_lst.size(); i++) {
          prsm_file_lst[i] = file_util::basename(spec_file_lst[i]) + "_" + input_prsm_name; 
        }
        //toppic_prsm 
        int N = 1000000;
        prsm_util::mergePrsmFiles(prsm_file_lst, N , full_combined_name + "_ms2." + output_prsm_name);
        std::cout << "Merging identification files finished." << std::endl;
        std::cout << "Merging files - finished." << std::endl;

        std::string sp_file_name = full_combined_name + "_ms2.msalign";
        arguments["spectrumFileName"] = sp_file_name;

        TopMerge_post(arguments);
      }
      bool keep_temp_files = (arguments["keepTempFiles"] == "true"); 

      if (!keep_temp_files) {
        std::cout << "Deleting temporary files - started." << std::endl;
        std::string ori_db_file_name = arguments["oriDatabaseFileName"];

        /*for (size_t k = 0; k < spec_file_lst.size(); k++) {
          std::string sp_file_name = spec_file_lst[k];
          cleanTopMergeDir(ori_db_file_name, sp_file_name, keep_temp_files);
        }*/

        if (arguments["combinedOutputName"] != "") {
          std::string sp_file_name = full_combined_name + "_ms2.msalign";
          cleanTopMergeDir(ori_db_file_name, sp_file_name, keep_temp_files, arguments["toolName"]);
        }
        std::cout << "Deleting temporary files - finished." << std::endl; 
      }
      base_data::release();

      std::cout << "TopMerge - finished." << std::endl;
    } catch (const char* e) {
        std::cout << "[Exception]" << std::endl;
        std::cout << e << std::endl;
        exit(EXIT_FAILURE);
    }
  }
}
}
