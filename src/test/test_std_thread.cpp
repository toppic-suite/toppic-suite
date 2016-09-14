// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include <iostream>
#include <iomanip>
#include <thread>

#include "base/fasta_reader.hpp"
#include "base/fasta_util.hpp"
#include "base/base_data.hpp"
#include "base/web_logger.hpp"

#include "spec/msalign_reader.hpp"
#include "spec/msalign_util.hpp"

#include "prsm/prsm_para.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_top_selector.hpp"
#include "prsm/prsm_cutoff_selector.hpp"
#include "prsm/prsm_species.hpp"
#include "prsm/prsm_table_writer.hpp"
#include "prsm/prsm_fdr.hpp"

#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_filter_processor.hpp"

#include "zeroptmsearch/zero_ptm_search_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_processor.hpp"

#include "oneptmsearch/ptm_search_mng.hpp"
#include "oneptmsearch/one_ptm_search.hpp"

#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_processor.hpp"

#include "ptmsearch/ptm_search_processor.hpp"

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"

#include "local/local_mng.hpp"
#include "local/local_processor.hpp"

#include "prsmview/xml_generator.hpp"
#include "prsmview/transformer.hpp"

#include "console/argument.hpp"

namespace prot {

void spec_split(int nthread, const std::string & fname) {
  std::vector<std::string> split_files(nthread);
  std::vector<std::shared_ptr<std::ofstream>> out;
  for (size_t i = 0; i < split_files.size(); i++) {
    split_files[i] = FileUtil::basename(fname) + "_" + std::to_string(i) + ".msalign";
    out.push_back(std::make_shared<std::ofstream>(split_files[i]));
  }

  MsAlignReader msreader(fname, 1);
  std::vector<std::string> spectrum = msreader.readOneSpectrum();
  int count = 0;
  while (spectrum.size() != 0) {
    *out[count / nthread] << std::endl;
    for (size_t i = 0; i < spectrum.size(); i++) {
      *out[count / nthread] << spectrum[i] << std::endl;
    }
    *out[count / nthread] << std::endl;
    spectrum = msreader.readOneSpectrum();
  }
  for (size_t i = 0; i < out.size(); i++) {
    out[i]->close();
  }
}

void toppic_thread(int i, std::map<std::string, std::string> arguments) {

  std::string exe_dir = arguments["executiveDir"];
  LOG_DEBUG("Executive file directory is: " << exe_dir);

  BaseData::init(exe_dir);

  LOG_DEBUG("Init base data completed");
  arguments["spectrumFileName"] = FileUtil::basename(arguments["spectrumFileName"]) + "_" + std::to_string(i) + ".msalign";

  std::string db_file_name = arguments["databaseFileName"];
  std::string sp_file_name = arguments["spectrumFileName"];
  std::string ori_db_file_name = arguments["oriDatabaseFileName"];
  std::string log_file_name = arguments["logFileName"];

  int n_top = std::stoi(arguments["numOfTopPrsms"]);
  int ptm_num = std::stoi(arguments["ptmNumber"]);
  double max_ptm_mass = std::stod(arguments["maxPtmMass"]);
  int filter_result_num = std::stoi(arguments["filteringResultNumber"]);

  bool use_gf = false; 
  if (arguments["useGf"] == "true") {
    use_gf = true;
  }
  // initialize log file 
  WebLog::init(log_file_name, use_gf, ptm_num);
  LOG_DEBUG("web log inited");

  PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments));
  LOG_DEBUG("prsm para inited");

  bool decoy = false;
  if (arguments["searchType"] == "TARGET+DECOY") {
    decoy = true;
  }
  LOG_DEBUG("block size " << arguments["databaseBlockSize"]);
  int db_block_size = std::stoi(arguments["databaseBlockSize"]);

  FastaUtil::dbPreprocess (ori_db_file_name, db_file_name, decoy, db_block_size);
  MsAlignUtil::geneSpIndex(sp_file_name);

  std::vector<std::string> input_exts ;

  LOG_DEBUG("Zero PTM filtering started.");
  ZeroPtmFilterMngPtr zero_filter_mng_ptr = ZeroPtmFilterMngPtr(new ZeroPtmFilterMng (prsm_para_ptr, "ZERO_FILTER"));
  ZeroPtmFilterProcessorPtr zero_filter_processor = ZeroPtmFilterProcessorPtr(new ZeroPtmFilterProcessor(zero_filter_mng_ptr));
  zero_filter_processor->process();
  //WebLog::completeFunction(WebLog::ZeroPtmTime());
  LOG_DEBUG("Zero PTM filtering finished.");

  LOG_DEBUG("Zero PTM search started.");
  ZeroPtmSearchMngPtr zero_search_mng_ptr = ZeroPtmSearchMngPtr(new ZeroPtmSearchMng (prsm_para_ptr, "ZERO_FILTER", "ZERO_PTM"));
  ZeroPtmSearch::process(zero_search_mng_ptr);
  LOG_DEBUG("Zero PTM search finished.");

  input_exts.push_back("ZERO_PTM_COMPLETE");
  input_exts.push_back("ZERO_PTM_PREFIX");
  input_exts.push_back("ZERO_PTM_SUFFIX");
  input_exts.push_back("ZERO_PTM_INTERNAL");

  if (ptm_num >= 1) {

    LOG_DEBUG("One PTM filtering started.");
    OnePtmFilterMngPtr one_ptm_filter_mng_ptr = OnePtmFilterMngPtr(new OnePtmFilterMng (prsm_para_ptr, "ONE_PTM_FILTER"));
    OnePtmFilterProcessorPtr one_filter_processor = OnePtmFilterProcessorPtr(new OnePtmFilterProcessor(one_ptm_filter_mng_ptr));
    one_filter_processor->process();
    LOG_DEBUG("One PTM filtering finished.");

    LOG_DEBUG("One PTM search started.");
    int shift_num = 1;
    PtmSearchMngPtr one_search_mng_ptr 
        = PtmSearchMngPtr(new PtmSearchMng (prsm_para_ptr, n_top, max_ptm_mass, shift_num, "ONE_PTM_FILTER", "ONE_PTM"));
    OnePtmSearch::process(one_search_mng_ptr);
    LOG_DEBUG("One PTM search finished.");

    input_exts.push_back("ONE_PTM_COMPLETE");
    input_exts.push_back("ONE_PTM_PREFIX");
    input_exts.push_back("ONE_PTM_SUFFIX");
    input_exts.push_back("ONE_PTM_INTERNAL");
  }

  if (ptm_num >= 2) {
    LOG_DEBUG("Diagonal PTM filtering started.");
    DiagFilterMngPtr diag_filter_mng_ptr 
        = DiagFilterMngPtr(new DiagFilterMng (prsm_para_ptr, filter_result_num, "DIAG_FILTER"));
    DiagFilterProcessorPtr diag_filter_processor 
        = DiagFilterProcessorPtr(new DiagFilterProcessor(diag_filter_mng_ptr));
    diag_filter_processor->process();

    LOG_DEBUG("Diagonal filtering finished.");

    LOG_DEBUG("Two PTM search started.");
    PtmSearchMngPtr two_search_mng_ptr 
        = PtmSearchMngPtr(new PtmSearchMng (prsm_para_ptr, n_top, max_ptm_mass, ptm_num,
                                            "DIAG_FILTER", "PTM"));
    PtmSearchProcessorPtr processor = PtmSearchProcessorPtr(new PtmSearchProcessor(two_search_mng_ptr));
    processor->process();
    LOG_DEBUG("Two PTM search finished.");
    input_exts.push_back("PTM");
  }

  LOG_DEBUG("Combining PRSMs started.");
  int prsm_top_num = (ptm_num + 1) * 4;
  PrsmStrCombinePtr combine_ptr(new PrsmStrCombine(sp_file_name, input_exts, "RAW_RESULT", prsm_top_num));
  combine_ptr->process();
  combine_ptr = nullptr;
  LOG_DEBUG("Combining PRSMs finished.");

  LOG_DEBUG("E-value computation started.");
  bool variable_ptm = false;
  TdgfMngPtr tdgf_mng_ptr = TdgfMngPtr(new TdgfMng (prsm_para_ptr, ptm_num, max_ptm_mass, use_gf,
                                                    variable_ptm, "RAW_RESULT", "EVALUE"));
  EValueProcessorPtr processor = EValueProcessorPtr(new EValueProcessor(tdgf_mng_ptr));
  processor->init();
  // compute E-value for a set of prsms each run 
  processor->process(false);
  processor = nullptr;
  LOG_DEBUG("E-value computation finished.");

  if (arguments["searchType"]=="TARGET") { 
    LOG_DEBUG("Top PRSM selecting started");
    PrsmTopSelectorPtr selector = PrsmTopSelectorPtr(new PrsmTopSelector(db_file_name, sp_file_name, "EVALUE", "TOP", n_top));
    selector->process();
    selector = nullptr;
    LOG_DEBUG("Top PRSM selecting finished.");
  }
  else {
    std::cout << "Top PRSM selecting started " << std::endl;
    PrsmTopSelectorPtr selector = PrsmTopSelectorPtr(new PrsmTopSelector(db_file_name, sp_file_name, "EVALUE", "TOP_PRE", n_top));
    selector->process();
    selector = nullptr;
    LOG_DEBUG("Top PRSM selecting finished.");

    LOG_DEBUG("FDR computation started. ");
    PrsmFdrPtr fdr = PrsmFdrPtr(new PrsmFdr(db_file_name, sp_file_name, "TOP_PRE", "TOP"));
    fdr->process();
    fdr = nullptr;
    LOG_DEBUG("FDR computation finished.");
  }

  LOG_DEBUG("PRSM selecting by cutoff started.");
  std::string cutoff_type = arguments["cutoffType"];
  double cutoff_value;
  std::istringstream (arguments["cutoffValue"]) >> cutoff_value;
  PrsmCutoffSelectorPtr cutoff_selector = PrsmCutoffSelectorPtr(
      new PrsmCutoffSelector(db_file_name, sp_file_name, "TOP", "CUTOFF_RESULT", 
                             cutoff_type, cutoff_value));
  cutoff_selector->process();
  cutoff_selector = nullptr;
  LOG_DEBUG("PRSM selecting by cutoff finished.");

  LOG_DEBUG("Finding protein species started.");
  double ppo;
  std::istringstream(arguments["errorTolerance"]) >> ppo;
  ppo = ppo / 1000000.0;
  ModPtrVec fix_mod_list = prsm_para_ptr->getFixModPtrVec();
  PrsmSpeciesPtr prsm_species = PrsmSpeciesPtr(
      new PrsmSpecies(db_file_name, sp_file_name, "CUTOFF_RESULT",
                      fix_mod_list, "OUTPUT_RESULT",ppo));
  prsm_species->process();
  prsm_species = nullptr;
  LOG_DEBUG("Finding protein species finished.");

  LOG_DEBUG("PTM localization started.");
  LocalMngPtr local_mng = LocalMngPtr(
      new LocalMng(prsm_para_ptr, arguments["local_threshold"],
                   arguments["residueModFileName"], max_ptm_mass,
                   "OUTPUT_RESULT", "LOCAL_RESULT"));
  LocalProcessorPtr local_ptr = LocalProcessorPtr(new LocalProcessor(local_mng));
  local_ptr->process();
  local_ptr = nullptr;
  LOG_DEBUG("PTM localization finished.");

}

int toppic_combine(int num_threads, std::map<std::string, std::string> arguments) {
  try {

    PrsmParaPtr prsm_para_ptr = PrsmParaPtr(new PrsmPara(arguments)); 
    LOG_DEBUG("Outputting table starts ");
    PrsmTableWriterPtr table_out = PrsmTableWriterPtr(
        new PrsmTableWriter(prsm_para_ptr, "LOCAL_RESULT", "OUTPUT_TABLE"));
    table_out->write();
    table_out = nullptr;
    LOG_DEBUG("Outputting table finished.");

    LOG_DEBUG("Generating view xml files started.");
    XmlGeneratorPtr xml_gene = XmlGeneratorPtr(new XmlGenerator(prsm_para_ptr, arguments["executiveDir"], "LOCAL_RESULT"));
    xml_gene->process();
    xml_gene = nullptr;
    LOG_DEBUG("Generating view xml files finished.");

    LOG_DEBUG("Converting xml files to html files started.");
    translate(arguments);
    LOG_DEBUG("Converting xml files to html files finished.");

  } catch (const char* e) {
    std::cout << "[Exception]" << std::endl;
    std::cout << e << std::endl;
  }
  std::cout << "TopPC finished." << std::endl;
  return 0;

}

}

int main(int argc, char* argv[]) {
  prot::log_level = 2;
  int num_threads = 4;
  std::thread t[num_threads];
  std::cout << std::setprecision(10);

  prot::Argument argu_processor;
  bool success = argu_processor.parse(argc, argv);
  if (!success) {
    return 1;
  }
  std::map<std::string, std::string> arguments = argu_processor.getArguments();

  std::string ori_spec_filename = arguments["spectrumFileName"];

  prot::spec_split(num_threads, arguments["spectrumFileName"]);

  for (int i = 0; i < num_threads; ++i) 
    t[i] = std::thread(prot::toppic_thread, i, arguments);

  for (int i = 0; i < num_threads; ++i) 
    t[i].join();

  arguments["spectrumFileName"] = ori_spec_filename;

  return prot::toppic_combine(num_threads, arguments);
}
