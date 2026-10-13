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

#include "console/topmg_argument.hpp"

#include <exception>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
// #include <set>

#include "common/base/mod_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"
#include "para/sp_para.hpp"
#include "seq/fasta_util.hpp"

namespace toppic {

namespace {

// A built-in PTM's mass as ptm_base.xml gives it: 15 significant digits
// write each mass of the file back exactly
std::string ptmMassStr(const PtmPtr& ptm_ptr) {
  std::ostringstream stream;
  stream << std::setprecision(15) << ptm_ptr->getMonoMass();
  return stream.str();
}

}  // namespace

TopmgArgument::TopmgArgument() { arguments_ = initArguments(); }

std::map<std::string, std::string> TopmgArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["oriDatabaseFileName"] = "";
  arguments["databaseFileName"] = "";
  arguments["databaseBlockSize"] = "60000000";
  arguments["maxFragmentLength"] = "500";
  arguments["minBlockNum"] = "20";
  arguments["spectrumFileName"] = "";
  arguments["combinedOutputName"] = "";
  arguments["activation"] = "FILE";
  arguments["searchType"] = "TARGET";
  arguments["fixedMod"] = "";
  // the unexpected shifts of the graph alignment, TopMG's own number,
  // kept apart from toppic's shiftNumber (PrsmPara::getTopmgShiftNum)
  arguments["topmgShiftNumber"] = "1";
  arguments["precErrorTolerance"] = "10";
  arguments["fragErrorTolerance"] = "10";
  arguments["envCnnCutoff"] = "0.4";
  arguments["proteoformPpmError"] = "false";
  arguments["proteoformErrorTolerance"] = "1.2";
  // the error tolerance (Da) of the graph alignment
  // (GraphAlignMng::topmg_align_error_tolerance_), set by the hidden
  // --align-error-tolerance; toppic's hidden --topmg-align-error-tolerance
  // sets the same key
  arguments["topmgAlignErrorTolerance"] = "0.1";
  arguments["cutoffSpectralType"] = "EVALUE";
  arguments["cutoffSpectralValue"] = "0.01";
  arguments["cutoffProteoformType"] = "EVALUE";
  arguments["cutoffProteoformValue"] = "0.01";
  arguments["cutoffProteinType"] = "EVALUE";
  arguments["cutoffProteinValue"] = "0.01";
  arguments["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments["allowProtType"] = "COMPLETE,PREFIX,SUFFIX,INTERNAL";
  arguments["numOfTopPrsms"] = "1";
  arguments["maxShiftMass"] = "500";
  arguments["executiveDir"] = ".";
  arguments["resourceDir"] = "";
  arguments["keepTempFiles"] = "false";
  arguments["groupSpectrumNumber"] = "1";
  // the number of results the ASF filtering keeps per spectrum. TopMG's own
  // number: toppic's filteringResultNumber belongs to toppic's stages. The
  // default is the 15 that the ASF-Diagonal filtering used to hardcode.
  arguments["asfFilteringResultNumber"] = "15";
  arguments["varModFileName"] = "";
  arguments["threadNumber"] = "1";
  arguments["useFeatureFile"] = "true";
  // --no-sql: no SQLite file is written or changed: no <base>_topmg.sqlite
  // (toppic_test experiments 164 and 170)
  arguments["noSql"] = "false";
  // the hidden --disable-zstd: <base>_topmg.sqlite is written uncompressed
  // (by default it is compressed at the end of the run, sql/zstd_db;
  // toppic_test experiment 172)
  arguments["disableZstd"] = "false";
  // topmg's database of the file being reported, <base>_topmg.sqlite
  // (TopMGProgress_multi_file); "" for the combined result of -c
  arguments["topmgDbFileName"] = "";
  arguments["proteoGraphGap"] = "40";
  arguments["useAsfDiag"] = "false";
  arguments["varPtmNumber"] = "5";
  arguments["varPtmNumInGap"] = "5";
  arguments["keepDecoyResults"] = "false";
  // the filters' mass list (hidden --filter-mass-list): merge, fragment or
  // single (SpPara::FilterMassListMode)
  arguments["filterMassList"] = "fragment";
  arguments["version"] = "";
  return arguments;
}

void TopmgArgument::outputArguments(
    std::ostream& output, const std::string& sep,
    std::map<std::string, std::string> arguments) {
  int gap = 47;
  output << "********************** Parameters ***********************"
         << std::endl;
  output << std::setw(gap) << std::left << "Protein database file:" << sep
         << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum file:" << sep
         << arguments["spectrumFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Fragmentation method:" << sep
         << arguments["activation"] << std::endl;
  output << std::setw(gap) << std::left << "Search type:" << sep
         << arguments["searchType"] << std::endl;

  if (arguments["fixedMod"] != "") {
    // add fixed PTM information
    if (arguments["fixedMod"] == "C57") {
      output << std::setw(gap) << std::left << "Fixed modifications BEGIN"
             << std::endl;
      output << std::setw(gap) << std::left << "Carbamidomethylation" << sep
             << ptmMassStr(PtmBase::getPtmPtr_C57()) << sep << "C" << std::endl;
      output << std::setw(gap) << std::left << "Fixed modifications END"
             << std::endl;
    } else if (arguments["fixedMod"] == "C58") {
      output << std::setw(gap) << std::left << "Fixed modifications BEGIN"
             << std::endl;
      output << std::setw(gap) << std::left << "Carboxymethylation" << sep
             << ptmMassStr(PtmBase::getPtmPtr_C58()) << sep << "C" << std::endl;
      output << std::setw(gap) << std::left << "Fixed modifications END"
             << std::endl;
    } else {
      output << std::setw(gap) << std::left
             << "Fixed modifications file name:" << sep << arguments["fixedMod"]
             << std::endl;
      output << std::setw(gap) << std::left << "Fixed modifications BEGIN"
             << std::endl;
      std::vector<std::vector<std::string>> mod_data =
          mod_util::readModTxtForTsv(arguments["fixedMod"]);
      for (size_t i = 0; i < mod_data.size(); i++) {
        output << std::setw(gap) << std::left << mod_data[i][0] << sep
               << mod_data[i][1] << sep << mod_data[i][2] << std::endl;
      }
    }
  }
  output << std::setw(gap) << std::left << "Allowed N-terminal forms:" << sep
         << arguments["allowProtMod"] << std::endl;

  output << std::setw(gap) << std::left
         << "Maximum number of variable PTMs:" << sep
         << arguments["varPtmNumber"] << std::endl;
  output << std::setw(gap) << std::left
         << "Variable modification file name:" << sep
         << arguments["varModFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Variable modifications BEGIN"
         << std::endl;
  std::vector<std::vector<std::string>> mod_data =
      mod_util::readModTxtForTsv(arguments["varModFileName"]);
  for (size_t i = 0; i < mod_data.size(); i++) {
    output << std::setw(gap) << std::left << mod_data[i][0] << sep
           << mod_data[i][1] << sep << mod_data[i][2] << std::endl;
  }
  output << std::setw(gap) << std::left << "Variable modifications END"
         << std::endl;

  output << std::setw(gap) << std::left
         << "Maximum number of unexpected modifications:" << sep
         << arguments["topmgShiftNumber"] << std::endl;
  output << std::setw(gap) << std::left
         << "Maximum mass shift of modifications:" << sep
         << arguments["maxShiftMass"] << " Da" << std::endl;
  output << std::setw(gap) << std::left << "Gap in proteoform graph:" << sep
         << arguments["proteoGraphGap"] << std::endl;
  output << std::setw(gap) << std::left
         << "Maximum number of variable PTMs in a graph gap:" << sep
         << arguments["varPtmNumInGap"] << std::endl;

  // whole proteins only is the proteoform type COMPLETE alone
  if (arguments["allowProtType"] == "COMPLETE") {
    output << std::setw(gap) << std::left
           << "Report only proteoforms from whole proteins:" << sep << "True"
           << std::endl;
  } else {
    output << std::setw(gap) << std::left
           << "Report only proteoforms from whole proteins:" << sep << "False"
           << std::endl;
  }

  output << std::setw(gap) << std::left << "Spectrum-level cutoff type:" << sep
         << arguments["cutoffSpectralType"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum-level cutoff value:" << sep
         << arguments["cutoffSpectralValue"] << std::endl;
  output << std::setw(gap) << std::left
         << "Proteoform-level cutoff type:" << sep
         << arguments["cutoffProteoformType"] << std::endl;
  output << std::setw(gap) << std::left
         << "Proteoform-level cutoff value:" << sep
         << arguments["cutoffProteoformValue"] << std::endl;
  output << std::setw(gap) << std::left << "Protein-level cutoff type:" << sep
         << arguments["cutoffProteinType"] << std::endl;
  output << std::setw(gap) << std::left << "Protein-level cutoff value:" << sep
         << arguments["cutoffProteinValue"] << std::endl;

  output << std::setw(gap) << std::left
         << "Error tolerance for precursor masses:" << sep
         << arguments["precErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(gap) << std::left
         << "Error tolerance for fragment masses:" << sep
         << arguments["fragErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(gap) << std::left << "Confidence score cutoff:" << sep
         << arguments["envCnnCutoff"] << std::endl;
  output << std::setw(gap) << std::left
         << "Error tolerance for identifying PrSM clusters:" << sep
         << arguments["proteoformErrorTolerance"] << " Da" << std::endl;
  output << std::setw(gap) << std::left
         << "Error tolerance for graph alignment:" << sep
         << arguments["topmgAlignErrorTolerance"] << " Da" << std::endl;
  // printed only when it is not the default, so that the default parameter
  // block is unchanged
  if (arguments["filterMassList"] != "" &&
      arguments["filterMassList"] != "fragment") {
    output << std::setw(gap) << std::left << "Filter mass list:" << sep
           << arguments["filterMassList"] << std::endl;
  }

  if (arguments["useFeatureFile"] == "true") {
    output << std::setw(gap) << std::left << "Use TopFD feature file:" << sep
           << "True" << std::endl;
  } else {
    output << std::setw(gap) << std::left << "Use TopFD feature file:" << sep
           << "False" << std::endl;
  }
  // printed only when set, so that the default output files do not change
  if (arguments["noSql"] == "true") {
    output << std::setw(gap) << std::left << "SQLite database output:" << sep
           << "No" << std::endl;
  }

  output << std::setw(gap) << std::left << "Thread number:" << sep
         << arguments["threadNumber"] << std::endl;
  output << std::setw(gap) << std::left << "Executable file directory:" << sep
         << arguments["executiveDir"] << std::endl;
  output << std::setw(gap) << std::left << "Start time:" << sep
         << arguments["startTime"] << std::endl;
  if (arguments["endTime"] != "") {
    output << std::setw(gap) << std::left << "End time:" << sep
           << arguments["endTime"] << std::endl;
  }
  output << std::setw(gap) << std::left << "Version:" << sep
         << arguments["version"] << std::endl;
  output << "********************** Parameters ***********************"
         << std::endl;
}

std::string TopmgArgument::outputTsvArguments(
    std::map<std::string, std::string> arguments) {
  std::stringstream output;
  outputArguments(output, "\t", arguments);
  return output.str();
}

void TopmgArgument::showUsage(
    boost::program_options::options_description& desc) {
  std::cout << "Usage: topmg [options] -i mod-file-name database-file-name "
               "spectrum-file-names"
            << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool TopmgArgument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string activation = "";
  std::string fixed_mod = "";
  std::string allow_mod = "";
  std::string shift_num = "";
  std::string frag_error_tole = "";
  std::string prec_error_tole = "";
  std::string env_cnn_cutoff = "";
  std::string form_error_tole = "";
  std::string align_error_tole = "";
  std::string max_shift_mass = "";
  std::string cutoff_spectral_type = "";
  std::string cutoff_spectral_value = "";
  std::string cutoff_proteoform_type = "";
  std::string cutoff_proteoform_value = "";
  std::string cutoff_protein_type = "";
  std::string cutoff_protein_value = "";
  std::string use_table = "";
  std::string use_asf_diag = "";
  std::string filtering_result_num = "";
  std::string var_mod_file_name = "";
  std::string proteo_graph_gap = "";
  std::string thread_number = "";
  std::string var_ptm_num = "";
  std::string var_ptm_in_gap = "";
  std::string combined_output_name = "";
  std::string filter_mass_list = "";

  // Define and parse the progra/m options
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()("help,h", "Print the help message.")(
        "activation,a", po::value<std::string>(&activation),
        "<CID|ETD|ETCID|ETHCD|HCD|MPD|UVPD|FILE>. Fragmentation method of "
        "MS/MS spectra. When FILE is used, fragmentation methods of spectra "
        "are given in the input spectral data file. ETCID (ETciD) and ETHCD "
        "(EThcD) are searched with the ion types of ETD (c and z-dot). Default "
        "value: FILE.")(
        "fixed-mod,f", po::value<std::string>(&fixed_mod),
        "<C57|C58|a fixed modification file>. Fixed modifications. Three "
        "available options: C57, C58, or the name of a text file containing "
        "the information of fixed modifications. When C57 is selected, "
        "carbamidomethylation on cysteine is the only fixed modification. When "
        "C58 is selected, carboxymethylation on cysteine is the only fixed "
        "modification.")(
        "n-terminal-form,n", po::value<std::string>(&allow_mod),
        "<a list of allowed N-terminal forms>. N-terminal forms of proteins. "
        "Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, "
        "and M_ACETYLATION. NONE stands for no modifications, NME for "
        "N-terminal methionine excision, NME_ACETYLATION for N-terminal "
        "acetylation after the initiator methionine is removed, and "
        "M_ACETYLATION for N-terminal methionine acetylation. When multiple "
        "forms are allowed, they are separated by commas. Default value: "
        "NONE,NME,NME_ACETYLATION,M_ACETYLATION.")(
        "decoy,d",
        "Use a shuffled decoy protein database to estimate false discovery "
        "rates.")("frag-error-tolerance,e",
                  po::value<std::string>(&frag_error_tole),
                  "<a positive integer>. Error tolerance for fragment masses "
                  "in PPM. Default value: 10.")(
        "prec-error-tolerance,E", po::value<std::string>(&prec_error_tole),
        "<a positive integer>. Error tolerance for precursor masses in PPM. "
        "Default value: 10.")(
        "filter-by-score,F", po::value<std::string>(&env_cnn_cutoff),
        "<a number between 0 and 1>. Remove the deconvoluted masses whose "
        "confidence score (e.g., EnvCNN score) is below the cutoff before "
        "spectral filtering, database search and E-value computation. "
        "Default value: 0.4.")(
        "proteoform-error-tolerance,p",
        po::value<std::string>(&form_error_tole),
        "<a positive number>. Error tolerance for identifying PrSM clusters "
        "(in Dalton). Default value: 1.2 Dalton.")(
        "max-shift,M", po::value<std::string>(&max_shift_mass),
        "<a positive number>. Maximum absolute value of an unexpected mass "
        "shift (in Dalton). Default value: 500.")(
        "spectrum-cutoff-type,t", po::value<std::string>(&cutoff_spectral_type),
        "<EVALUE|FDR>. Spectrum-level cutoff type for filtering identified "
        "proteoform spectrum-matches. Default value: EVALUE.")(
        "spectrum-cutoff-value,v",
        po::value<std::string>(&cutoff_spectral_value),
        "<a positive number>. Spectrum-level cutoff value for filtering "
        "identified proteoform spectrum-matches. Default value: 0.01.")(
        "proteoform-cutoff-type,T",
        po::value<std::string>(&cutoff_proteoform_type),
        "<EVALUE|FDR>. Proteoform-level cutoff type for filtering identified "
        "proteoform spectrum-matches. Default value: EVALUE.")(
        "proteoform-cutoff-value,V",
        po::value<std::string>(&cutoff_proteoform_value),
        "<a positive number>. Proteoform-level cutoff value for filtering "
        "identified proteoform spectrum-matches. Default value: 0.01.")(
        "protein-cutoff-type,y", po::value<std::string>(&cutoff_protein_type),
        "<EVALUE|FDR>. Protein-level cutoff type for filtering identified "
        "proteoform spectrum-matches. With FDR, a protein is represented by "
        "its best proteoform (the lowest E-value). Default value: EVALUE.")(
        "protein-cutoff-value,Y", po::value<std::string>(&cutoff_protein_value),
        "<a positive number>. Protein-level cutoff value for filtering "
        "identified proteoform spectrum-matches. Default value: 0.01.")(
        "mod-file-name,i", po::value<std::string>(&var_mod_file_name),
        "<a modification file>. Specify a text file containing the information "
        "of common PTMs for constructing proteoform graphs. Required.")(
        "thread-number,u", po::value<std::string>(&thread_number),
        "<a positive integer>. Number of threads used in the computation. "
        "Default value: 1.")(
        "no-topfd-feature,x",
        "No TopFD feature file for proteoform identification.")(
        "proteo-graph-gap,j", po::value<std::string>(&proteo_graph_gap),
        "<a positive integer>. Gap in constructing proteoform graph: the "
        "maximum number of residues between two vertices whose mass "
        "difference is matched. Default value: 40.")("var-ptm-in-gap,G",
                      po::value<std::string>(&var_ptm_in_gap),
                      "<a positive integer>. Maximum number of variable PTMs "
                      "in a proteoform graph gap, at most --var-ptm. Default "
                      "value: 5.")(
        "use-asf-diagonal,D",
        "Use the ASF-DIAGONAL method for protein sequence filtering.")(
        "var-ptm,P", po::value<std::string>(&var_ptm_num),
        "<a positive integer>. Maximum number of variable PTMs in a "
        "proteoform. Default value: 5.")(
        "num-shift,s", po::value<std::string>(&shift_num),
        "<0|1|2>. Maximum number of unexpected modifications in a "
        "proteoform spectrum-match. Default value: 1.")(
        "whole-protein-only,w",
        "Report only proteoforms from whole proteins, that is, allow the "
        "proteoform type COMPLETE alone.")(
        "combined-file-name,c", po::value<std::string>(&combined_output_name),
        "Specify a file name for the combined spectrum data file and analysis "
        "results.")("keep-temp-files,k", "Keep intermediate files.")(
        "keep-decoy-ids,K", "Keep decoy identifications.")(
        "no-sql,z",
        "Do not write an SQLite database. By default, topmg writes "
        "<name>_topmg.sqlite with the MS/MS envelope tables of the TopFD "
        "database <name>.sqlite and the PrSM, proteoform and protein tables. "
        "With this option, no SQLite file is written or changed.");

    po::options_description desc("Options");

    desc.add_options()("help,h", "")("activation,a",
                                     po::value<std::string>(&activation), "")(
        "fixed-mod,f", po::value<std::string>(&fixed_mod), "")(
        "n-terminal-form,n", po::value<std::string>(&allow_mod), "")(
        "decoy,d", "")("frag-error-tolerance,e",
                       po::value<std::string>(&frag_error_tole), "")(
        "prec-error-tolerance,E", po::value<std::string>(&prec_error_tole), "")(
        "filter-by-score,F", po::value<std::string>(&env_cnn_cutoff), "")(
        "proteoform-ppm-error,Q", "")("proteoform-error-tolerance,p",
                                      po::value<std::string>(&form_error_tole),
                                      "")(
        // hidden (not in display_desc): the graph alignment's error tolerance
        "align-error-tolerance", po::value<std::string>(&align_error_tole),
        "")(
        "max-shift,M", po::value<std::string>(&max_shift_mass), "")(
        "spectrum-cutoff-type,t", po::value<std::string>(&cutoff_spectral_type),
        "")("spectrum-cutoff-value,v",
            po::value<std::string>(&cutoff_spectral_value),
            "")("proteoform-cutoff-type,T",
                po::value<std::string>(&cutoff_proteoform_type),
                "")("proteoform-cutoff-value,V",
                    po::value<std::string>(&cutoff_proteoform_value),
                    "")("protein-cutoff-type,y",
                        po::value<std::string>(&cutoff_protein_type), "")(
        "protein-cutoff-value,Y", po::value<std::string>(&cutoff_protein_value),
        "")("filtering-result-number",
            po::value<std::string>(&filtering_result_num),
            "ASF-Diagonal filtering result number. Default value: 15.")(
        "full-binary-path,b",
                                                           "Full binary path.")(
        "mod-file-name,i", po::value<std::string>(&var_mod_file_name), "")(
        "thread-number,u", po::value<std::string>(&thread_number), "")(
        "no-topfd-feature,x", "")(
        "combined-file-name,c", po::value<std::string>(&combined_output_name),
        "")("proteo-graph-gap,j", po::value<std::string>(&proteo_graph_gap),
            "")("var-ptm-in-gap,G", po::value<std::string>(&var_ptm_in_gap),
                "")("use-asf-diagonal,D", "")(
        "filter-mass-list", po::value<std::string>(&filter_mass_list), "")(
        "var-ptm,P", po::value<std::string>(&var_ptm_num), "")(
        "num-shift,s", po::value<std::string>(&shift_num), "")(
        "whole-protein-only,w", "")("keep-temp-files,k", "")("keep-decoy-ids,K",
                                                             "")("no-sql,z", "")(
        "disable-zstd", "")(
        "database-file-name",
        po::value<std::string>(&database_file_name)->required(),
        "Database file name with its path.")(
        "spectrum-file-name",
        po::value<std::vector<std::string>>()->multitoken()->required(),
        "Spectrum file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv)
                    .options(desc)
                    .positional(positional_options)
                    .run(),
                vm);
      if (vm.count("help")) {
        showUsage(display_desc);
        return false;
      }
      po::notify(vm);
      // throws on error, so do after help in case there are any problems
    } catch (boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    } catch (boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }

    std::string argv_0(argv[0]);
    if (vm.count("full-binary-path")) {
      arguments_["executiveDir"] = argv[0];
    } else {
      arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
      if (file_util::checkSpace(arguments_["executiveDir"])) {
        LOG_ERROR("Current directory "
                  << arguments_["executiveDir"]
                  << " contains space and will cause errors in the program!");
        exit(EXIT_FAILURE);
      }
    }
    LOG_DEBUG("Executive Dir " << arguments_["executiveDir"]);

    arguments_["resourceDir"] =
        file_util::getResourceDir(arguments_["executiveDir"]);

    arguments_["oriDatabaseFileName"] = database_file_name;

    if (vm.count("spectrum-file-name")) {
      spec_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string>>();
    }

    if (vm.count("combined-file-name")) {
      arguments_["combinedOutputName"] = combined_output_name;
    }

    if (vm.count("activation")) {
      arguments_["activation"] = activation;
    }

    if (vm.count("decoy")) {
      arguments_["searchType"] = "TARGET+DECOY";
    }


    if (vm.count("fixed-mod")) {
      arguments_["fixedMod"] = fixed_mod;
    }

    if (vm.count("n-terminal-form")) {
      arguments_["allowProtMod"] = allow_mod;
    }

    if (vm.count("frag-error-tolerance")) {
      arguments_["fragErrorTolerance"] = frag_error_tole;
    }

    if (vm.count("prec-error-tolerance")) {
      arguments_["precErrorTolerance"] = prec_error_tole;
    }

    if (vm.count("filter-by-score")) {
      arguments_["envCnnCutoff"] = env_cnn_cutoff;
    }

    if (vm.count("proteoform-ppm-error")) {
      arguments_["proteoformPpmError"] = "true";
    }

    if (vm.count("proteoform-error-tolerance")) {
      arguments_["proteoformErrorTolerance"] = form_error_tole;
    }

    if (vm.count("align-error-tolerance")) {
      arguments_["topmgAlignErrorTolerance"] = align_error_tole;
    }

    if (vm.count("max-shift")) {
      arguments_["maxShiftMass"] = max_shift_mass;
    }

    if (vm.count("spectrum-cutoff-type")) {
      arguments_["cutoffSpectralType"] = cutoff_spectral_type;
    }

    if (vm.count("spectrum-cutoff-value")) {
      arguments_["cutoffSpectralValue"] = cutoff_spectral_value;
    }

    if (vm.count("proteoform-cutoff-type")) {
      arguments_["cutoffProteoformType"] = cutoff_proteoform_type;
    }

    if (vm.count("proteoform-cutoff-value")) {
      arguments_["cutoffProteoformValue"] = cutoff_proteoform_value;
    }

    if (vm.count("protein-cutoff-type")) {
      arguments_["cutoffProteinType"] = cutoff_protein_type;
    }

    if (vm.count("protein-cutoff-value")) {
      arguments_["cutoffProteinValue"] = cutoff_protein_value;
    }

    if (vm.count("filtering-result-number")) {
      arguments_["asfFilteringResultNumber"] = filtering_result_num;
    }

    if (vm.count("mod-file-name")) {
      arguments_["varModFileName"] = var_mod_file_name;
    }

    if (vm.count("proteo-graph-gap")) {
      arguments_["proteoGraphGap"] = proteo_graph_gap;
    }

    if (vm.count("thread-number")) {
      arguments_["threadNumber"] = thread_number;
    }

    if (vm.count("no-topfd-feature")) {
      arguments_["useFeatureFile"] = "false";
    }

    if (vm.count("use-asf-diagonal")) {
      arguments_["useAsfDiag"] = "true";
    }

    // hidden: the entries of the filters' mass list
    // (prm_ms_util::getIntMassErrorList) of the one-shift filter (and the
    // multi-shift filter with --use-asf-diagonal): "merge", one per integer
    // mass (main up to 7b5112c2); "fragment", one per integer mass and base
    // type (the default); "single", one per peak (toppic_test experiment
    // 142)
    if (vm.count("filter-mass-list")) {
      arguments_["filterMassList"] = filter_mass_list;
    }

    if (vm.count("var-ptm")) {
      arguments_["varPtmNumber"] = var_ptm_num;
    }

    if (vm.count("num-shift")) {
      arguments_["topmgShiftNumber"] = shift_num;
    }

    if (vm.count("var-ptm-in-gap")) {
      arguments_["varPtmNumInGap"] = var_ptm_in_gap;
    }

    // --whole-protein-only is the proteoform type COMPLETE alone, from which
    // TopMG_identify reads it; there is no separate wholeProteinOnly argument
    if (vm.count("whole-protein-only")) {
      arguments_["allowProtType"] = "COMPLETE";
    }

    if (vm.count("keep-temp-files")) {
      arguments_["keepTempFiles"] = "true";
    }

    if (vm.count("keep-decoy-ids")) {
      arguments_["keepDecoyResults"] = "true";
    }

    if (vm.count("no-sql")) {
      arguments_["noSql"] = "true";
    }

    if (vm.count("disable-zstd")) {
      arguments_["disableZstd"] = "true";
    }
  } catch (std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line " << e.what()
              << ", application will now exit." << std::endl;
    return false;
  }

  return validateArguments();
}

bool TopmgArgument::validateArguments() {
  if (!std::filesystem::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource directory "
              << arguments_["resourceDir"] << " does not exist!\n"
              << "Please check if file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  if (!std::filesystem::exists(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file "
              << arguments_["oriDatabaseFileName"] << " does not exist!\n"
              << "Please check if file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  if (!str_util::endsWith(arguments_["oriDatabaseFileName"], ".fasta") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".fa") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".FASTA") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".FA")) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"]
                               << " is not a fasta file!");
    return false;
  }

  if (arguments_["oriDatabaseFileName"].length() > 200) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"]
                               << " path is too long!");
    return false;
  }

  // A database whose protein ids (the header text up to the first
  // whitespace) repeat cannot be indexed by htslib: search its
  // _dedup.fasta copy, which keeps the first record of every id, instead.
  arguments_["oriDatabaseFileName"] =
      fasta_util::removeDuplicatedIds(arguments_["oriDatabaseFileName"]);

  if (arguments_["searchType"] == "TARGET+DECOY") {
    arguments_["databaseFileName"] =
        arguments_["oriDatabaseFileName"] + "_target_decoy";
  } else {
    arguments_["databaseFileName"] =
        arguments_["oriDatabaseFileName"] + "_target";
  }

  for (size_t k = 0; k < spec_file_list_.size(); k++) {
    if (!std::filesystem::exists(spec_file_list_[k])) {
      LOG_ERROR("Spectrum file "
                << spec_file_list_[k]
                << " does not exist!\nPlease check if file directory contains "
                   "unproper characters such as spaces/quotation makrks");
      return false;
    }

    if (!str_util::endsWith(spec_file_list_[k], ".msalign")) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k]
                                 << " is not a msalign file!");
      return false;
    }

    if (spec_file_list_[k].length() > 200) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k] << " path is too long!");
      return false;
    }

    if (str_util::endsWith(spec_file_list_[k], "_ms1.msalign")) {
      std::cerr << "Warning: Please make sure " << spec_file_list_[k]
                << " is the ms2 spectral file." << std::endl;
    }
  }

  if (arguments_["varModFileName"] == "") {
    LOG_ERROR(
        "Modification file is not set. Please use the -i argument to specify "
        "the modification file!");
    return false;
  }

  if (!SpPara::isFilterMassListModeName(arguments_["filterMassList"])) {
    LOG_ERROR("--filter-mass-list " << arguments_["filterMassList"]
                                    << " error! The value should be merge, "
                                       "fragment or single!");
    return false;
  }

  if (!std::filesystem::exists(arguments_["varModFileName"])) {
    LOG_ERROR("Modification file "
              << arguments_["varModFileName"]
              << " does not exist!\nPlease check if file directory contains "
                 "unproper characters such as spaces/quotation makrks");
    return false;
  }

  std::string activation = arguments_["activation"];
  if (activation != "CID" && activation != "HCD" && activation != "MPD" &&
      activation != "ETD" && activation != "ETCID" && activation != "ETHCD" &&
      activation != "FILE" && activation != "UVPD") {
    LOG_ERROR("Activation type "
              << activation
              << " error! The value should be "
                 "CID|ETD|ETCID|ETHCD|HCD|MPD|UVPD|FILE!");
    return false;
  }

  std::string search_type = arguments_["searchType"];
  if (search_type != "TARGET" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Search type "
              << search_type
              << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
  }
  std::string shift_number = arguments_["topmgShiftNumber"];
  if (shift_number != "0" && shift_number != "1" && shift_number != "2") {
    LOG_ERROR("Shift number " << shift_number
                              << " error! The value should be 0|1|2!");
    return false;
  }

  std::string allow_mod = arguments_["allowProtMod"];
  std::vector<std::string> strs = str_util::split(allow_mod, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    if (strs[i] != "NONE" && strs[i] != "M_ACETYLATION" && strs[i] != "NME" &&
        strs[i] != "NME_ACETYLATION") {
      LOG_ERROR(
          "N-Terminal Variable PTM can only be NONE, M_ACETYLATION, NME or "
          "NME_ACETYLATION.");
      return false;
    }
  }

  std::string cutoff_spectral_type = arguments_["cutoffSpectralType"];
  if (cutoff_spectral_type != "EVALUE" && cutoff_spectral_type != "FDR") {
    LOG_ERROR("Spectrum-level cutoff type "
              << cutoff_spectral_type
              << " error! The value should be EVALUE|FDR");
    return false;
  }

  std::string cutoff_proteoform_type = arguments_["cutoffProteoformType"];
  if (cutoff_proteoform_type != "EVALUE" && cutoff_proteoform_type != "FDR") {
    LOG_ERROR("Proteoform-level cutoff type "
              << cutoff_proteoform_type
              << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_spectral_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Spectrum-level cutoff type "
              << cutoff_spectral_type
              << " error! FDR cutoff cannot be used when no decoy database is "
                 "used! Please add argument '-d' in the command.");
    return false;
  }

  if (cutoff_proteoform_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Proteoform-level cutoff type "
              << cutoff_proteoform_type
              << " error! FDR cutoff cannot be used when no decoy database is "
                 "used! Please add argument '-d' in the command.");
    return false;
  }

  std::string cutoff_protein_type = arguments_["cutoffProteinType"];
  if (cutoff_protein_type != "EVALUE" && cutoff_protein_type != "FDR") {
    LOG_ERROR("Protein-level cutoff type "
              << cutoff_protein_type
              << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_protein_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Protein-level cutoff type "
              << cutoff_protein_type
              << " error! FDR cutoff cannot be used when no decoy database is "
                 "used! Please add argument '-d' in the command.");
    return false;
  }

  std::string align_error_tole = arguments_["topmgAlignErrorTolerance"];
  try {
    if (std::stod(align_error_tole) <= 0) {
      LOG_ERROR("Graph alignment error tolerance "
                << align_error_tole << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Graph alignment error tolerance " << align_error_tole
                                                 << " should be a number.");
    return false;
  }

  std::string max_shift_mass = arguments_["maxShiftMass"];
  try {
    double mass = std::stod(max_shift_mass.c_str());
    if (mass <= 0.0) {
      LOG_ERROR("Maximum shift mass "
                << max_shift_mass << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Maximum shift mass " << max_shift_mass
                                    << " should be a number.");
    return false;
  }

  std::string cutoff_spectral_value = arguments_["cutoffSpectralValue"];
  try {
    double th = std::stod(cutoff_spectral_value);
    if (th < 0) {
      LOG_ERROR("Spectrum-level cutoff value "
                << cutoff_spectral_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Spectrum-level cutoff value " << cutoff_spectral_value
                                             << " should be a number.");
    return false;
  }

  std::string cutoff_proteoform_value = arguments_["cutoffProteoformValue"];
  try {
    double th = std::stod(cutoff_proteoform_value);
    if (th < 0) {
      LOG_ERROR("Proteoform-level cutoff value "
                << cutoff_proteoform_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Proteoform-level cutoff value " << cutoff_proteoform_value
                                               << " should be a number.");
    return false;
  }

  std::string cutoff_protein_value = arguments_["cutoffProteinValue"];
  try {
    double th = std::stod(cutoff_protein_value);
    if (th < 0) {
      LOG_ERROR("Protein-level cutoff value "
                << cutoff_protein_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Protein-level cutoff value " << cutoff_protein_value
                                            << " should be a number.");
    return false;
  }

  std::string thread_number = arguments_["threadNumber"];
  try {
    int num = std::stoi(thread_number.c_str());
    bool valid = mem_check::checkThreadNum(num, "topmg");
    if (!valid) {
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }

  for (const auto& [key, name] :
       {std::pair<const char*, const char*>{"precErrorTolerance",
                                            "Precursor error tolerance"},
        std::pair<const char*, const char*>{"fragErrorTolerance",
                                            "Fragment error tolerance"}}) {
    std::string tole_value = arguments_[key];
    try {
      double tole = std::stoi(tole_value);
      if (tole <= 0) {
        LOG_ERROR(name << ": " << tole_value
                       << " error! The value should be positive.");
        return false;
      }
    } catch (const std::exception& ex) {
      LOG_ERROR(name << ": " << tole_value << " should be a number.");
      return false;
    }
  }

  std::string env_cnn_cutoff_value = arguments_["envCnnCutoff"];
  try {
    double cutoff = std::stod(env_cnn_cutoff_value);
    if (cutoff < 0 || cutoff > 1) {
      LOG_ERROR("Confidence score cutoff: "
                << env_cnn_cutoff_value
                << " error! The value should be between 0 and 1.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("Confidence score cutoff: " << env_cnn_cutoff_value
                                      << " should be a number.");
    return false;
  }

  std::string form_error_tole_value = arguments_["proteoformErrorTolerance"];
  try {
    double tole = std::stod(form_error_tole_value);
    if (tole <= 0) {
      LOG_ERROR("PrSM clustering error tolerance: "
                << form_error_tole_value
                << " error! The value should be positive.");
      return false;
    }
  } catch (const std::exception& ex) {
    LOG_ERROR("PrSM clustering error tolerance: " << form_error_tole_value
                                                  << " should be a number.");
    return false;
  }

  // -P, -j and -G are read with std::stoi by TopMG_identify
  for (const auto& [key, name] :
       {std::pair<const char*, const char*>{"varPtmNumber",
                                            "Maximum variable PTM number"},
        std::pair<const char*, const char*>{"proteoGraphGap",
                                            "Proteoform graph gap"},
        std::pair<const char*, const char*>{"varPtmNumInGap",
                                            "Maximum variable PTM number in "
                                            "a gap"}}) {
    std::string value = arguments_[key];
    try {
      if (std::stoi(value) <= 0) {
        LOG_ERROR(name << " " << value
                       << " error! The value should be a positive integer.");
        return false;
      }
    } catch (const std::exception& ex) {
      LOG_ERROR(name << " " << value << " should be a positive integer.");
      return false;
    }
  }
  return true;
}
} /* namespace toppic */
