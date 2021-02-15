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

#ifndef BOOST_SYSTEM_NO_DEPRECATED
#define BOOST_SYSTEM_NO_DEPRECATED 1
#endif

#include <iomanip>
#include <vector>
#include <string>
#include <set>
#include <map>

#include "boost/thread/thread.hpp"

#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/util/str_util.hpp"

#include "console/topmg_argument.hpp"

namespace toppic {

Argument::Argument() {
  initArguments();
}

void Argument::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "400000000";
  arguments_["maxFragmentLength"] = "1000";
  arguments_["spectrumFileName"] = "";
  arguments_["combinedOutputName"] = "";
  arguments_["activation"] = "FILE";
  arguments_["searchType"] = "TARGET";
  arguments_["fixedMod"] = "";
  arguments_["ptmNumber"] = "0";
  arguments_["massErrorTolerance"] = "15";
  arguments_["proteoformErrorTolerance"] = "1.2";
  arguments_["cutoffSpectralType"] = "EVALUE";
  arguments_["cutoffSpectralValue"] = "0.01";
  arguments_["cutoffProteoformType"] = "EVALUE";
  arguments_["cutoffProteoformValue"] = "0.01";
  arguments_["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments_["numOfTopPrsms"] = "1";
  arguments_["maxPtmMass"] = "500";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["keepTempFiles"] = "false";
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["filteringResultNumber"] = "20";
  arguments_["varModFileName"] = "";
  arguments_["threadNumber"] = "1";
  arguments_["useFeatureFile"] = "true";
  arguments_["proteoGraphGap"] = "40";
  arguments_["useAsfDiag"] = "false";
  arguments_["varPtmNumber"] = "5";
  arguments_["varPtmNumInGap"] = "5";
  arguments_["geneHTMLFolder"] = "true";
  arguments_["wholeProteinOnly"] = "false";
}

void Argument::outputArguments(std::ostream &output, std::map<std::string, std::string> arguments) {
  output << "********************** Parameters **********************" << std::endl;
  output << std::setw(50) << std::left << "Protein database file: " << "\t" << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(50) << std::left << "Spectrum file: " << "\t" << arguments["spectrumFileName"] << std::endl;
  output << std::setw(50) << std::left << "Fragmentation method: " << "\t" << arguments["activation"] << std::endl;
  output << std::setw(50) << std::left << "Search type: " << "\t" << arguments["searchType"] << std::endl;

  if (arguments["fixedMod"] == "") {
    output << std::setw(50) << std::left << "Fixed modifications: " << "\t" << "None" << std::endl;
  } 
  else if (arguments["fixedMod"] == "C57") {
    output << std::setw(50) << std::left << "Fixed modifications: " << "\t" << "C57:carbamidomethylation on cysteine" << std::endl;
  }
  else if (arguments["fixedMod"] == "C58") {
    output << std::setw(50) << std::left << "Fixed modifications: " << "\t" << "C58:carboxymethylation on cysteine" << std::endl;
  }
  else {
    output << std::setw(50) << std::left << "Fixed modifications:," << arguments["fixedMod"] << std::endl;
  }

  if (arguments["useFeatureFile"] == "true") {
    output << std::setw(50) << std::left << "Use TopFD feature file: " << "\t" << "True" << std::endl;
  }
  else {
    output << std::setw(50) << std::left << "Use TopFD feature file: " << "\t" << "False" << std::endl;
  }

  output << std::setw(50) << std::left << "Error tolerance for matching masses: " << "\t" << arguments["massErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(50) << std::left << "Error tolerance for identifying PrSM clusters: " << "\t" << arguments["proteoformErrorTolerance"]  
      << " Da" << std::endl;
  output << std::setw(50) << std::left << "Spectrum-level cutoff type: " << "\t" << arguments["cutoffSpectralType"] << std::endl;
  output << std::setw(50) << std::left << "Spectrum-level cutoff value: " << "\t" << arguments["cutoffSpectralValue"] << std::endl;
  output << std::setw(50) << std::left << "Proteoform-level cutoff type: " << "\t" << arguments["cutoffProteoformType"] << std::endl;
  output << std::setw(50) << std::left << "Proteoform-level cutoff value: " << "\t" << arguments["cutoffProteoformValue"] << std::endl;
  output << std::setw(50) << std::left << "Allowed N-terminal forms: " << "\t" << arguments["allowProtMod"] << std::endl;
  output << std::setw(50) << std::left << "Maximum mass shift of modifications: " << "\t" << arguments["maxPtmMass"] << " Da" << std::endl;
  output << std::setw(50) << std::left << "Thread number: " << "\t" << arguments["threadNumber"] << std::endl;
  output << std::setw(50) << std::left << "Modification file name: " << "\t" << arguments["varModFileName"] << std::endl;
  output << std::setw(50) << std::left << "Gap in proteoform graph: " << "\t" << arguments["proteoGraphGap"] << std::endl;
  output << std::setw(50) << std::left << "Maximum number of variable PTMs: " << "\t" << arguments["varPtmNumber"] << std::endl;
  output << std::setw(50) << std::left << "Maximum number of variable PTMs in a graph gap: " << "\t" << arguments["varPtmNumInGap"] << std::endl;
  output << std::setw(50) << std::left << "Maximum number of unexpected modifications: " << "\t" << arguments["ptmNumber"] << std::endl;
  output << std::setw(50) << std::left << "Report only proteoforms from whole proteins: " << "\t" << arguments["wholeProteinOnly"] << std::endl;
  output << std::setw(50) << std::left << "Executable file directory: " << "\t" << arguments["executiveDir"] << std::endl;
  output << std::setw(50) << std::left << "Start time: " << "\t" << arguments["startTime"] << std::endl;
  if (arguments["endTime"] != "") {
    output << std::setw(50) << std::left << "End time: " << "\t" << arguments["endTime"] << std::endl;
  }
  output << "********************** Parameters **********************" << std::endl;
}

std::string Argument::outputTsvArguments(std::map<std::string, std::string> arguments) {
  std::stringstream output;
  outputArguments(output, arguments);
  return output.str();
}

void Argument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: topmg [options] database-file-name spectrum-file-name" << std::endl;
  std::cout << desc << std::endl;
}

bool Argument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string activation = "";
  std::string fixed_mod = "";
  std::string allow_mod = "";
  std::string ptm_num = "";
  std::string mass_error_tole = "";
  std::string form_error_tole = "";
  std::string max_ptm_mass = "";
  std::string cutoff_spectral_type = "";
  std::string cutoff_spectral_value = "";
  std::string cutoff_proteoform_type = "";
  std::string cutoff_proteoform_value = "";
  std::string use_table = "";
  std::string use_asf_diag = "";
  std::string filtering_result_num = "";
  std::string var_mod_file_name = "";
  std::string proteo_graph_gap = "";
  std::string thread_number = "";
  std::string var_ptm_num = "";
  std::string var_ptm_in_gap = "";
  std::string combined_output_name = "";

  // Define and parse the program options
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()
        ("help,h", "Print the help message.")
        ("activation,a", po::value<std::string>(&activation),
         "<CID|HCD|ETD|UVPD|FILE>. Fragmentation method of tandem mass spectra. When FILE is used, fragmentation methods of spectra are given in the input spectral data file. Default value: FILE.")
        ("fixed-mod,f", po::value<std::string> (&fixed_mod),
         "<C57|C58|a fixed modification file>. Fixed modifications. Three available options: C57, C58, or the name of a text file containing the information of fixed modifications. When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.")
        ("n-terminal-form,n", po::value<std::string> (&allow_mod),
         "<a list of allowed N-terminal forms>. N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,NME,NME_ACETYLATION,M_ACETYLATION.")
        ("decoy,d", "Use a decoy protein database to estimate false discovery rates.")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "<a positive integer>. Error tolerance for precursor and fragment masses in PPM. Default value: 15.")
        ("proteoform-error-tolerance,p", po::value<std::string> (&form_error_tole), "<a positive number>. Error tolerance for identifying PrSM clusters (in Dalton). Default value: 1.2 Dalton.")
        ("max-shift,M", po::value<std::string> (&max_ptm_mass), "<a positive number>. Maximum absolute value of the mass shift (in Dalton). Default value: 500.")
        ("spectrum-cutoff-type,t", po::value<std::string> (&cutoff_spectral_type), "<EVALUE|FDR>. Spectrum-level cutoff type for filtering identified proteoform spectrum-matches. Default value: EVALUE.")
        ("spectrum-cutoff-value,v", po::value<std::string> (&cutoff_spectral_value), "<a positive number>. Spectrum-level cutoff value for filtering identified proteoform spectrum-matches. Default value: 0.01.")
        ("proteoform-cutoff-type,T", po::value<std::string> (&cutoff_proteoform_type), "<EVALUE|FDR>. Proteoform-level cutoff type for filtering identified proteoform spectrum-matches. Default value: EVALUE.")
        ("proteoform-cutoff-value,V", po::value<std::string> (&cutoff_proteoform_value), "<a positive number>. Proteoform-level cutoff value for filtering identified proteoform spectrum-matches. Default value: 0.01.")
        ("mod-file-name,i", po::value<std::string>(&var_mod_file_name), "<a modification file>. Specify a text file containing the information of common PTMs for constructing proteoform graphs.")
        ("thread-number,u", po::value<std::string> (&thread_number), "<a positive integer>. Number of threads used in the computation. Default value: 1.")
        ("no-topfd-feature,x", "No TopFD feature file for proteoform identification.")
        ("proteo-graph-gap,j", po::value<std::string> (&proteo_graph_gap), "<a positive number>. Gap in constructing proteoform graph. Default value: 40.")
        ("var-ptm-in-gap,G", po::value<std::string>(&var_ptm_in_gap) , "<a positive number>. Maximum number of variable PTMs in a proteform graph gap. Default value: 5.")
        ("use-asf-diagonal,D", "Use the ASF-DIAGONAL method for protein sequence filtering.")
        ("var-ptm,P", po::value<std::string>(&var_ptm_num) , "<a positive number>. Maximum number of variable PTMs. Default value: 5.")
        ("num-shift,s", po::value<std::string> (&ptm_num), "<0|1|2>. Maximum number of unexpected modifications in a proteoform spectrum-match. Default value: 0.")
        ("whole-protein-only,w", "Report only proteoforms from whole proteins.")
        ("combined-file-name,c", po::value<std::string>(&combined_output_name) , "Specify a file name for the combined spectrum data file and analysis results.")
        ("keep-temp-files,k", "Keep temporary files.")
        ("skip-html-folder,g", "Skip generating an html folder containing TopView and spectrum data for visualization.");
  
    po::options_description desc("Options");

    desc.add_options()
        ("help,h", "")
        ("activation,a", po::value<std::string>(&activation), "")
        ("fixed-mod,f", po::value<std::string> (&fixed_mod), "")
        ("n-terminal-form,n", po::value<std::string> (&allow_mod), "")
        ("decoy,d", "")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "")
        ("proteoform-error-tolerance,p", po::value<std::string> (&form_error_tole), "")
        ("max-shift,M", po::value<std::string> (&max_ptm_mass), "")
        ("spectrum-cutoff-type,t", po::value<std::string> (&cutoff_spectral_type), "")
        ("spectrum-cutoff-value,v", po::value<std::string> (&cutoff_spectral_value), "")
        ("proteoform-cutoff-type,T", po::value<std::string> (&cutoff_proteoform_type), "")
        ("proteoform-cutoff-value,V", po::value<std::string> (&cutoff_proteoform_value), "")
        ("filtering-result-number", po::value<std::string>(&filtering_result_num), "Filtering result number. Default value: 20.")
        ("full-binary-path,b", "Full binary path.")
        ("mod-file-name,i", po::value<std::string>(&var_mod_file_name), "")
        ("thread-number,u", po::value<std::string> (&thread_number), "")
        ("no-topfd-feature,x", "")
        ("combined-file-name,c", po::value<std::string>(&combined_output_name) , "")
        ("proteo-graph-gap,j", po::value<std::string> (&proteo_graph_gap), "")
        ("var-ptm-in-gap,G", po::value<std::string>(&var_ptm_in_gap) , "")
        ("use-asf-diagonal,D", "")
        ("var-ptm,P", po::value<std::string>(&var_ptm_num) , "")
        ("num-shift,s", po::value<std::string> (&ptm_num), "")
        ("whole-protein-only,w", "")
        ("keep-temp-files,k", "")
        ("skip-html-folder,g", "")
        ("database-file-name", po::value<std::string>(&database_file_name)->required(), "Database file name with its path.")
        ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(), "Spectrum file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), vm);
      if (vm.count("help")) {
        showUsage(display_desc);
        return false;
      }
      po::notify(vm);
      // throws on error, so do after help in case there are any problems
    }
    catch(boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    catch(boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }

    std::string argv_0(argv[0]);
    if (vm.count("full-binary-path")) {
      arguments_["executiveDir"] = argv[0];
    } else {
      arguments_["executiveDir"] = file_util::getExecutiveDir(argv_0);
    }
    LOG_DEBUG("Executive Dir " << arguments_["executiveDir"]);

    arguments_["resourceDir"] = file_util::getResourceDir(arguments_["executiveDir"]);

    arguments_["oriDatabaseFileName"] = database_file_name;

    if (vm.count("spectrum-file-name")) {
      spec_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string> >(); 
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

    if (arguments_["searchType"] == "TARGET+DECOY") {
      arguments_["databaseFileName"] = arguments_["oriDatabaseFileName"] + "_target_decoy";
    } else {
      arguments_["databaseFileName"] = arguments_["oriDatabaseFileName"] + "_target";
    }

    if (vm.count("fixed-mod")) {
      arguments_["fixedMod"] = fixed_mod;
    }

    if (vm.count("n-terminal-form")) {
      arguments_["allowProtMod"] = allow_mod;
    }

    if (vm.count("mass-error-tolerance")) {
      arguments_["massErrorTolerance"] = mass_error_tole;
    }

    if (vm.count("proteoform-error-tolerance")) {
      arguments_["proteoformErrorTolerance"] = form_error_tole;
    }

    if (vm.count("max-shift")) {
      arguments_["maxPtmMass"] = max_ptm_mass;
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


    if (vm.count("filtering-result-number")) {
      arguments_["filteringResultNumber"] = filtering_result_num;
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

    if (vm.count("var-ptm")) {
      arguments_["varPtmNumber"] = var_ptm_num;
    }

    if (vm.count("num-shift")) {
      arguments_["ptmNumber"] = ptm_num;
    }

    if (vm.count("var-ptm-in-gap")) {
      arguments_["varPtmNumInGap"] = var_ptm_in_gap;
    }

    if (vm.count("whole-protein-only")) {
      arguments_["wholeProteinOnly"] = "true";
    }

    if (vm.count("keep-temp-files")) {
      arguments_["keepTempFiles"] = "true";
    }

    if (vm.count("skip-html-folder")) {
      arguments_["geneHTMLFolder"] = "false";
    }
  }
  catch(std::exception & e) {
    std::cerr << "Unhandled Exception in parsing command line " << e.what() << ", application will now exit." << std::endl;
    return false;
  }

  return validateArguments();
}

bool Argument::validateArguments() {
  if (!file_util::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource direcotry " << arguments_["resourceDir"] << " does not exist!");
    return false;
  }

  if (!file_util::exists(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file " << arguments_["databaseFileName"] << " does not exist!");
    return false;
  }

  if (!str_util::endsWith(arguments_["oriDatabaseFileName"], ".fasta") &&
      !str_util::endsWith(arguments_["oriDatabaseFileName"], ".fa")) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"] << " is not a fasta file!");
    return false;
  }

  if (arguments_["oriDatabaseFileName"].length() > 200) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"] << " path is too long!");
    return false;
  }

  for (size_t k = 0; k < spec_file_list_.size(); k++) {
    if (!file_util::exists(spec_file_list_[k])) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k] << " does not exist!");
      return false;
    }

    if (!str_util::endsWith(spec_file_list_[k], ".msalign")) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k] << " is not a msalign file!");
      return false;
    }

    if (spec_file_list_[k].length() > 200) {
      LOG_ERROR("Spectrum file " << spec_file_list_[k] << " path is too long!");
      return false;
    }

    if (str_util::endsWith(spec_file_list_[k], "_ms1.msalign")) {
      std::cerr << "Warning: Please make sure " << spec_file_list_[k] << " is the ms2 spectral file." << std::endl;
    }
  }

  if (!file_util::exists(arguments_["varModFileName"])) {
    LOG_ERROR("Modification file " << arguments_["varModFileName"] << " does not exist!");
    return false;
  }

  std::string activation = arguments_["activation"];
  if (activation != "CID" && activation != "HCD"
      && activation != "ETD" && activation != "FILE" && activation != "UVPD") {
    LOG_ERROR("Activation type " << activation << " error! The value should be CID|HCD|ETD|UVPD|FILE!");
    return false;
  }

  std::string search_type = arguments_["searchType"];
  if (search_type != "TARGET" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Search type " << search_type << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
  }
  std::string ptm_number = arguments_["ptmNumber"];
  if (ptm_number != "0" && ptm_number != "1" && ptm_number != "2") {
    LOG_ERROR("PTM number "<< ptm_number <<" error! The value should be 0|1|2!");
    return false;
  }

  std::string allow_mod = arguments_["allowProtMod"];
  std::vector<std::string> strs = str_util::split(allow_mod, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    if (strs[i] != "NONE" && strs[i] != "M_ACETYLATION" && strs[i] != "NME" && strs[i] != "NME_ACETYLATION") {
      LOG_ERROR("N-Terminal Variable PTM can only be NONE, M_ACETYLATION, NME or NME_ACETYLATION.");
      return false;
    }
  }

  std::string cutoff_spectral_type = arguments_["cutoffSpectralType"];
  if (cutoff_spectral_type != "EVALUE" && cutoff_spectral_type != "FDR") {
    LOG_ERROR("Spectrum-level cutoff type " << cutoff_spectral_type << " error! The value should be EVALUE|FDR");
    return false;
  }

  std::string cutoff_proteoform_type = arguments_["cutoffProteoformType"];
  if (cutoff_proteoform_type != "EVALUE" && cutoff_proteoform_type != "FDR") {
    LOG_ERROR("Proteoform-level cutoff type " << cutoff_proteoform_type << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_spectral_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Spectrum-level cutoff type "<< cutoff_spectral_type << " error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.");
    return false;
  }

  if (cutoff_proteoform_type == "FDR" && search_type != "TARGET+DECOY") {
    LOG_ERROR("Proteoform-level cutoff type "<< cutoff_proteoform_type << " error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.");
    return false;
  }

  std::string max_ptm_mass = arguments_["maxPtmMass"];
  try {
    double mass = std::stod(max_ptm_mass.c_str());
    if (mass <= 0.0) {
      LOG_ERROR("Maximum PTM mass " << max_ptm_mass << " error! The value should be positive.");
      return false;
    }
  }
  catch (int e) {
    LOG_ERROR("Maximum ptm mass " << max_ptm_mass << " should be a number.");
    return false;
  }

  std::string cutoff_spectral_value = arguments_["cutoffSpectralValue"];
  try {
    double th = std::stod(cutoff_spectral_value);
    if (th < 0) {
      LOG_ERROR("Spectrum-level cutoff value " << cutoff_spectral_value << " error! The value should be positive.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("Spectrum-level cutoff value " << cutoff_spectral_value << " should be a number.");
    return false;
  }

  std::string cutoff_proteoform_value = arguments_["cutoffSpectralValue"];
  try {
    double th = std::stod(cutoff_proteoform_value);
    if (th < 0) {
      LOG_ERROR("Proteoform-level cutoff value " << cutoff_proteoform_value << " error! The value should be positive.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("Proteoform-level cutoff value " << cutoff_proteoform_value << " should be a number.");
    return false;
  }

  std::string thread_number = arguments_["threadNumber"];
  try {
    int num = std::stoi(thread_number.c_str());
    if (num <= 0) {
      LOG_ERROR("Thread number " << thread_number << " error! The value should be positive.");
      return false;
    }
    int n = static_cast<int>(boost::thread::hardware_concurrency());
    if (num > n) {
      LOG_ERROR("Thread number " << thread_number << " error! The value is too large. Only " << n << " threads are supported.");
      return false;
    }
  }
  catch (int e) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }

  std::string mass_error_tole_value = arguments_["massErrorTolerance"];
  try {
    double tole = std::stoi(mass_error_tole_value);
    if (tole <= 0) {
      LOG_ERROR("Mass error tolerance: " << mass_error_tole_value << " error! The value should be positive.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("Mass error tolerance: " << mass_error_tole_value << " should be a number.");
    return false;
  }

  std::string form_error_tole_value = arguments_["proteoformErrorTolerance"];
  try {
    double tole = std::stod(form_error_tole_value);
    if (tole <= 0) {
      LOG_ERROR("PrSM clustering error tolerance: " << form_error_tole_value << " error! The value should be positive.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("PrSM clustering error tolerance: " << form_error_tole_value << " should be a number.");
    return false;
  }
  
  std::string ptm_num = arguments_["varPtmNumber"];
  try {
    double ptm_n = std::stod(ptm_num);
    if (ptm_n <= 0.0) {
      LOG_ERROR("Maximum PTM sites " << ptm_num << " error! The value should be positive.");
      return false;
    }
  }
  catch (int e) {
    LOG_ERROR("Maximum ptm sites " << ptm_num << " should be a number.");
    return false;
  }
  return true;
}
} /* namespace toppic */
