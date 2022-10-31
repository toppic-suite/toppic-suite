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

#include "boost/thread/thread.hpp"

#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/util/str_util.hpp"
#include "common/base/mod_util.hpp"
#include "common/util/version.hpp"
#include "common/util/mem_check.hpp"

#include "console/toppic_argument.hpp"

namespace toppic {

ToppicArgument::ToppicArgument() {
  arguments_ = initArguments();
}

std::map<std::string, std::string> ToppicArgument::initArguments() {
  std::map<std::string, std::string> arguments;
  arguments["oriDatabaseFileName"]="";
  arguments["databaseFileName"] = "";
  arguments["databaseBlockSize"] = "250000000";
  arguments["maxFragmentLength"] = "500";
  arguments["spectrumFileName"] = "";
  arguments["combinedOutputName"] = "";
  arguments["activation"] = "FILE";
  arguments["searchType"] = "TARGET";
  arguments["fixedMod"] = "";
  arguments["ptmNumber"] = "1";
  arguments["massErrorTolerance"] = "15";
  arguments["proteoformErrorTolerance"] = "1.2";
  arguments["cutoffSpectralType"] = "EVALUE";
  arguments["cutoffSpectralValue"] = "0.01";
  arguments["cutoffProteoformType"] = "EVALUE";
  arguments["cutoffProteoformValue"] = "0.01";
  arguments["allowProtMod"] = "NONE,NME,NME_ACETYLATION,M_ACETYLATION";
  arguments["numOfTopPrsms"] = "1";
  arguments["maxPtmMass"] = "500";
  arguments["minPtmMass"] = "-500";
  arguments["useLookupTable"] = "false";
  arguments["executiveDir"] = ".";
  arguments["resourceDir"] = "";
  arguments["keepTempFiles"] = "false";
  arguments["localThreshold"] = "0.15";
  arguments["groupSpectrumNumber"] = "1";
  arguments["filteringResultNumber"] = "20";
  arguments["residueModFileName"] = "";
  arguments["threadNumber"] = "1";
  arguments["useFeatureFile"] = "true";
  arguments["geneHTMLFolder"] = "true";
  arguments["keepDecoyResults"] = "false";
  arguments["nTermLabelMass"] = "0";
  arguments["version"] = "";
  return arguments;
}

void ToppicArgument::outputArguments(std::ostream &output, 
				     const std::string &sep,
                                     std::map<std::string, std::string> arguments) {
  int gap = 46;
  output << "********************** Parameters **********************" << std::endl;
  output << std::setw(gap) << std::left << "Protein database file:" << sep << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum file:" << sep << arguments["spectrumFileName"] << std::endl;
  output << std::setw(gap) << std::left << "Number of combined spectra:" << sep << arguments["groupSpectrumNumber"] << std::endl;
  output << std::setw(gap) << std::left << "Fragmentation method:" << sep << arguments["activation"] << std::endl;
  output << std::setw(gap) << std::left << "Search type:" << sep << arguments["searchType"] << std::endl;

  if (arguments["fixedMod"] != "") {
    //add fixed PTM information 
    if (arguments["fixedMod"] == "C57") {
      output << std::setw(gap) << std::left << "Fixed PTMs BEGIN" << std::endl;
      output << std::setw(gap) << std::left << "Carbamidomethylation" << sep << 57.021464 << sep << "C" << std::endl;
      output << std::setw(gap) << std::left << "Fixed PTMs END" << std::endl;
    }
    else if (arguments["fixedMod"] == "C58") {
      output << std::setw(gap) << std::left << "Fixed PTMs BEGIN" << std::endl;
      output << std::setw(gap) << std::left << "Carboxymethylation" << sep << 58.005479 << sep << "C" << std::endl;
      output << std::setw(gap) << std::left << "Fixed PTMs END" << std::endl;
    }
    else {
      output << std::setw(gap) << std::left << "Fixed PTMs file name:" << sep << arguments["fixedMod"] << std::endl;
      output << std::setw(gap) << std::left << "Fixed PTMs BEGIN" << std::endl;
      std::vector<std::vector<std::string>> mod_data = mod_util::readModTxtForTsv(arguments["fixedMod"]);
      for (size_t i = 0; i < mod_data.size(); i++) {
        output << std::setw(gap) << std::left << mod_data[i][0] << sep << mod_data[i][1] << sep << mod_data[i][2] << std::endl;
      }
      output << std::setw(gap) << std::left << "Fixed PTMs END" << std::endl;
    }
  }

  if (arguments["useFeatureFile"] == "true") {
    output << std::setw(gap) << std::left << "Use TopFD feature file:" << sep << "True" << std::endl;
  }
  else {
    output << std::setw(gap) << std::left << "Use TopFD feature file:" << sep << "False" << std::endl;
  }

  output << std::setw(gap) << std::left << "Maximum number of unexpected modifications:" << sep << arguments["ptmNumber"] << std::endl;
  output << std::setw(gap) << std::left << "Error tolerance for matching masses:" << sep << arguments["massErrorTolerance"] << " ppm" << std::endl;
  output << std::setw(gap) << std::left << "Error tolerance for identifying PrSM clusters:" << sep << arguments["proteoformErrorTolerance"] 
      << " Da" << std::endl;
  output << std::setw(gap) << std::left << "Spectrum-level cutoff type:" << sep << arguments["cutoffSpectralType"] << std::endl;
  output << std::setw(gap) << std::left << "Spectrum-level cutoff value:" << sep << arguments["cutoffSpectralValue"] << std::endl;
  output << std::setw(gap) << std::left << "Proteoform-level cutoff type:" << sep << arguments["cutoffProteoformType"] << std::endl;
  output << std::setw(gap) << std::left << "Proteoform-level cutoff value:" << sep << arguments["cutoffProteoformValue"] << std::endl;
  output << std::setw(gap) << std::left << "Allowed N-terminal forms:" << sep <<  arguments["allowProtMod"] << std::endl;
  output << std::setw(gap) << std::left << "Maximum mass shift of modifications:" << sep << arguments["maxPtmMass"] << " Da" << std::endl;
  output << std::setw(gap) << std::left << "Minimum mass shift of modifications:" << sep << arguments["minPtmMass"] << " Da" << std::endl;
  output << std::setw(gap) << std::left << "Thread number: " << sep << arguments["threadNumber"] << std::endl;

  if (arguments["useLookupTable"] == "true") {
    output << std::setw(gap) << std::left << "E-value computation:" << sep << "Lookup table" << std::endl;
  } else {
    output << std::setw(gap) << std::left << "E-value computation:" << sep << "Generating function" << std::endl;
  }

  if (arguments["residueModFileName"] != "") {
    output << std::setw(gap) << std::left << "Common modification file name:" << sep << arguments["residueModFileName"] << std::endl;
    output << std::setw(gap) << std::left <<  "PTMs for MIScore BEGIN" << std::endl;
    std::vector<std::vector<std::string>> mod_data = mod_util::readModTxtForTsv(arguments["residueModFileName"]);
    for (size_t i = 0; i < mod_data.size(); i++) {
      output << std::setw(gap) << std::left << mod_data[i][0] << sep << mod_data[i][1] << sep << mod_data[i][2] << std::endl;
    }
    output << std::setw(gap) << std::left <<  "PTMs for MIScore END" << std::endl;
    output << std::setw(gap) << std::left << "MIScore threshold:" << sep << arguments["localThreshold"] << std::endl;
  }

  output << std::setw(gap) << std::left << "Executable file directory:" << sep << arguments["executiveDir"] << std::endl;
  output << std::setw(gap) << std::left << "Start time:" << sep << arguments["startTime"] << std::endl;
  if (arguments["endTime"] != "") {
    output << std::setw(gap) << std::left << "End time:" << sep << arguments["endTime"] << std::endl;
  }
  output << std::setw(gap) << std::left << "Version:" << sep << arguments["version"] << std::endl;
  output << "********************** Parameters **********************" << std::endl;

}

std::string ToppicArgument::outputTsvArguments(std::map<std::string, 
                                               std::string> arguments) {
  std::stringstream output;
  outputArguments(output, "\t", arguments); 
  return output.str();
}

void ToppicArgument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: toppic [options] database-file-name spectrum-file-name" << std::endl; 
  std::cout << desc << std::endl; 
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool ToppicArgument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string argument_file_name = "";
  std::string activation = "";
  std::string fixed_mod = "";
  std::string allow_mod = "";
  std::string ptm_num = "";
  std::string mass_error_tole = "";
  std::string form_error_tole = "";
  std::string max_ptm_mass = "";
  std::string min_ptm_mass = "";
  std::string cutoff_spectral_type = "";
  std::string cutoff_spectral_value = "";
  std::string cutoff_proteoform_type = "";
  std::string cutoff_proteoform_value = "";
  std::string use_table = "";
  std::string group_num = "";
  std::string local_threshold = "";
  std::string filtering_result_num = "";
  std::string residue_mod_file_name = "";
  std::string thread_number = "";
  std::string combined_output_name = "";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options() 
        ("help,h", "Print the help message.") 
        ("activation,a", po::value<std::string>(&activation),
         "<CID|ETD|HCD|MPD|UVPD|FILE>. Fragmentation method of MS/MS spectra. When FILE is used, fragmentation methods of spectra are given in the input spectral data file. Default value: FILE.")
        ("fixed-mod,f", po::value<std::string> (&fixed_mod), 
         "<C57|C58|a fixed modification file>. Fixed modifications. Three available options: C57, C58, or the name of a text file containing the information of fixed modifications. When C57 is selected, carbamidomethylation on cysteine is the only fixed modification. When C58 is selected, carboxymethylation on cysteine is the only fixed modification.")
        ("n-terminal-form,n", po::value<std::string> (&allow_mod), 
         "<a list of allowed N-terminal forms>. N-terminal forms of proteins. Four N-terminal forms can be selected: NONE, NME, NME_ACETYLATION, and M_ACETYLATION. NONE stands for no modifications, NME for N-terminal methionine excision, NME_ACETYLATION for N-terminal acetylation after the initiator methionine is removed, and M_ACETYLATION for N-terminal methionine acetylation. When multiple forms are allowed, they are separated by commas. Default value: NONE,NME,NME_ACETYLATION,M_ACETYLATION.")
        ("decoy,d", "Use a shuffled decoy protein database to estimate false discovery rates.")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "<a positive integer>. Error tolerance for precursor and fragment masses in PPM. Default value: 15.")
        ("proteoform-error-tolerance,p", po::value<std::string> (&form_error_tole), "<a positive number>. Error tolerance for identifying PrSM clusters. Default value: 1.2 Dalton.")
        ("min-shift,m", po::value<std::string> (&min_ptm_mass), "Minimum value of the mass shift of an unexpected modification. Default value: -500 Dalton.")
        ("max-shift,M", po::value<std::string> (&max_ptm_mass), "Maximum value of the mass shift of an unexpected modification. Default value: 500 Dalton.")
        ("num-shift,s", po::value<std::string> (&ptm_num), "<0|1|2>. Maximum number of unexpected modifications in a proteoform-spectrum-match. Default value: 1.")
        ("spectrum-cutoff-type,t", po::value<std::string> (&cutoff_spectral_type), "<EVALUE|FDR>. Spectrum-level cutoff type for filtering identified proteoform-spectrum-matches. Default value: EVALUE.")
        ("spectrum-cutoff-value,v", po::value<std::string> (&cutoff_spectral_value), "<a positive number>. Spectrum-level cutoff value for filtering identified proteoform-spectrum-matches. Default value: 0.01.")
        ("proteoform-cutoff-type,T", po::value<std::string> (&cutoff_proteoform_type), "<EVALUE|FDR>. Proteoform-level cutoff type for filtering identified proteoform-spectrum-matches. Default value: EVALUE.")
        ("proteoform-cutoff-value,V", po::value<std::string> (&cutoff_proteoform_value), "<a positive number>. Proteoform-level cutoff value for filtering identified proteoform-spectrum-matches. Default value: 0.01.")
        ("lookup-table,l", "Use a lookup table method for computing p-values and E-values.")
        ("num-combined-spectra,r", po::value<std::string> (&group_num), 
         "<a positive integer>. Number of combined spectra. The parameter is set to 2 (or 3) for combining spectral pairs (or triplets) generated by the alternating fragmentation mode. Default value: 1")
        ("mod-file-name,i", po::value<std::string>(&residue_mod_file_name), "<a common modification file>. Specify a text file containing the information of common PTMs for characterization of PTMs in proteoform-spectrum-matches.")
        ("miscore-threshold,H", po::value<std::string> (&local_threshold), "<a positive number between 0 and 1>. Score threshold (modification identification score) for filtering results of PTM characterization. Default value: 0.15.")
        ("thread-number,u", po::value<std::string> (&thread_number), "<a positive integer>. Number of threads used in the computation. Default value: 1.")
        ("no-topfd-feature,x", "No TopFD feature file for proteoform identification.")
        ("combined-file-name,c", po::value<std::string>(&combined_output_name) , "Specify a file name for the combined spectrum data file and analysis results.")
        ("keep-temp-files,k", "Keep intermediate files.")
        ("keep-decoy-ids,K", "Keep decoy identifications.")
        ("skip-html-folder,g", "Skip the generation of HTML files for visualization.");

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
        ("min-shift,m", po::value<std::string> (&min_ptm_mass), "")
        ("num-shift,s", po::value<std::string> (&ptm_num), "")
        ("spectrum-cutoff-type,t", po::value<std::string> (&cutoff_spectral_type), "")
        ("spectrum-cutoff-value,v", po::value<std::string> (&cutoff_spectral_value), "")
        ("proteoform-cutoff-type,T", po::value<std::string> (&cutoff_proteoform_type), "")
        ("proteoform-cutoff-value,V", po::value<std::string> (&cutoff_proteoform_value), "")
        ("keep-temp-files,k", "")
        ("keep-decoy-ids,K", "")
        ("lookup-table,l", "")
        ("num-combined-spectra,r", po::value<std::string> (&group_num), "")
        ("thread-number,u", po::value<std::string> (&thread_number), "")
        ("no-topfd-feature,x", "")
        ("combined-file-name,c", po::value<std::string>(&combined_output_name) , "")
        ("mod-file-name,i", po::value<std::string>(&residue_mod_file_name), "")
        ("miscore-threshold,H", po::value<std::string> (&local_threshold), "")
        ("filtering-result-number", po::value<std::string>(&filtering_result_num), "Filtering result number. Default value: 20.")
        ("skip-html-folder,g","")
        ("database-file-name", po::value<std::string>(&database_file_name)->required(), "Database file name with its path.")
        ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(), "Spectrum file name with its path.");

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(),vm); 
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
    std::string argv_0 (argv[0]);
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
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target_decoy";
    } else {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target";
    }

    if (vm.count("fixed-mod")) {
      arguments_["fixedMod"] = fixed_mod;
    }

    if (vm.count("n-terminal-form")) {
      arguments_["allowProtMod"] = allow_mod;
    }    

    if (vm.count("num-shift")) {
      arguments_["ptmNumber"] = ptm_num;
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

    if (vm.count("min-shift")) {
      arguments_["minPtmMass"] = min_ptm_mass;
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

    if (vm.count("keep-temp-files")) {
      arguments_["keepTempFiles"] = "true";
    }
    if (vm.count("keep-decoy-ids")) {
      arguments_["keepDecoyResults"] = "true";
    }
    if (vm.count("lookup-table")) {
      arguments_["useLookupTable"] = "true";
    }

    if (vm.count("miscore-threshold")) {
      arguments_["localThreshold"] = local_threshold;
    }

    if (vm.count("num-combined-spectra")) {
      arguments_["groupSpectrumNumber"] = group_num;
    }

    if (vm.count("filtering-result-number")) {
      arguments_["filteringResultNumber"] = filtering_result_num;
    }

    if (vm.count("mod-file-name")) {
      arguments_["residueModFileName"] = residue_mod_file_name;
    }

    if (vm.count("thread-number")) {
      arguments_["threadNumber"] = thread_number;
    }

    if (vm.count("no-topfd-feature")) {
      arguments_["useFeatureFile"] = "false";
    }
    if (vm.count("skip-html-folder")) {
      arguments_["geneHTMLFolder"] = "false";
    }   
  }
  catch(std::exception&e ) {
    std::cerr << "Unhandled Exception in parsing command line" << e.what() << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool ToppicArgument::validateArguments() {
  if (!file_util::exists(arguments_["resourceDir"])) {
    LOG_ERROR("Resource direcotry " << arguments_["resourceDir"] << " does not exist!");
    return false;
  }

  if (!file_util::exists(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file " << arguments_["oriDatabaseFileName"] << " does not exist!");
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
      LOG_ERROR(spec_file_list_[k] << " does not exist!");
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

  std::string activation = arguments_["activation"];
  if(activation != "CID" && activation != "HCD" && activation != "MPD" 
     && activation != "ETD" && activation != "FILE" && activation != "UVPD") {
    LOG_ERROR("Activation type " << activation << " error! The value should be CID|ETD|HCD|MPD|UVPD|FILE!");
    return false;
  }

  std::string search_type = arguments_["searchType"];
  if(search_type != "TARGET" && search_type != "TARGET+DECOY"){
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

  if (cutoff_spectral_type == "FDR" && search_type != "TARGET+DECOY"){
    LOG_ERROR("Spectrum-level cutoff type "<< cutoff_spectral_type << " error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.");
    return false;
  }

  if (cutoff_proteoform_type == "FDR" && search_type != "TARGET+DECOY"){
    LOG_ERROR("Proteoform-level cutoff type "<< cutoff_proteoform_type << " error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.");
    return false;
  }

  std::string use_lookup_table = arguments_["useLookupTable"];
  if(use_lookup_table != "true" && use_lookup_table != "false"){
    LOG_ERROR("Use lookup_table " << use_lookup_table << " error! The value should be true|false!");
    return false;
  }

  if(use_lookup_table == "true" && arguments_["massErrorTolerance"] !="5" && arguments_["massErrorTolerance"]!="10" && arguments_["massErrorTolerance"]!="15"){
    LOG_ERROR("Error tolerance can only be 5, 10 or 15 when the lookup table approach for E-value computation is selected!");
    return false;
  }

  std::string max_ptm_mass = arguments_["maxPtmMass"];
  try {
    double mass = std::stod(max_ptm_mass.c_str());
    if(mass <= 0.0){
      LOG_ERROR("Maximum PTM mass " << max_ptm_mass << " error! The value should be positive.");
      return false;
    }
  }
  catch (const std::exception& ex) {
    LOG_ERROR("Maximum ptm mass " << max_ptm_mass << " should be a number.");
    return false;
  }

  std::string mass_error_tole_value = arguments_["massErrorTolerance"];
  try {
    double tole = std::stoi(mass_error_tole_value);
    if (tole <= 0) {
      LOG_ERROR("Mass error tolerance: " << mass_error_tole_value << " error! The value should be positive.");
      return false;
    }
  } 
  catch (const std::exception& ex) {
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
  } 
  catch (const std::exception& ex) {
    LOG_ERROR("PrSM clustering error tolerance: " << form_error_tole_value << " should be a number.");
    return false;
  }


  std::string cutoff_spectral_value = arguments_["cutoffSpectralValue"];
  try {
    double th = std::stod(cutoff_spectral_value);
    if (th < 0) {
      LOG_ERROR("Spectrum-level cutoff value " << cutoff_spectral_value << " error! The value should be positive.");
      return false;
    }
  } 
  catch (const std::exception& ex) {
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
  } 
  catch (const std::exception& ex) {
    LOG_ERROR("Proteoform-level cutoff value " << cutoff_proteoform_value << " should be a number.");
    return false;
  }

  std::string thread_number = arguments_["threadNumber"];
  try {
    int num = std::stoi(thread_number.c_str());
    bool valid = mem_check::checkThreadNum(num, "toppic");
    if (!valid) {
      return false;
    }
  } 
  catch (const std::exception& ex) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }

  return true;
}

}  // namespace toppic
