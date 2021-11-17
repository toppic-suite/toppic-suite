//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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

#include "common/base/mod_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "common/util/version.hpp"
#include "common/util/mem_check.hpp"

#include "console/topmerge_argument.hpp"

namespace toppic {

Argument::Argument() {
  initArguments();
}

void Argument::initArguments() {
  arguments_["spectrumFileName"] = "";
  arguments_["combinedOutputName"] = "combined";
  arguments_["oriDatabaseFileName"] = "";
  arguments_["databaseFileName"] = "";
  arguments_["executiveDir"] = ".";
  arguments_["resourceDir"] = "";
  arguments_["maxPtmMass"] = "500";
  arguments_["minPtmMass"] = "-500";
  arguments_["proteoformErrorTolerance"] = "1.2";
  arguments_["useFeatureFile"] = "true";
  arguments_["searchType"] = "TARGET";
  arguments_["cutoffSpectralType"] = "EVALUE";
  arguments_["cutoffSpectralValue"] = "0.01";
  arguments_["cutoffProteoformType"] = "EVALUE";
  arguments_["cutoffProteoformValue"] = "0.01";
  arguments_["keepTempFiles"] = "false";
  arguments_["keepDecoyResults"] = "false";
  arguments_["localThreshold"] = "0.15";
  arguments_["residueModFileName"] = "";
  arguments_["version"] = "";
  arguments_["geneHTMLFolder"] = "true";
  // the following arguments are used for the initializatio of prsm para
  arguments_["groupSpectrumNumber"] = "1";
  arguments_["activation"] = "FILE";
  arguments_["massErrorTolerance"] = "0";
  arguments_["threadNumber"] = "1";
}

void Argument::outputArguments(std::ostream &output, 
                               std::map<std::string, std::string> arguments) {
  output << "********************** Parameters **********************" << std::endl;
  output << std::setw(44) << std::left << "Protein database file: " << "\t" << arguments["oriDatabaseFileName"] << std::endl;
  output << std::setw(44) << std::left << "Spectrum file: " << "\t" << arguments["spectrumFileName"] << std::endl;
  output << std::setw(44) << std::left << "Search type: " << "\t" << arguments["searchType"] << std::endl;

  if (arguments["useFeatureFile"] == "true") {
    output << std::setw(44) << std::left << "Use TopFD feature file: " << "\t" << "True" << std::endl;
  }
  else {
    output << std::setw(44) << std::left << "Use TopFD feature file: " << "\t" << "False" << std::endl;
  }
  output << std::setw(44) << std::left << "Error tolerance for identifying PrSM clusters: " << "\t" << arguments["proteoformErrorTolerance"] 
  << " Da" << std::endl;
  output << std::setw(44) << std::left << "Spectrum-level cutoff type: " << "\t" << arguments["cutoffSpectralType"] << std::endl;
  output << std::setw(44) << std::left << "Spectrum-level cutoff value: " << "\t" << arguments["cutoffSpectralValue"] << std::endl;
  output << std::setw(44) << std::left << "Proteoform-level cutoff type: " << "\t" << arguments["cutoffProteoformType"] << std::endl;
  output << std::setw(44) << std::left << "Proteoform-level cutoff value: " << "\t" << arguments["cutoffProteoformValue"] << std::endl;
  output << std::setw(44) << std::left << "Maximum mass shift of modifications: " << "\t" << arguments["maxPtmMass"] << " Da" << std::endl;
  output << std::setw(44) << std::left << "Minimum mass shift of modifications: " << "\t" << arguments["minPtmMass"] << " Da" << std::endl;
  if (arguments["residueModFileName"] != "") {
    output << std::setw(44) << std::left << "Common modification file name: " << "\t" << arguments["residueModFileName"] << std::endl;
    output << std::setw(44) << std::left <<  "PTMs for MIScore BEGIN" << std::endl;
    std::vector<std::vector<std::string>> mod_data = mod_util::readModTxtForTsv(arguments["residueModFileName"]);
    for (size_t i = 0; i < mod_data.size(); i++) {
      output << std::setw(44) << std::left << mod_data[i][0] << "\t" << mod_data[i][1] << "\t" << mod_data[i][2] << std::endl;
    }
    output << std::setw(44) << std::left <<  "PTMs for MIScore END" << std::endl;
    output << std::setw(44) << std::left << "MIScore threshold: " << "\t" << arguments["localThreshold"] << std::endl;
  }
  output << std::setw(44) << std::left << "Executable file directory: " << "\t" << arguments["executiveDir"] << std::endl;
  output << std::setw(44) << std::left << "Start time: " << "\t" << arguments["startTime"] << std::endl;
  if (arguments["endTime"] != "") {
    output << std::setw(44) << std::left << "End time: " << "\t" << arguments["endTime"] << std::endl;
  }
  output << std::setw(44) << std::left << "Version: " << "\t" << arguments["version"] << std::endl;

  output << "********************** Parameters **********************" << std::endl;
}

std::string Argument::outputTsvArguments(std::map<std::string, std::string> arguments) {
  std::stringstream output;
  outputArguments(output, arguments); 
  return output.str();
}

void Argument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: topmerge [options] database-file-name spectrum-file-name-1 spectrum-file-name-2 spectrum-file-name-3-or-more" << std::endl; 
  std::cout << desc << std::endl; 
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool Argument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string argument_file_name = "";
  std::string mass_error_tole = "";
  std::string form_error_tole = "";
  std::string max_ptm_mass = "";
  std::string min_ptm_mass = "";
  std::string cutoff_spectral_type = "";
  std::string cutoff_spectral_value = "";
  std::string cutoff_proteoform_type = "";
  std::string cutoff_proteoform_value = "";
  std::string local_threshold = "";
  std::string residue_mod_file_name = "";
  std::string combined_output_name = "";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options() 
        ("help,h", "Print the help message.") 
        ("decoy,d", "Use a shuffled decoy protein database to estimate false discovery rates.")
        ("proteoform-error-tolerance,p", po::value<std::string> (&form_error_tole), "<a positive number>. Error tolerance for identifying PrSM clusters. Default value: 1.2 Dalton.")
        ("min-shift,m", po::value<std::string> (&min_ptm_mass), "Minimum value of the mass shift of an unexpected modification. Default value: -500 Dalton.")
        ("max-shift,M", po::value<std::string> (&max_ptm_mass), "Maximum value of the mass shift of an unexpected modification. Default value: 500 Dalton.")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "<a positive integer>. Error tolerance for precursor and fragment masses in PPM. Default value: 15.")
        ("proteoform-error-tolerance,p", po::value<std::string> (&form_error_tole), "<a positive number>. Error tolerance for identifying PrSM clusters. Default value: 1.2 Dalton.")
        ("spectrum-cutoff-type,t", po::value<std::string> (&cutoff_spectral_type), "<EVALUE|FDR>. Spectrum-level cutoff type for filtering identified proteoform-spectrum-matches. Default value: EVALUE.")
        ("spectrum-cutoff-value,v", po::value<std::string> (&cutoff_spectral_value), "<a positive number>. Spectrum-level cutoff value for filtering identified proteoform-spectrum-matches. Default value: 0.01.")
        ("proteoform-cutoff-type,T", po::value<std::string> (&cutoff_proteoform_type), "<EVALUE|FDR>. Proteoform-level cutoff type for filtering identified proteoform-spectrum-matches. Default value: EVALUE.")
        ("proteoform-cutoff-value,V", po::value<std::string> (&cutoff_proteoform_value), "<a positive number>. Proteoform-level cutoff value for filtering identified proteoform-spectrum-matches. Default value: 0.01.")
        ("mod-file-name,i", po::value<std::string>(&residue_mod_file_name), "<a common modification file>. Specify a text file containing the information of common PTMs for characterization of PTMs in proteoform-spectrum-matches.")
        ("miscore-threshold,H", po::value<std::string> (&local_threshold), "<a positive number between 0 and 1>. Score threshold (modification identification score) for filtering results of PTM characterization. Default value: 0.15.")
        ("no-topfd-feature,x", "No TopFD feature file for proteoform identification.")
        ("combined-file-name,c", po::value<std::string>(&combined_output_name) , "Specify a file name for the combined spectrum data file and analysis results.")
        ("keep-temp-files,k", "Keep intermediate files.")        
        ("keep-decoy-ids,K", "Keep decoy identifications.")       
        ("skip-html-folder,g", "Skip the generation of HTML files for visualization.");

    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "") 
        ("no-topfd-feature,x", "")
        ("decoy,d", "")
        ("mass-error-tolerance,e", po::value<std::string> (&mass_error_tole), "")
        ("proteoform-error-tolerance,p", po::value<std::string> (&form_error_tole), "")
        ("max-shift,M", po::value<std::string> (&max_ptm_mass), "")
        ("min-shift,m", po::value<std::string> (&min_ptm_mass), "")
        ("spectrum-cutoff-type,t", po::value<std::string> (&cutoff_spectral_type), "")
        ("spectrum-cutoff-value,v", po::value<std::string> (&cutoff_spectral_value), "")
        ("proteoform-cutoff-type,T", po::value<std::string> (&cutoff_proteoform_type), "")
        ("proteoform-cutoff-value,V", po::value<std::string> (&cutoff_proteoform_value), "")
        ("keep-temp-files,k", "")             
        ("keep-decoy-ids,K", "")
        ("mod-file-name,i", po::value<std::string>(&residue_mod_file_name), "")
        ("miscore-threshold,H", po::value<std::string> (&local_threshold), "")
        ("database-file-name", po::value<std::string>(&database_file_name)->required(), "Database file name with its path.")
        ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(), "Spectrum file name with its path.")
        ("combined-file-name,c", po::value<std::string>(&combined_output_name) , "")
        ("skip-html-folder,g","");

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
    if (vm.count("decoy")) {
      arguments_["searchType"] = "TARGET+DECOY";
    }
    if (arguments_["searchType"] == "TARGET+DECOY") {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target_decoy";
    } else {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target";
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
    if (vm.count("keep-decoy-ids")) {
      arguments_["keepDecoyResults"] = "true";
    }

    if (vm.count("keep-temp-files")) {
      arguments_["keepTempFiles"] = "true";
    }
    if (vm.count("miscore-threshold")) {
      arguments_["localThreshold"] = local_threshold;
    }
    if (vm.count("mod-file-name")) {
      arguments_["residueModFileName"] = residue_mod_file_name;
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

bool Argument::validateArguments() {
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
  std::string search_type = arguments_["searchType"];
  if(search_type != "TARGET" && search_type != "TARGET+DECOY"){
    LOG_ERROR("Search type " << search_type << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
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
  return true;
}

}  // namespace toppic
