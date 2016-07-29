#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "base/file_util.hpp"
#include "base/xml_dom_util.hpp"
#include "feature/deconv_argu.hpp"

namespace prot {

DeconvArgument::DeconvArgument() {
  initArguments();
}

void DeconvArgument::initArguments() {
  arguments_["executiveDir"] = "";
  arguments_["spectrumFileName"] = "";
  arguments_["inputType"] = "mzXML";
  arguments_["outputType"] = "msalign";
  arguments_["refinePrecMass"]="true";
  arguments_["msLevel"] = "1";
  arguments_["missingLevelOne"] = "false";
  arguments_["maxCharge"] = "30";
  arguments_["maxMass"] = "100000";
  arguments_["mzError"] = "0.02";
  arguments_["snRatio"] = "1.0";
  arguments_["keepUnusedPeaks"] = "false";
  arguments_["outMultipleMass"] = "false";
  arguments_["precWindow"] = "2.0";
}

void DeconvArgument::outputArguments(std::ostream &output, 
                               std::map<std::string, std::string> arguments) {
  output << "********************** Parameters **********************" << std::endl;
  output << std::setw(44) << std::left << "Spectral file: " << "\t" << arguments["spectrumFileName"] << std::endl;
  output << "********************** Parameters **********************" << std::endl;
}

void DeconvArgument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: toppfd [options] spectrum-file-name" << std::endl; 
  std::cout << desc << std::endl; 
}

bool DeconvArgument::parse(int argc, char* argv[]) {
  std::string spectrum_file_name = "";
  std::string output_type = "";
  std::string max_charge = "";
  std::string max_mass = "";
  std::string mz_error = "";
  std::string sn_ratio = "";

  // Define and parse the program options
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options() 
        ("help,h", "Print the help message.") 
        ;
    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "") 
        ("output,o",po::value<std::string>(&output_type),"")
        ("level-one,l", "")
        ("keep,k", "")
        ("max-charge,c", po::value<std::string> (&max_charge), "")
        ("max_mass,m", po::value<std::string> (&max_mass), "")
        ("mz-error,e", po::value<std::string> (&mz_error), "")
        ("sn-ratio,s", po::value<std::string> (&sn_ratio), "")
        ("missing-level-one,n","")
        ("spectrum-file-name", po::value<std::string>(&spectrum_file_name)->required(), "Spectrum file name with its path.")
        ;

    po::positional_options_description positional_options;
    positional_options.add("spectrum-file-name", 1);


    std::string app_name;
    //= boost::filesystem::basename(argv[0]);
    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(),vm); 
      if ( vm.count("help") ) {
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
    arguments_["executiveDir"] = FileUtil::getExecutiveDir(argv_0);
    arguments_["spectrumFileName"] = spectrum_file_name;
    if (vm.count("output_type")) {
      arguments_["outputType"] = output_type;
    }

  }
  catch(std::exception&e ) {
    std::cerr << "Unhandled Exception in parsing command line"<<e.what()<<", application will now exit"<<std::endl;
    return false;
  }

  return validateArguments();
}

bool DeconvArgument::validateArguments() {
  return true;
}

} /* namespace prot */
