#include "base/file_util.hpp"
#include "console/argument.hpp"

namespace prot {

Argument::Argument() {
  initArguments();
}

void Argument::initArguments() {
  arguments_["oriDatabaseFileName"]="";
  arguments_["databaseFileName"] = "";
  arguments_["databaseBlockSize"] = "1000000";
  arguments_["spectrumFileName"] = "";
  arguments_["activation"] = "FILE";
  arguments_["searchType"] = "TARGET";
  arguments_["cysteineProtection"] = "C0";
  arguments_["shiftNumber"] = "2";
  arguments_["errorTolerance"] = "15";
  arguments_["cutoffType"] = "EVALUE";
  arguments_["cutoffValue"] = "0.01";
  arguments_["allowProtMod"] = "NONE,NME,NME_ACETYLATION";
  arguments_["numOfTopPrsms"] = "1";
  arguments_["maxPtmMass"] = "1000000";
  arguments_["executiveDir"] = ".";
  arguments_["logFileName"] = "";
  arguments_["keepTempFiles"] = "false";
  arguments_["fullBinaryPath"] = "false";
  arguments_["useTable"] = "true";
}

void Argument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: toppic [options] database-file-name spectrum-file-name" << std::endl; 
  std::cout << desc << std::endl; 
}

void Argument::setArgumentsByConfigFile(const std::string &filename){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if(parser){
    XmlDOMDocument* doc = new XmlDOMDocument(parser, filename.c_str());
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      arguments_["oriDatabaseFileName"]=getChildValue(root,"database_file_name",0);
      arguments_["spectrumFileName"]=getChildValue(root,"spectrum_file_name",0);
      arguments_["logFileName"]=getChildValue(root,"log_file_name",0);
      arguments_["activation"]=getChildValue(root,"fragmentation_method",0);
      arguments_["cysteineProtection"]=getChildValue(root,"cysteine_protecting_group",0);
      arguments_["searchType"]=getChildValue(root,"search_type",0);
      arguments_["shiftNumber"]=getChildValue(root,"shift_number",0);
      arguments_["errorTolerance"]=getChildValue(root,"error_tolerance",0);
      arguments_["cutoffType"]=getChildValue(root,"cutoff_type",0);
      arguments_["cutoffValue"]=getChildValue(root,"cutoff_value",0);
      arguments_["maxPtmMass"]=getChildValue(root,"max_ptm_mass",0);
      arguments_["useTable"]=getChildValue(root,"use_table",0);

      xercesc::DOMElement* prot_mod_list = getChildElement(root,"protein_variable_ptm_list",0);
      int allow_prot_node_number = getChildCount(prot_mod_list,"protein_variable_ptm");
      std::string allow_mod="";
      for(int i=0;i<allow_prot_node_number;i++){
        if(i==0){
          allow_mod = getChildValue(prot_mod_list,"protein_variable_ptm",i);
        }
        else{
          allow_mod = allow_mod+","+getChildValue(prot_mod_list,"protein_variable_ptm",i);
        }
      }
      arguments_["allowProtMod"]=allow_mod;
    }
    delete doc;
  }
}

bool Argument::parse(int argc, char* argv[]) {
  std::string database_file_name = "";
  std::string spectrum_file_name = "";
  std::string argument_file_name = "";
  std::string activation = "";
  std::string protection = "";
  std::string shift_num = "";
  std::string error_tole = "";
  std::string max_ptm_mass = "";
  std::string cutoff_type = "";
  std::string cutoff_value = "";
  std::string log_file_name = "";
  std::string use_table = "";

  /** Define and parse the program options*/
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options() 
        ("help,h", "Print the help message.") 
        ("activation,a", po::value<std::string>(&activation),
         "<CID|HCD|ETD|FILE>. Activation type of tandem mass spectra. When FILE is used, the activation type information is given in the input spectral data file. Default value: FILE.")
        ("cysteine-protection,c", po::value<std::string> (&protection), 
         "<C0|C57|C58>. Cysteine protecting group: C0: no modification, C57: Carbamidoemetylation, or C58:Carboxymethylation. Default value: C0.")
        ("decoy,d", "Use a decoy protein database to estimate false discovery rates.")
        ("error-tolerance,e", po::value<std::string> (&error_tole), "<positive integer value>. Error tolerance for precursor and fragment masses in PPM. Default value: 15.")
        ("max-ptm,m", po::value<std::string> (&max_ptm_mass), "<positive double value>. Maximum absolute value of masses (in Dalton) of unexpected post-translational modifications in proteoforms. Default value: 1000000.")
        ("ptm-number,p", po::value<std::string> (&shift_num), "<0|1|2>. Maximum number of unexpected post-translational modifications in a proteoform-spectrum-match. Default value: 2.")
        ("cutoff-type,t", po::value<std::string> (&cutoff_type), "<EVALUE|FDR>. Cutoff type for reporting protein-spectrum-matches. Default value: EVALUE.")
        ("cutoff-value,v", po::value<std::string> (&cutoff_value), "<positive double value>. Cutoff value for reporting protein-spectrum-matches. Default value: 0.01.")
        ("use-table,u", po::value<std::string> (&use_table), "<true|false>. Use precomputed tables to estimate E-values. When the tables are used, the valid error tolerances are only 5, 10 and 15. Default value: true.");
    po::options_description desc("Options");

    desc.add_options() 
        ("help,h", "Print the help message.") 
        ("argument-file,f",po::value<std::string>(&argument_file_name),"Argument file name.")
        ("activation,a", po::value<std::string>(&activation),
         "<CID|HCD|ETD|FILE>. The activation type of tandem mass spectra. When FILE is used, the activation type information is given in spectral data file. Default value: FILE.")
        ("cysteine-protection,c", po::value<std::string> (&protection), 
         "<C0|C57|C58>. Cysteine protecting group: C0: no modification, C57: Carbamidoemetylation, or C58:Carboxymethylation. Default value: C0.")
        ("decoy,d", "Use a decoy protein database to estimate false discovery rates.")
        ("error-tolerance,e", po::value<std::string> (&error_tole), "<int value>. Error tolerance of precursor and fragment masses in PPM. Default value: 15.")
        ("max-ptm,m", po::value<std::string> (&max_ptm_mass), "<positive double value>. Maximum absolute value (in Dalton) of the masses of unexpected PTMs in the identified proteoform. Default value: 1000000.")
        ("ptm-number,p", po::value<std::string> (&shift_num), "<0|1|2>. Maximum number of unexpected PTMs. Default value: 2.")
        ("cutoff-type,t", po::value<std::string> (&cutoff_type), "<EVALUE|FDR>. Cutoff value type for reporting protein-spectrum-matches. Default value: EVALUE.")
        ("cutoff-value,v", po::value<std::string> (&cutoff_value), "<positive double value>. Cutoff value for reporting protein-spectrum-matches. Default value: 0.01.")
        ("log-file-name,l", po::value<std::string>(&log_file_name), "Log file name with its path.")
        ("keep-temp-files,k", "Keep temporary files.")
        ("use-table,u", po::value<std::string> (&use_table), "<true|false>. Use precomputed tables to estimate p-values. If used, the error tolerance can only by 5, 10 or 15. Default value: true.")
        ("full-binary-path,b", "Full binary path.")
        ("database-file-name", po::value<std::string>(&database_file_name)->required(), "Database file name with its path.")
        ("spectrum-file-name", po::value<std::string>(&spectrum_file_name)->required(), "Spectrum file name with its path.");
        

    po::positional_options_description positional_options;
    positional_options.add("database-file-name", 1);
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
    if (vm.count("full-binary-path")) {
      arguments_["executiveDir"] = argv[0];
    }
    else {
      arguments_["executiveDir"] = getExecutiveDir(argv_0);
    }
    LOG_DEBUG("Executive Dir " << arguments_["ExecutiveDir"]);
    if (vm.count("argument-file")) {
      setArgumentsByConfigFile(argument_file_name);
    }
    arguments_["oriDatabaseFileName"] = database_file_name;
    arguments_["spectrumFileName"] = spectrum_file_name;
    if (vm.count("activation")) {
      arguments_["activation"] = activation;
    }
    if (vm.count("decoy")) {
      arguments_["searchType"] = "TARGET+DECOY";
    }
    if (arguments_["searchType"] == "TARGET+DECOY") {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target_decoy";
    }
    else {
      arguments_["databaseFileName"]=arguments_["oriDatabaseFileName"] + "_target";
    }
    if (vm.count("cysteine-protection")) {
      arguments_["cysteineProtection"] = protection;
    }
    if (vm.count("shift-number")) {
      arguments_["shiftNumber"] = shift_num;
    }
    if (vm.count("error-tolerance")) {
      arguments_["errorTolerance"] = error_tole;
    }
    if (vm.count("max-ptm")) {
      arguments_["maxPtmMass"] = max_ptm_mass;
    }
    if (vm.count("cutoff-type")) {
      arguments_["cutoffType"] = cutoff_type;
    }
    if (vm.count("cutoff-value")) {
      arguments_["cutoffValue"] = cutoff_value;
    }
    if (vm.count("log-file-name")) {
      arguments_["logFileName"] = log_file_name;
    }
    if (vm.count("keep-temp-files")) {
      arguments_["keepTempFiles"] = "true";
    }
    if (vm.count("full-binary-path")) {
      arguments_["fullBinaryPath"] = "true";
    }
    if (vm.count("use-table")) {
      arguments_["useTable"] = use_table;
    }
  }
  catch(std::exception&e ) {
    std::cerr << "Unhandled Exception in parsing command line"<<e.what()<<", application will now exit"<<std::endl;
    return false;
  }
  
  /*
  std::cout <<"*** Parameters begin ***" << std::endl;
  for(std::map<std::string, std::string>::const_iterator it = arguments_.begin();
      it != arguments_.end(); ++it) {
        std::cout << it->first << " " << it->second << std::endl;
  }
  std::cout <<"*** Parameters end ***" << std::endl;
  */

  return validateArguments();
}

bool Argument::validateArguments() {
  if (!testFile(arguments_["oriDatabaseFileName"])) {
    LOG_ERROR("Database file " << arguments_["databaseFileName"] << " does not exist!");
    return false;
  }

  if (!testFile(arguments_["spectrumFileName"])) {
    LOG_ERROR("Spectrum file " << arguments_["spectrumFileName"] << " does not exist!");
    return false;
  }
  std::string activation = arguments_["activation"];
  if(activation != "CID" && activation != "HCD" 
     && activation != "ETD" && activation != "FILE") {
    LOG_ERROR("Activation type " << activation << " error! The value should be CID|HCD|ETD|FILE!");
    return false;
  }
  std::string protection = arguments_["cysteineProtection"];
  if(protection != "C0" && protection != "C57" && protection != "C58") {
    LOG_ERROR("Cysteine protection group " << protection <<" error! The value should be C0|C57|C58!");
    return false;
  }
  std::string search_type = arguments_["searchType"];
  if(search_type != "TARGET" && search_type != "TARGET+DECOY"){
    LOG_ERROR("Search type " << search_type << " error! The value should be TARGET|TARGET+DECOY!");
    return false;
  }
  std::string shift_number = arguments_["shiftNumber"];
  if (shift_number != "0" && shift_number != "1" && shift_number != "2") {
    LOG_ERROR("PTM number "<< shift_number <<" error! The value should be 0|1|2!");
    return false;
  }
  
  std::string cutoff_type = arguments_["cutoffType"];
  if (cutoff_type != "EVALUE" && cutoff_type != "FDR") {
    LOG_ERROR("Cutoff type " << cutoff_type << " error! The value should be EVALUE|FDR");
    return false;
  }

  if (cutoff_type == "FDR" && search_type != "TARGET+DECOY"){
    LOG_ERROR("Cutoff type "<< cutoff_type << " error! FDR cutoff cannot be used when no decoy database is used! Please add argument '-d' in the command.");
    return false;
  }

  std::string use_table = arguments_["useTable"];
  if(use_table != "true" && use_table != "false"){
    LOG_ERROR("Use table " << use_table << " error! The value should be true|false!");
    return false;
  }

  if(use_table == "true" && arguments_["errorTolerance"] !="5" && arguments_["errorTolerance"]!="10" && arguments_["errorTolerance"]!="15"){
    LOG_ERROR("Error tolerance can only be 5, 10 or 15 when using precomputed tables!");
    return false;
  }

  std::string max_ptm_mass = arguments_["maxPtmMass"];
  try {
    double mass = atof(max_ptm_mass.c_str());
    if(mass <= 0.0){
      LOG_ERROR("Maximum PTM mass " << max_ptm_mass << " error! The value should be positive.");
      return false;
    }
  }
  catch (int e) {
    LOG_ERROR("Maximum ptm mass " << max_ptm_mass << " should be a number.");
    return false;
  }
  return true;

  std::string cutoff_value = arguments_["cutoffValue"];
  try {
    double th = atof(cutoff_value.c_str());
    if(th<0){
      LOG_ERROR("Cutoff value " << cutoff_value << " error! The value should be positive.");
      return false;
    }
  }
  catch (int e) {
    LOG_ERROR("Cutoff value " << cutoff_value << " should be a number.");
    return false;
  }
  return true;
}

bool testFile(std::string filename){
  std::fstream file;
    file.open(filename, std::ios::in);
    if (!file) {
      file.close();
      return false;
    }
    else{
      file.close();
      return true;
    }
}

} /* namespace prot */
