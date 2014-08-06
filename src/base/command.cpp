
#include "base/command.hpp"

namespace prot {

command::command() {
  initOption();
  initArguments();
}

command::~command() {
}

int command::setCommandValue(std::string command,std::string value){
  if(validateCommandValue(command,value)){
    return -1;
  }
  std::string argument_key = options_[command];
  arguments_[argument_key]=value;
  return 0;
}

void command::setArgumentsWithConfigFile(std::string filename){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
    if(parser){
      XmlDOMDocument* doc = new XmlDOMDocument(parser, filename.c_str());
      if (doc) {
        xercesc::DOMElement* root = doc->getDocumentElement();
        arguments_["activation"]=prot::getChildValue(root,"Fragmentation_Method",0);
        arguments_["databaseFileName"]=prot::getChildValue(root,"Database_File",0);
        arguments_["spectrumFileName"]=prot::getChildValue(root,"Spectrum_File",0);
        arguments_["errorTolerance"]=prot::getChildValue(root,"Parent_and_Fragment_Mass_Error_Tolerance",0);
        arguments_["cutoffType"]=prot::getChildValue(root,"Cutoff_Type",0);
        arguments_["cutoff"]=prot::getChildValue(root,"Cutoff_Value",0);
        arguments_["cysteineProtection"]=prot::getChildValue(root,"Cysteine_Protecting_Group",0);
        arguments_["searchType"]=prot::getChildValue(root,"Decoy",0);
        arguments_["shiftNumber"]=prot::getChildValue(root,"Number_of_Unexpected_PTMs",0);

        xercesc::DOMElement* allow_prot_mod_list = prot::getChildElement(root,"Protein_N_Terminal_Variable_PTM_List",0);
        int allow_prot_node_number = prot::getChildCount(allow_prot_mod_list,"Protein_N_Terminal_Variable_PTM");
        std::string allow_mode="";
        for(int i=0;i<allow_prot_node_number;i++){
          if(i==0){
            allow_mode = prot::getChildValue(allow_prot_mod_list,"Protein_N_Terminal_Variable_PTM",i);
          }
          else{
            allow_mode = allow_mode+","+prot::getChildValue(allow_prot_mod_list,"Protein_N_Terminal_Variable_PTM",i);
          }
        }
        arguments_["allProtMode"]=allow_mode;
      }
      delete doc;
    }
}

bool command::validateCommandValue(std::string command,std::string value){
  std::string argument = options_[command];
  if(argument.compare("")==0){
    LOG_ERROR("Can not find the command!");
    return false;
  }
  if (argument.compare("argumentFileName") == 0
      || argument.compare("configuration") == 0
      || argument.compare("databaseFileName") == 0
      || argument.compare("spectrumFileName") == 0) {
    if (!existFile(value)) {
      LOG_ERROR("Argument "<<command<<" error! The file :"<<value<<" not exist!");
      return false;
    }
  }
  else if(argument.compare("activation") == 0){
    if(value.compare("CID")!=0 && value.compare("HCD")!=0
        && value.compare("ETD")!=0 && value.compare("FILE")!=0){
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be CID|HCD|ETD|FILE!");
      return false;
    }
  }
  else if(argument.compare("searchType") == 0){
    if(value.compare("TARGET")!=0 && value.compare("TARGET+DECOY")!=0){
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be TARGET|TARGET+DECOY!");
      return false;
    }
  }
  else if(argument.compare("cysteineProtection") == 0){
    if(value.compare("C0")!=0 && value.compare("C57")!=0
        && value.compare("C58")!=0){
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be C0|C57|C58!");
      return false;
    }
  }
  else if(argument.compare("cutoffType") == 0){
    if (value.compare("EVALUE") != 0 && value.compare("FDR") != 0) {
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be EVALUE|FDR");
      return false;
    }
    if (value.compare("FDR") == 0 && arguments_["searchType"].compare("TARGET+DECOY") != 0){
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be EVALUE beacuse of the value"
          <<" of s|search-type is TARGET. please check the argument s|search-type's value");
      return false;
    }
  }
  else if(argument.compare("shiftNumber") == 0){
    if (value.compare("0") != 0 && value.compare("1") != 0
        && value.compare("2") != 0) {
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be 0|1|2!");
      return false;
    }
  }
  else if(argument.compare("errorTolerance") == 0){
    int ppm = std::stoi(value);
    if (ppm < 0 || ppm > 100) {
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be [0,100]");
      return false;
    }
    std::string result = convertToString(ppm);
    if(result.compare(value)!=0){
      LOG_ERROR("Argument "<<command<<" error! The value should be int");
      return false;
    }
  }
  else if(argument.compare("cutoff") == 0){
    double th = std::stod(value);
    if(th<0){
      LOG_ERROR("Argument "<<command<<" error! The value :"<<value<<" should be positive");
      return false;
    }
    std::string result = convertToString(th);
    if(result.compare(value)!=0){
      LOG_ERROR("Argument "<<command<<" error! The value should be double");
      return false;
    }
  }
  return true;
}

void command::initOption() {
  options_["-i"] = "argumentFileName";
  options_["--arguments-filename"] = "argumentFileName";
  options_["-o"] = "configuration";
  options_["--configuration"] = "configuration";
  options_["-D"] = "databaseFileName";
  options_["--database-filename"] = "databaseFileName";
  options_["-S"] = "spectrumFileName";
  options_["--spectra-filename"] = "spectrumFileName";
  options_["-a"] = "activation";
  options_["--activation"] = "activation";
  options_["-s"] = "searchType";
  options_["--search-type"] = "searchType";
  options_["-c"] = "cysteineProtection";
  options_["--cysteineprotection"] = "cysteineProtection";
  options_["-n"] = "shiftNumber";
  options_["--shift-number"] = "shiftNumber";
  options_["-e"] = "errorTolerance";
  options_["--error-tolerance"] = "errorTolerance";
  options_["-t"] = "cutoffType";
  options_["--cutoff-type"] = "cutoffType";
  options_["-v"] = "cutoff";
  options_["--cutoff-value"] = "cutoff";
}

void command::initArguments() {
  arguments_["argumentFileName"] = "conf/arguments.xml";
  arguments_["configuration"] = "conf/configuration.xml";
  arguments_["databaseFileName"] = "in/prot.fasta";
  arguments_["spectrumFileName"] = "in/spectra.msalign";
  arguments_["activation"] = "FILE";
  arguments_["searchType"] = "TARGET";
  arguments_["cysteineProtection"] = "C0";
  arguments_["shiftNumber"] = "2";
  arguments_["errorTolerance"] = "15";
  arguments_["cutoffType"] = "EVALUE";
  arguments_["cutoff"] = "0.01";
  arguments_["allProtMode"] = "NONE";
}

void showCommondList(){
    std::cout<<"USAGE:"<<std::endl;
//    std::cout<<"-i, --arguments-filename <xml file name with path>"<<std::endl;
//    std::cout<<"                         Deauft arguments file in which set the "<<
//                                     "running arguments value"<<std::endl;
//    std::cout<<"                         default value: conf/arguments.xml"<<std::endl;
//    std::cout<<"-o, --configuration      <xml file name with path>"<<std::endl;
//    std::cout<<"                         Deauft configuration file in which set the "<<
//                                     "configuration file list"<<std::endl;
//    std::cout<<"                         default value: conf/configuration.xml"<<std::endl;
    std::cout<<"-D, --database-filename  <fasta file name with path>"<<std::endl;
    std::cout<<"                         The database file that record unknown proteins"<<std::endl;
    std::cout<<"                         default value: in/prot.fasta"<<std::endl;
    std::cout<<"-S, --spectra-filename   <msalign file name with path>"<<std::endl;
    std::cout<<"                         The spectrum file that record the spectral data"<<
                                     " we want to search"<<std::endl;
    std::cout<<"                         default value: in/spectra.msalign"<<std::endl;
    std::cout<<"-a, --activation         <CID|HCD|ETD|FILE>"<<std::endl;
    std::cout<<"                         The activation type and FILE means the activation "<<
                                     "type is given in spectral data file"<<std::endl;
    std::cout<<"                         default value: FILE"<<std::endl;
    std::cout<<"-s, --search-type        <TARGET|TARGET+DECOY>"<<std::endl;
    std::cout<<"                         When TARGET+DECOY is used, MS-Align+ will generate a "<<
                                     "scramble database from the protein database, and search"<<
                                     " spectra against the TARGET+DECOY database."<<std::endl;
    std::cout<<"                         default value: TARGET"<<std::endl;
    std::cout<<"-c, --cysteineprotection <C0|C57|C58>"<<std::endl;
    std::cout<<"                         Cysteine protection group C0: no modification "<<
                                     "C57: Carbamidoemetylation or C58:Carboxymethylation"<<std::endl;
    std::cout<<"                         default value: C0"<<std::endl;
    std::cout<<"-n, --shift-number       <int value>"<<std::endl;
    std::cout<<"                         Maximum number of unknown modifications "<<
                                     "configuration file list"<<std::endl;
    std::cout<<"                         default value: 2"<<std::endl;
    std::cout<<"-e, --error-tolerance    <int value>"<<std::endl;
    std::cout<<"                         "<<
                                     "Error tolerance in PPM."<<std::endl;
    std::cout<<"                         default value: 15"<<std::endl;
    std::cout<<"-t, --cutoff-type        <EVALUE|FDR>"<<std::endl;
    std::cout<<"                         'FDR' can be used only "<<
                                     "if searchtype=TARGET+DECOY."<<std::endl;
    std::cout<<"                         default value: EVALUE"<<std::endl;
    std::cout<<"-v, --cutoff-value       <EVALUE|FDR>"<<std::endl;
    std::cout<<"                         'When cutoffType = EValue, this cutoff value is for EVALUE."<<
                                     "When cutoff = FDR, this cutoff value is for FDR."<<std::endl;
    std::cout<<"                         default value: EVALUE"<<std::endl;

}

bool existFile(std::string filename){
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

bool checkPath(std::string path){
  std::vector<std::string> file_to_check;
  file_to_check.push_back("/conf/acid.xml");
  file_to_check.push_back("/conf/activation.xml");
  file_to_check.push_back("/conf/configuration.xml");
  file_to_check.push_back("/conf/ion_type.xml");
  file_to_check.push_back("/conf/neutral_loss.xml");
  file_to_check.push_back("/conf/prot_mod.xml");
  file_to_check.push_back("/conf/ptm.xml");
  file_to_check.push_back("/conf/residue.xml");
  file_to_check.push_back("/conf/support_peak_type.xml");
  file_to_check.push_back("/conf/trunc.xml");
  bool flag = true;
  for(size_t i=0;i<file_to_check.size();i++){
    if(!existFile(path+file_to_check[i])){
      LOG_ERROR("Can not find "<<file_to_check[i]<<" in path "<<path);
      flag=false;
    }
  }
  return flag;
}

std::string getExePath(std::string & command_path) {
  std::string run_path = "########";
  size_t pos = command_path.find_last_of("/");
  if (pos == std::string::npos) {
    command_path = "";
  } else {
    command_path = command_path.substr(0, pos);
  }
  std::fstream file;
  file.open(command_path + "/conf/acid.xml", std::ios::in);
  if (!file) {
    std::string os_path = std::getenv("PATH");
    size_t pos = os_path.find(";", 0);
    if (pos == std::string::npos) {
      //it's unix
      std::vector<std::string> paths = prot::split(os_path, ':');
      for (size_t i = 0; i < paths.size(); i++) {
        if (paths[i].find_last_of('/') == paths[i].length() - 1) {
          paths[i] = paths[i].substr(0, paths[i].length() - 1);
        }
        file.close();
        file.open(paths[i] + "/conf/acid.xml", std::ios::in);
        if (file) {
          if (checkPath(paths[i])) {
            run_path = paths[i];
            break;
          }
        }
        file.close();
      }

    } else {
      //it's windows
      std::vector<std::string> paths = prot::split(os_path, ';');
      for (size_t i = 0; i < paths.size(); i++) {
        if (paths[i].find_last_of('/') == paths[i].length() - 1) {
          paths[i] = paths[i].substr(0, paths[i].length() - 1);
        }
        if (paths[i].find_last_of('\\') == paths[i].length() - 1) {
          paths[i] = paths[i].substr(0, paths[i].length() - 1);
        }
        file.open(paths[i] + "/conf/acid.xml", std::ios::in);
        if (file) {
          if (checkPath(paths[i])) {
            run_path = paths[i];
            break;
          }
        }
        file.close();
      }
    }
  } else {
    run_path = command_path;
  }
  //Debug in eclipse
  if(command_path.compare("")!=0)
  {
    std::string up_path = command_path;
    size_t up_pos = up_path.find_last_of("/");
    if (up_pos == std::string::npos) {
      up_path = "";
    } else {
      up_path = up_path.substr(0, up_pos);
    }
//    std::cout<<command_path<<std::endl;
//    std::cout<<up_path<<std::endl;
//    std::cout<<up_path<<"/conf/acid.xml"<<std::endl;
    if(existFile(up_path+"/conf/acid.xml")){
      run_path=up_path;
    }
  }
  return run_path;
}

int getOS(){
  std::string os_path = std::getenv("PATH");
//  std::cout<<os_path<<std::endl;
  size_t pos_win_symbol = os_path.find(";", 0);
  size_t pos_win_home = os_path.find("\\Windows\\System32\\", 0);

  size_t pos_uni_symbol = os_path.find(";", 0);
  size_t pos_uni_home = os_path.find("/usr/", 0);

  if(pos_win_symbol != std::string::npos && pos_win_home!= std::string::npos){
    return 1;
  }
  if(pos_uni_symbol == std::string::npos && pos_uni_home!= std::string::npos){
    return 2;
  }
  return 0;
}

int runCommand(std::string cmd){
  char psBuffer[128];
  FILE *pPipe;
  if( (pPipe = popen( cmd.c_str(),"r" )) == NULL )
  return 0;
  while(fgets(psBuffer, 128, pPipe))
  printf("%s",psBuffer);
  int rc = pclose(pPipe);
  return rc;
}

int runCommand(std::string cmd,std::string mod){
  char psBuffer[128];
  FILE *pPipe;
  if( (pPipe = popen( cmd.c_str(), mod.c_str() )) == NULL )
  return 0;
  while(fgets(psBuffer, 128, pPipe))
  printf("%s",psBuffer);
  int rc = pclose(pPipe);
  return rc;
}

std::string SystemInfo::exe_path_;
std::string SystemInfo::os_;
std::string SystemInfo::html_path_;
std::string SystemInfo::xml_path_;

void SystemInfo::initSystemInfo(std::string argument_zero){
  std::string command_path = argument_zero;
    std::string run_path = getExePath(command_path);
    if(run_path.compare("########")==0){
      LOG_ERROR("Can't find the configuration folder: conf or "+command_path+"conf or &PATH/conf;Please check if the folder exist!");
      return;
    }
    SystemInfo::exe_path_=run_path;
    int type = getOS();
    if(type==1){
      os_="win";
    }
    else if(type==2){
      os_="uni";
    }
    else{
      os_="unk";
    }
}

} /* namespace prot */
