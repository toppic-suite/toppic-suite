/*
 * main.hpp
 *
 *  Created on: Jan 28, 2014
 *      Author: xunlikun
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "base/fasta_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "prsm/table_writer.hpp"
#include "prsm/output_selector.hpp"
#include "prsm/prsm_fdr.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"
#include "xpp/xml_generator.hpp"
#include "xpp/transformer.hpp"
#include "xpp/view_mng.hpp"
#include "xpp/folder_file.hpp"

using namespace prot;
using namespace std;
void showCommondList(){
//    cout <<"************************************************"<<endl;
//    cout <<"**                   usage                    **"<<endl;
//    cout <<"************************************************"<<endl;
    cout<<"USAGE:"<<endl;
    cout <<"-i, --arguments-filename <xml file name with path>"<<endl;
    cout <<"                         Deauft arguments file in which set the "<<
                                     "running arguments value"<<endl;
    cout <<"                         default value: conf/arguments.xml"<<endl;
    cout <<"-o, --configuration      <xml file name with path>"<<endl;
    cout <<"                         Deauft configuration file in which set the "<<
                                     "configuration file list"<<endl;
    cout <<"                         default value: conf/configuration.xml"<<endl;
//    cout <<"-o conf/configuration.xml or --configuration=conf/configuration.xml"<<endl;
    cout <<"-D, --database-filename  <fasta file name with path>"<<endl;
    cout <<"                         The database file that record unknown proteins"<<endl;
    cout <<"                         default value: in/prot.fasta"<<endl;
//    cout <<"-D in/prot.fasta or --databaseFileName=in/prot.fasta"<<endl;
    cout <<"-S, --spectra-filename   <msalign file name with path>"<<endl;
    cout <<"                         The spectrum file that record the spectral data"<<
                                     " we want to search"<<endl;
    cout <<"                         default value: in/spectra.msalign"<<endl;
//    cout <<"-S in/spectra.msalign or --spectrumFileName=in/spectra.msalign"<<endl;
    cout <<"-a, --activation         <CID|HCD|ETD|FILE>"<<endl;
    cout <<"                         The activation type and FILE means the activation "<<
                                     "type is given in spectral data file"<<endl;
    cout <<"                         default value: FILE"<<endl;
//    cout <<"# Activation can be CID, HCD, ETD or FILE (FILE means that the "
//        "activation type is given in the spectral data file)."<<endl;
//    cout <<"-a FILE or --activation=FILE"<<endl;
    cout <<"-s, --search-type        <TARGET|TARGET+DECOY>"<<endl;
    cout <<"                         When TARGET+DECOY is used, MS-Align+ will generate a "<<
                                     "scramble database from the protein database, and search"<<
                                     " spectra against the TARGET+DECOY database."<<endl;
    cout <<"                         default value: TARGET"<<endl;
//    cout <<"# SearchType can be TARGET or TARGET+DECOY. When TARGET+DECOY is "
//        "used, MS-Align+ will generate a scramble database from the protein "
//        "database, and search spectra against the TARGET+DECOY database."<<endl;
//    cout <<"-s TARGET or --searchType=TARGET"<<endl;
    cout <<"-c, --cysteineprotection <C0|C57|C58>"<<endl;
    cout <<"                         Cysteine protection group C0: no modification "<<
                                     "C57: Carbamidoemetylation or C58:Carboxymethylation"<<endl;
    cout <<"                         default value: C0"<<endl;
//    cout <<"#Cysteine protection group can be C0, C57 or C58. C0: no modification, "
//        "C57: Carbamidoemetylation or C58:Carboxymethylation"<<endl;
//    cout <<"-c C0 or --cysteineProtection=C0"<<endl;
    cout <<"-n, --shift-number       <int value>"<<endl;
    cout <<"                         Maximum number of unknown modifications "<<
                                     "configuration file list"<<endl;
    cout <<"                         default value: 2"<<endl;
//    cout <<"# Maximum number of modifications"<<endl;
//    cout <<"-n 2 or --shiftNumber=2"<<endl;
    cout <<"-e, --error-tolerance    <int value>"<<endl;
    cout <<"                         "<<
                                     "Error tolerance in PPM."<<endl;
    cout <<"                         default value: 15"<<endl;
//    cout <<"# Error tolerance in PPM"<<endl;
//    cout <<"-e 15 or --errorTolerance=15"<<endl;
    cout <<"-t, --cutoff-type        <EVALUE|FDR>"<<endl;
    cout <<"                         'FDR' can be used only "<<
                                     "if searchtype=TARGET+DECOY."<<endl;
    cout <<"                         default value: EVALUE"<<endl;
//    cout <<"#CutoffType can be EVALUE or FDR. 'FDR' can be used only if "
//        "searchtype=TARGET+DECOY. "<<endl;
//    cout <<"-u EVALUE or --cutoffType=EVALUE"<<endl;
    cout <<"-v, --cutoff-value       <EVALUE|FDR>"<<endl;
    cout <<"                         'When cutoffType = EValue, this cutoff value is for EVALUE."<<
                                     "When cutoff = FDR, this cutoff value is for FDR."<<endl;
    cout <<"                         default value: EVALUE"<<endl;
//    cout <<"# cutoff value. When cutoffType = EValue, this cutoff value "
//        "is for EVALUE. When cutoff = FDR, this cutoff value is for FDR."<<endl;
//    cout <<"-v 0.01 or --cutoff=0.01"<<endl;

}

void split(std::string& s, std::string& delim,std::vector< std::string >& ret)
{
	size_t last = 0;
	size_t index=s.find_first_of(delim,last);
	while (index!=std::string::npos){
		ret.push_back(s.substr(last,index-last));
		last=index+1;
		index=s.find_first_of(delim,last);
	}
	if (index-last>0)
	{
		ret.push_back(s.substr(last,index-last));
	}
}

void setArgumentsWithConfigFile(std::string filename,map<string,string> & arguments){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
    if(parser){
      XmlDOMDocument* doc = new XmlDOMDocument(parser, filename.c_str());
      if (doc) {
        xercesc::DOMElement* root = doc->getDocumentElement();
        arguments["activation"]=prot::getChildValue(root,"Fragmentation_Method",0);
        arguments["databaseFileName"]=prot::getChildValue(root,"Database_File",0);
        arguments["spectrumFileName"]=prot::getChildValue(root,"Spectrum_File",0);
        arguments["errorTolerance"]=prot::getChildValue(root,"Parent_and_Fragment_Mass_Error_Tolerance",0);
        arguments["cutoffType"]=prot::getChildValue(root,"Cutoff_Type",0);
        arguments["cutoff"]=prot::getChildValue(root,"Cutoff_Value",0);
        arguments["cysteineProtection"]=prot::getChildValue(root,"Cysteine_Protecting_Group",0);
        arguments["searchType"]=prot::getChildValue(root,"Decoy",0);
        arguments["shiftNumber"]=prot::getChildValue(root,"Number_of_Unexpected_PTMs",0);

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
        arguments["allProtMode"]=allow_mode;

//        std::cout<<arguments["activation"]<<
//            arguments["databaseFileName"]<<
//            arguments["spectrumFileName"]<<
//            arguments["errorTolerance"]<<
//            arguments["cutoffType"]<<
//            arguments["cutoff"]<<
//            arguments["cysteineProtection"]<<
//            arguments["searchType"]<<
//            arguments["shiftNumber"]<<
//            arguments["allProtMode"]<<endl;
      }
      delete doc;
    }
//	std::ifstream file;
//	file.open(filename.c_str(), std::ios::in);
//	std::string line;
//	string spliter = "=";
//	while (std::getline(file, line)) {
//		line = trim(line);
//		vector<string> commond_arr;
//		std::string temp = line;
//		split(temp,spliter,commond_arr);
//		std::cout<<commond_arr[0]<<"|"<<commond_arr[1]<<std::endl;
//		arguments[commond_arr[0]]=commond_arr[1];
//	}
}

map<string,string> getArugmentMap(){
	map<string,string> arguments;
	arguments["databaseFileName"]="in/prot.fasta";
	arguments["spectrumFileName"]="in/spectra.msalign";
	arguments["activation"]="FILE";
	arguments["searchType"]="TARGET";
	arguments["cysteineProtection"]="C0";
	arguments["shiftNumber"]="2";
	arguments["errorTolerance"]="15";
	arguments["cutoffType"]="EVALUE";
	arguments["cutoff"]="0.01";
	arguments["doOneDaltonCorrection"]="false";
	arguments["doChargeCorrection"]="false";
	return arguments;
}
bool findArgument(map<string,string> &arguments,string name){
	if(arguments[name].compare("")==0){
		return true;
	}
	return false;
}

void preProcess(map<string,string> arguments){
  if(arguments["searchType"].compare("TARGET+DECOY")==0){
      arguments["databaseFileName"]=arguments["databaseFileName"]+arguments["searchType"];
      ofstream onf;
      onf.open(arguments["databaseFileName"].c_str(), ios::out);
      FastaReader reader(arguments["databaseFileName"]);
//      std::vector<std::string> seq_info = reader.getNextSeq();
      FastaSeqPtr seq_info = reader.getNextSeq();
      while (seq_info!=nullptr) {
        std::string name = seq_info->getName();
        std::string seq = seq_info->getSeq();
        std::string temp = seq;
        std::string decoy_name = "DECOY_" + name;
        std::string decoy_seq = temp.substr(0, 2);
        temp = temp.substr(2, temp.length() - 2);
        for (unsigned int i = temp.length(); i > 0; i--) {
          int n = std::rand() % i;
          decoy_seq += temp.substr(n, 1);
          temp = temp.substr(0, n) + temp.substr(n + 1, temp.length() - n - 1);
        }
        onf << decoy_name << std::endl;
        onf << decoy_seq << std::endl;
        onf << name << std::endl;
        onf << seq << std::endl;
        seq_info = reader.getNextSeq();
      }
    }
}

void MsAlignPipeline(map<string,string> arguments){

  //pre_process
//  preProcess(arguments);
//
//  prot::ZeroPtmMngPtr mng_ptr = prot::ZeroPtmMngPtr(new prot::ZeroPtmMng (arguments));
//  mng_ptr->output_file_ext_ = "ZERO_PTM_SEARCH";
//  prot::zeroPtmSearchProcess(mng_ptr);
//
//  PtmFastFilterMngPtr fast_filter_mng = PtmFastFilterMngPtr(new PtmFastFilterMng(arguments));
//  PtmFastFilterProcessorPtr processorb = PtmFastFilterProcessorPtr(new PtmFastFilterProcessor(fast_filter_mng));
//  processorb->process();
//
//  PtmMngPtr ptm_search_mng = PtmMngPtr(new PtmMng(arguments));
//  ptm_search_mng->input_file_ext_ ="_COMBINED";
//  ptm_search_mng->output_file_ext_="PTM_SEARCH_RESULT";
//  PtmProcessorPtr processor = PtmProcessorPtr(new PtmProcessor(ptm_search_mng));
//  processor->process();
//
//  std::vector<std::string> input_exts ;
//  input_exts.push_back("PTM_SEARCH_RESULT");
//  input_exts.push_back("ZERO_PTM_SEARCH");
//  PrSMCombinePtr combine = PrSMCombinePtr(new PrSMCombine(arguments,
//                                                          input_exts,
//                                                          "RAW_SEARCH_RESULT"));
//  combine->process();
//
//
//  prot::TdgfMngPtr tdgf_mng = prot::TdgfMngPtr(new prot::TdgfMng (arguments));
//  tdgf_mng->input_file_ext_ = "RAW_SEARCH_RESULT";
//  tdgf_mng->output_file_ext_ = "EVALUED_RESULT";
//  EValueProcessor evalue_processor(tdgf_mng);
//  evalue_processor.init();
//  /* compute E-value for a set of prsms each run */
//  evalue_processor.process(false);
//
//  if(arguments["searchType"].compare("TARGET+DECOY")==0){
//	PrSMSelector selector = PrSMSelector(arguments,
//	                                     "EVALUED_RESULT",
//	                                     "TOP",
//	                                     1);
//	selector.process();
//
//    prot::PrSMFdr fdr(arguments, "TOP", "TOP_RESULT");
//    fdr.process();
//  }
//  else{
//    PrSMSelector selector = PrSMSelector(arguments,
//                                       "EVALUED_RESULT",
//                                       "TOP_RESULT",
//                                       1);
//    selector.process();
//  }
//
//  OutputSelector output_selector = OutputSelector(arguments,
//                                                  "TOP_RESULT",
//                                                  "OUTPUT_RESULT");
//  output_selector.process();
//
//  TableWriterPtr table_out = TableWriterPtr(new TableWriter(arguments,"OUTPUT_RESULT","OUTPUT_TABLE"));
//  table_out->write();


  createFolder("xml/species");
  createFolder("xml/prsms");
  createFolder("xml/proteins");
  XmlGenerator xml_gene = XmlGenerator(arguments,"OUTPUT_RESULT");
  xml_gene.process();

  createFolder("html/species");
  createFolder("html/prsms");
  createFolder("html/proteins");
  copyFile("etc/FreeMono.ttf","html/FreeMono.ttf",true);
  copyFile("etc/sorttable.js","html/sorttable.js",true);

    TransformerPtr trans = TransformerPtr(new Transformer());
    trans->trans();

}

bool checkFile(std::string filename){
  fstream file;
    file.open(filename, ios::in);
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
  file_to_check.push_back("conf/acid.xml");
  file_to_check.push_back("conf/activation.xml");
  file_to_check.push_back("conf/arguments.xml");
  file_to_check.push_back("conf/configuration.xml");
  file_to_check.push_back("conf/fix_mod_residue.xml");
  file_to_check.push_back("conf/ion_type.xml");
  file_to_check.push_back("conf/neutral_loss.xml");
  file_to_check.push_back("conf/prot_mod.xml");
  file_to_check.push_back("conf/ptm.xml");
  file_to_check.push_back("conf/residue.xml");
  file_to_check.push_back("conf/support_peak_type.xml");
  file_to_check.push_back("conf/trunc.xml");
  bool flag = true;
  for(unsigned int i=0;i<file_to_check.size();i++){
    if(!checkFile(path+file_to_check[i])){
      LOG_DEBUG("Can not find "<<file_to_check[i]<<" in path "<<path);
      flag=false;
    }
  }
  return flag;
}

std::string getExePath(std::string & command_path) {
  std::string run_path = "";
  size_t pos = command_path.find_last_of("/");
  if (pos == std::string::npos) {
    command_path = "";
  } else {
    command_path = command_path.substr(0, pos + 1);
  }
  std::cout << command_path << std::endl;
  fstream file;
  file.open(command_path + "conf/acid.xml", ios::in);
  if (!file) {
    std::string os_path = std::getenv("PATH");
    size_t pos = os_path.find(";", 0);
    if (pos == std::string::npos) {
      //it's unix
      std::vector<std::string> paths = prot::split(os_path, ':');
      for (unsigned int i = 0; i < paths.size(); i++) {
        if (paths[i].find_last_of('/')==paths[i].length()-1) {
          paths[i] = paths[i].substr(0, paths[i].length() - 1);
        }
        std::cout<<paths[i]<<std::endl;
        file.close();
        file.open(paths[i] + "conf/acid.xml", ios::in);
        if (file) {
          if(checkPath(paths[i])){
            run_path = paths[i];
            break;
          }
        }
        file.close();
      }
    } else {
      //it's windows
      std::vector<std::string> paths = prot::split(os_path, ';');
      for (unsigned int i = 0; i < paths.size(); i++) {
        if (paths[i].find_last_of('/')==paths[i].length()-1) {
          paths[i] = paths[i].substr(0, paths[i].length() - 1);
        }
        if (paths[i].find_last_of('\\') == paths[i].length() - 1) {
          paths[i] = paths[i].substr(0, paths[i].length() - 1);
        }
        file.open(paths[i] + "conf/acid.xml", ios::in);
        if (file) {
          if (checkPath(paths[i])) {
            run_path = paths[i];
            break;
          }
        }
        file.close();
      }
    }
  }
  return run_path;
}

int main(int argc, char* argv[]){
//  std::cout<<argv[0]<<std::endl;
  std::string command_path = argv[0];
  std::string run_path = getExePath(command_path);
  if(run_path.compare("")==0){
    LOG_ERROR("Can't find the configuration folder: conf or "+command_path+"conf or &PATH/conf;Please check if the folder exist!");
  }

  map<string,string> arguments = getArugmentMap();
  for(int i=1;i<argc;i++){
    std::string command = argv[i];
    std::string argument = "-argumentsFileName";

    size_t index = command.find(argument);
    if(index!=std::string::npos){
      setArgumentsWithConfigFile(argv[i+1],arguments);
    }
  }
    string spliter_cpp = "=";
      for(int i=1;i<argc;i++){
          string command = argv[i];

          if(i==1 && (command.compare("--help")==0||command.compare("-h")==0)){
             showCommondList();
          }
          else{

          int last = 0;
          size_t index=command.find_first_of("-",last);
          if(index!=std::string::npos){
            command = command.substr(1,command.length()-1);
          }
          else{
            cout<< "argument error! argument name:"<<command<<" is wrong!"<<endl;
            cout<< "see all argument ,use '-h or --help'"<<endl;
            return 1;
          }
          index=command.find_first_of("-",last);
          if(index==std::string::npos){
            if(i<argc){
              i++;
              command = command+"=" + argv[i];
            }
            else{
              cout<< "argument error! number of argument is:"<<argc<<"!"<<endl;
              cout<< "see all argument ,use '-h or --help'"<<endl;
              return 1;
            }
          }
          else {
            command = command.substr(1,command.length()-1);
          }
            vector<string> command_arr;
            split(command,spliter_cpp,command_arr);
            std::cout<<command_arr[0]<<":"<<command_arr[1]<< std::endl;
            if(command_arr.size() == 1){
              cout<< "argument error! argument:"<<command_arr[0]<<" have no value!"<<endl;
              cout<< "see all argument ,use '-h or --help'"<<endl;
              return 1;
            }
            if(findArgument(arguments,command_arr[0])){
              cout<< "argument error! argument name:"<<command_arr[0]<<" is wrong!"<<endl;
              cout<< "see all argument ,use '-h or --help'"<<endl;
              return 1;
            }
            if(command_arr.size() >2){
              cout<< "argument error! argument:"<<command_arr[0]<<" have wrong format!"<<endl;
              cout<< "see all argument ,use '-h or --help'"<<endl;
              return 1;
            }
            arguments[command_arr[0]]=command_arr[1];
          }
      }

//              std::cout<<arguments["activation"]<<
//                  arguments["databaseFileName"]<<
//                  arguments["spectrumFileName"]<<
//                  arguments["errorTolerance"]<<
//                  arguments["cutoffType"]<<
//                  arguments["cutoff"]<<
//                  arguments["cysteineProtection"]<<
//                  arguments["searchType"]<<
//                  arguments["shiftNumber"]<<
//                  arguments["allProtMode"]<<endl;

//      MsAlignPipeline(arguments);

      return 0;
}




