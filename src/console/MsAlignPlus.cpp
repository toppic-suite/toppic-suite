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
#include "console/command.hpp"

using namespace prot;
using namespace std;

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

int main(int argc, char* argv[]){
//  std::cout<<argv[0]<<std::endl;
  cmdPtr cmd = cmdPtr(new command());
  std::string command_path = argv[0];
  std::string run_path = getExePath(command_path);
  if(run_path.compare("########")==0){
    LOG_ERROR("Can't find the configuration folder: conf or "+command_path+"conf or &PATH/conf;Please check if the folder exist!");
  }

  for(int i=1;i<argc;i++){
    int flag=0;
    std::string command = argv[i];
    std::vector<std::string> command_array = prot::split(command,'=');
    if(command_array[0].compare("-i")==0){
      flag = cmd->setCommandValue(command_array[0],argv[i+1]);
      cmd->setArgumentsWithConfigFile(argv[i+1]);
    }
    if(command_array[0].compare("--arguments-filename")==0){
      flag = cmd->setCommandValue(command_array[0],command_array[1]);
      cmd->setArgumentsWithConfigFile(command_array[1]);
    }
    if(flag==-1){
      return -1;
    }
  }
    string spliter_cpp = "=";
      for(int i=1;i<argc;i++){
          string command = argv[i];
          if (i == 1  && (command.compare("--help") == 0 || command.compare("-h") == 0)) {
              showCommondList();
          }
          else if(command.length()==2){
            cmd->setCommandValue(command,argv[i+1]);
            i++;
          }
          else{
            std::vector<std::string> command_array = prot::split(command,'=');
            if(command_array.size()>2){
              LOG_ERROR("argument error,may be it contains too much '='");
              return -1;
            }
            if(command_array.size()==2){
              int reuslt = cmd->setCommandValue(command_array[0],command_array[1]);
              if(reuslt == -1){
                return -1;
              }
            }
          }
      }

      std::map<std::string,std::string>::iterator it;
      for(it=cmd->getArgument().begin();it!=cmd->getArgument().end();++it){
        cout<<"key:"<<it->first<<" -> value:"<<it->second<<endl;
      }

//      MsAlignPipeline(arguments);

      return 0;
}




