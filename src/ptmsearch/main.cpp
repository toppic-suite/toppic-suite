///*
// * main.cpp
// *
// *  Created on: Jan 22, 2014
// *      Author: xunlikun
// */
//
//#include "base/fasta_reader.hpp"
//#include "prsm/prsm_writer.hpp"
//#include "zeroptmsearch/zero_ptm_mng.hpp"
//#include "zeroptmsearch/zero_ptm_search.hpp"
//#include "filterdiagonal/ptm_fast_filter_mng.hpp"
//#include "filterdiagonal/ptm_fast_filter_processor.hpp"
//#include "ptmsearch/ptm_mng.hpp"
//#include "ptmsearch/ptm_processor.hpp"
//#include "tdgf/tdgf_mng.hpp"
//
//using namespace prot;
//using namespace std;
//
//int main(){
//
//  std::string config_file_name = "conf/configuration.xml";
//  std::string search_db_file_name = "in/prot.fasta";
//  std::string spectrum_file_name_ = "in/sp.msalign";
//
//  PtmMngPtr ptm_search_mng = PtmMngPtr(new PtmMng(config_file_name));
//  ptm_search_mng->search_db_file_name_ = search_db_file_name;
//  ptm_search_mng->spectrum_file_name_ = spectrum_file_name_;
//  ptm_search_mng->input_file_ext_ ="_COMBINED";
//  ptm_search_mng->output_file_ext_="PTM_SEARCH_RESULT";
//  PtmProcessorPtr processor = PtmProcessorPtr(new PtmProcessor(ptm_search_mng));
//  processor->process();
//
//}
