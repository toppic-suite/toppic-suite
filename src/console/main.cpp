/*
 * main.cpp
 *
 *  Created on: Jan 22, 2014
 *      Author: xunlikun
 */

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm_writer.hpp"

using namespace prot;
using namespace std;

int main(){

		std::string config_file_name = "conf/configuration.xml";
		std::string search_db_file_name = "in/prot.fasta";
		std::string spectrum_file_name_ = "in/sp.msalign";

//		prot::ZeroPtmMngPtr mng_ptr = prot::ZeroPtmMngPtr(new prot::ZeroPtmMng (std::string(config_file_name)));
//		mng_ptr->search_db_file_name_ = search_db_file_name;
//		mng_ptr->spectrum_file_name_ = spectrum_file_name_;
//		mng_ptr->output_file_ext_ = "ZREO_PTM_SEARCH";
//		prot::zeroPtmSearchProcess(mng_ptr);
//
//		PtmFastFilterMngPtr fast_filter_mng = PtmFastFilterMngPtr(new PtmFastFilterMng(config_file_name));
//		fast_filter_mng->search_db_file_name_ = search_db_file_name;
//		fast_filter_mng->spectrum_file_name_ = spectrum_file_name_;
//		PtmFastFilterProcessorPtr processorb = PtmFastFilterProcessorPtr(new PtmFastFilterProcessor(fast_filter_mng));
//		processorb->process();
//
//		PtmMngPtr ptm_search_mng = PtmMngPtr(new PtmMng(config_file_name));
//		ptm_search_mng->search_db_file_name_ = search_db_file_name;
//		ptm_search_mng->spectrum_file_name_ = spectrum_file_name_;
//		ptm_search_mng->input_file_ext_ ="_COMBINED";
//		ptm_search_mng->output_file_ext_="PTM_SEARCH_RESULT";
//		PtmProcessorPtr processor = PtmProcessorPtr(new PtmProcessor(ptm_search_mng));
//		processor->process();



		BaseDataPtr basedata = BaseDataPtr(new BaseData(config_file_name));
		ProteoformPtrVec proteoforms = prot::readFastaToProteoform(search_db_file_name,
		                                                       ResidueFactory::getBaseResiduePtrVec());
		PrSMPtrVec prsms = prot::readPrsm("in/sp.msalign.PTM_SEARCH_RESULT",proteoforms);
//		std::cout<<prsms[0]->getProteoformPtr()->getChangePtrVec().size()<<std::endl;
		PrSMWriterPtr all_writer= PrSMWriterPtr(new PrSMWriter("in/sp.msalign.PTM_SEARCH_RESULT_CHECK"));
		all_writer->writeVector(prsms);

//		prot::TdgfMngPtr tdgf_mng = prot::TdgfMngPtr(new prot::TdgfMng (std::string(config_file_name)));
//		tdgf_mng->search_db_file_name_ = search_db_file_name;
//		tdgf_mng->spectrum_file_name_ = spectrum_file_name_;
//		tdgf_mng->input_file_ext_ = "PTM_SEARCH_RESULT";
//		tdgf_mng->output_file_ext_ = "EVALUED_RESULT";
//		EValueProcessor evalue_processor = new EValueProcessor(tdgf_mng);
//		evalue_processor.init();
//		/* compute E-value for a set of prsms each run */
//		evalue_processor.process(false);


}
