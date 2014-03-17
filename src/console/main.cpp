/*
 * main.cpp
 *
 *  Created on: Jan 22, 2014
 *      Author: xunlikun
 */

#include "base/fasta_reader.hpp"
#include "prsm/prsm_writer.hpp"
#include "prsm/prsm_combine.hpp"
#include "prsm/prsm_selector.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/evalue_processor.hpp"
#include "console/table_writer.hpp"
#include "console/output_selector.hpp"
#include "xpp/xml_generator.hpp"
#include "xpp/transformer.hpp"

using namespace prot;
using namespace std;

int main(){

		std::string config_file_name = "conf/configuration.xml";
		std::string search_db_file_name = "in/prot.fasta";
		std::string spectrum_file_name_ = "in/spectra_one.msalign";
//		std::string spectrum_file_name_ = "in/spectra.msalign";

		prot::ZeroPtmMngPtr mng_ptr = prot::ZeroPtmMngPtr(new prot::ZeroPtmMng (std::string(config_file_name)));

		mng_ptr->search_db_file_name_ = search_db_file_name;
		mng_ptr->spectrum_file_name_ = spectrum_file_name_;
		mng_ptr->output_file_ext_ = "ZERO_PTM_SEARCH";
		prot::zeroPtmSearchProcess(mng_ptr);

		PtmFastFilterMngPtr fast_filter_mng = PtmFastFilterMngPtr(new PtmFastFilterMng(config_file_name));
		fast_filter_mng->search_db_file_name_ = search_db_file_name;
		fast_filter_mng->spectrum_file_name_ = spectrum_file_name_;
		PtmFastFilterProcessorPtr processorb = PtmFastFilterProcessorPtr(new PtmFastFilterProcessor(fast_filter_mng));
		processorb->process();

		PtmMngPtr ptm_search_mng = PtmMngPtr(new PtmMng(config_file_name));
		ptm_search_mng->search_db_file_name_ = search_db_file_name;
		ptm_search_mng->spectrum_file_name_ = spectrum_file_name_;
		ptm_search_mng->input_file_ext_ ="_COMBINED";
		ptm_search_mng->output_file_ext_="PTM_SEARCH_RESULT";
		PtmProcessorPtr processor = PtmProcessorPtr(new PtmProcessor(ptm_search_mng));
		processor->process();

		std::vector<std::string> input_exts ;
		input_exts.push_back("PTM_SEARCH_RESULT");
		input_exts.push_back("ZERO_PTM_SEARCH");
		PrSMCombinePtr combine = PrSMCombinePtr(new PrSMCombine(search_db_file_name,
		                                                        spectrum_file_name_,
		                                                        input_exts,
		                                                        "RAW_SEARCH_RESULT"));
		combine->process();


		prot::TdgfMngPtr tdgf_mng = prot::TdgfMngPtr(new prot::TdgfMng (std::string(config_file_name)));
		tdgf_mng->search_db_file_name_ = search_db_file_name;
		tdgf_mng->spectrum_file_name_ = spectrum_file_name_;
		tdgf_mng->input_file_ext_ = "RAW_SEARCH_RESULT";
		tdgf_mng->output_file_ext_ = "EVALUED_RESULT";
		EValueProcessor evalue_processor(tdgf_mng);
		evalue_processor.init();
		/* compute E-value for a set of prsms each run */
		evalue_processor.process(false);

		PrSMSelector selector = PrSMSelector(search_db_file_name,
		                                     spectrum_file_name_,
		                                     "EVALUED_RESULT",
		                                     "TOP_RESULT",
		                                     1,
		                                     "e");
		selector.process();

		OutputSelector output_selector = OutputSelector(search_db_file_name,
		                                                spectrum_file_name_,
		                                                "TOP_RESULT",
		                                                "OUTPUT_RESULT",
		                                                "EVALUE",
		                                                0.01,
		                                                0.01,
		                                                mng_ptr->ppo_);
		output_selector.process();

		TableWriterPtr table_out = TableWriterPtr(new TableWriter(spectrum_file_name_,search_db_file_name,"OUTPUT_RESULT","OUTPUT_TABLE",mng_ptr->ppo_));
		table_out->write();




//		XmlGenerator xml_gene = XmlGenerator(spectrum_file_name_,search_db_file_name,"TOP_RESULT");
//		xml_gene.process();

//		TransformerPtr trans = TransformerPtr(new Transformer());
//		trans->trans();
//
}
