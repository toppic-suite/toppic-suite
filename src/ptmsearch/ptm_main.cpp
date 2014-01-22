/*
 * ptm_main.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_processor.hpp"

using namespace prot;
int main(){

//	PtmFastFilterMngPtr mngb = PtmFastFilterMngPtr(new PtmFastFilterMng());
//	mngb->search_db_file_name_ = "in/prot.fasta";
//	mngb->spectrum_file_name_ = "in/sp.msalign";
//	PtmFastFilterProcessorPtr processorb = PtmFastFilterProcessorPtr(new PtmFastFilterProcessor(mngb));
//	processorb->process();

	PtmMngPtr mng = PtmMngPtr(new PtmMng("conf/configuration.xml"));
	mng->search_db_file_name_ = "in/prot.fasta";
	mng->spectrum_file_name_ = "in/sp.msalign";
//	mng->res_file_name_ ="";
	mng->input_file_ext_ ="_COMBINED";
	mng->output_file_ext_="PTM_SEARCH_RESULT";

//	mng->sp_para_->setActivation(prot::getActivationPtrByName(mng->base_data->getActivationPtrVec(),""));

	PtmProcessorPtr processor = PtmProcessorPtr(new PtmProcessor(mng));
	processor->process();
}
