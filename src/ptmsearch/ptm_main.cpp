/*
 * ptm_main.cpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_processor.hpp"

namespace prot{
int main(){

	PtmMngPtr mng = PtmMngPtr(new PtmMng());
	mng->search_db_file_name_ = "in/prot.fasta";
	mng->spectrum_file_name_ = "in/sp.msalign";
//	mng->res_file_name_ ="";
	mng->input_file_ext_ ="FAST_FILTER_COMBINED";
	mng->output_file_ext_="PTM_SEARCH_RESULT";

//	mng->sp_para_->setActivation(prot::getActivationPtrByName(mng->base_data->getActivationPtrVec(),""));

	PtmProcessorPtr processor = PtmProcessorPtr(new PtmProcessor(mng));
	processor->process();
}
}
