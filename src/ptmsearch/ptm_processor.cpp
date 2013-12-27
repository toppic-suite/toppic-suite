/*
 * ptm_processor.cpp
 *
 *  Created on: Dec 20, 2013
 *      Author: xunlikun
 */

#include "ptm_processor.hpp"
#include "prsm/simple_prsm.hpp"
#include "base/fasta_reader.hpp"
#include "spec/deconv_ms.hpp"

namespace prot {

PtmProcessor::PtmProcessor(PtmMngPtr mng){
	mng_ = mng;
	init();
}

void PtmProcessor::init(){
	seqs_ = prot::readFastaToProteoform(mng_->search_db_file_name_,mng_->base_data->getAcidPtrVec(),mng_->base_data->getResiduePtrVec());
	std::string sp_file_name = mng_->spectrum_file_name_;
	std::string simplePrsmFileName = mng_->spectrum_file_name_ + "." + mng_->input_file_ext_;
	simplePrsms_  = prot::readSimplePrSM(simplePrsmFileName.c_str());
	prsmFindSeq(simplePrsms_,seqs_);
}

void PtmProcessor::prsmFindSeq(SimplePrSMPtrVec simple_prsms,ProteoformPtrVec seqs){
	for(unsigned int i =0;i<simple_prsms.size();i++){
		simple_prsms[i]->findSeq(seqs);
	}
}

void PtmProcessor::process(){
	PtmSearcherPtr searcher = PtmSearcherPtr(new PtmSearcher(mng_));
	processDatabase(searcher);
}

void PtmProcessor::processDatabase(PtmSearcherPtr searcher){
	std::string sp_file_name = mng_->spectrum_file_name_;
	std::string output_file_name = sp_file_name+"."+mng_->output_file_ext_;

	//reader & writer

//	SimplePrSMPtrVec2D prsms;
//
//	for(unsigned int i =0;i<mng_->n_unknown_shift_;i++){
//		for(unsigned int j=0;j<4;j++){
//
//		}
//	}

	DeconvMsPtr deconv_sp;

}

} /* namespace prot */
