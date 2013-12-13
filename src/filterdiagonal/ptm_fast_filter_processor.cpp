/*
 * ptm_fast_filter_processor.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <sstream>

#include "ptm_fast_filter_processor.hpp"
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_writer.hpp"

namespace prot {

PtmFastFilterProcessor::PtmFastFilterProcessor(PtmFastFilterMngPtr mng){
	mng_ = mng;
	//todo:xunlikun@config file name
	BaseDataPtr base_data = BaseDataPtr(new BaseData("config_file_name"));
	ProteoformPtrVec proteoforms = readFastaToProteoform(mng_->search_db_file_name_,base_data->getAcidPtrVec(),base_data->getResiduePtrVec());
	filter_ = PtmFastFilterBlockPtr(new PtmFastFilterBlock(proteoforms,mng_));
}

void PtmFastFilterProcessor::process(){
	//todo:xunlikun@config file name
	BaseDataPtr base_data = BaseDataPtr(new BaseData("config_file_name"));

	std::string sp_file_name = mng_->spectrum_file_name_;
	int n_spectrum = prot::countSpNum(sp_file_name.c_str(),base_data->getActivationPtrVec());
	for(int i=0;i<filter_->getBlockSize();i++){
		processBlock(i,sp_file_name,n_spectrum);
	}
	combineBlock(sp_file_name);
}

void PtmFastFilterProcessor::processBlock(int block,std::string sp_file_name,int n_spectra){
	//system.out
	filter_->initBlock(block);
	//todo:xunlikun@config file name
	BaseDataPtr base_data = BaseDataPtr(new BaseData("config_file_name"));
	MsAlignReader reader(sp_file_name.c_str(), base_data->getActivationPtrVec());
	std::stringstream block_s;
	block_s<<block;
	std::string output_file_name = mng_->spectrum_file_name_ + "." + mng_->output_file_ext_+"_"+block_s.str();
//	SimplePrSMWriter prsm_writer(sp_file_name.c_str());
	DeconvMsPtr deconv_sp;
	int cnt = 0;
	while((deconv_sp = reader.getNextMs()) != nullptr){
		cnt++;
//		for(unsigned int i =0;i<deconv_sp->size();i++){
			SpectrumSetPtr spectrum_set = prot::getSpectrumSet(deconv_sp,0,mng_->sp_para_,prot::getProtModPtrByName(base_data->getProtModPtrVec(),"ACETYLATION")->getPepShift());
			if(spectrum_set != nullptr){
				std::string scan = deconv_sp->getHeaderPtr()->getScansString();
				SimplePrSMPtrVec matches = filter_->getBestMathBatch(spectrum_set);
				//writer.write(matches);
			}
//		}
	}
	reader.close();
	//writer.close();
	//system.out
}

void PtmFastFilterProcessor::combineBlock(std::string sp_file_name){
	//system.out
	SimplePrSMPtrVec2D matches;
	//todo:xunlikun@config file name
	BaseDataPtr base_data = BaseDataPtr(new BaseData("config_file_name"));

	//pravate readsimplePrsm
	for(int i=0;i<filter_->getBlockSize();i++){
		std::stringstream block_s;
		block_s<<i;
		std::string block_file_name = mng_->spectrum_file_name_+ "." + mng_->output_file_ext_+"_"+block_s.str();
		matches.push_back(prot::readSimplePrSM(block_file_name.c_str()));
	}

	MsAlignReader reader(sp_file_name.c_str(), base_data->getActivationPtrVec());
	std::string output_file_name = mng_->spectrum_file_name_ + "." + mng_->output_file_ext_+"_COMBINED";
	//writer
	DeconvMsPtr deconv_sp;
	while((deconv_sp = reader.getNextMs()) != nullptr){
		//private getBestMatch
		SimplePrSMPtrVec selected_matche;
		for(unsigned int i=0;i<matches.size();i++){
			for(unsigned int j=0;i<matches[i].size();j++){
				if(matches[i][j]->isMatch(deconv_sp->getHeaderPtr())){
					selected_matche.push_back(matches[i][j]);
				}
			}
		}
		//writer;
	}

	reader.close();
	//writer.close
	//system.out
}

} /* namespace prot */
