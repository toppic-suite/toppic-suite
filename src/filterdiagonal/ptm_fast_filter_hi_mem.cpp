/*
 * ptm_fast_filter_hi_mem.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#include <algorithm>
#include <iostream>
#include "ptm_fast_filter_hi_mem.hpp"

namespace prot {

PtmFastFilterHiMem::PtmFastFilterHiMem(ProteoformPtrVec seqs,PtmFastFilterMngPtr mng){
	mng_ = mng;
	seqs_ = seqs;
	index_ = CompShiftHiMemPtr(new CompShiftHiMem(seqs,mng));
}
SimplePrSMPtrVec PtmFastFilterHiMem::getBestMatch(PrmMsPtr ms){
	SimplePrSMPtrVec2D matches = compute(ms);
//	std::cout<< matches.size() << std::endl;
	SimplePrSMPtrVec unique_match = sort(matches);
//	std::cout<< unique_match.size() << std::endl;
//	for(int i=0;i<unique_match.size();i++){
//		std::cout<<unique_match[i]->getSeqName()+"->"+convertToString(unique_match[i]->getScore())+" "<<std::endl;
//	}
	unsigned int num = mng_->ptm_fast_filter_result_num_;
	if(num > unique_match.size()){
		num = unique_match.size();
	}
	SimplePrSMPtrVec result;
	for(unsigned int i=0;i<num;i++){
		SimplePrSMPtr match = unique_match[i];
		if(match->getScore() > 0.0){
			result.push_back(match);
		}
		else{
			break;
		}
	}
//	std::cout<< result.size() << std::endl;
	return result;
}
SimplePrSMPtrVec2D PtmFastFilterHiMem::compute(PrmMsPtr ms){
	std::vector<std::vector<int>> masses = prot::getIntMassErrorList(ms,mng_->ptm_fast_filter_scale_,true,false);
	//logger
//	for(int i=0;i<masses[0].size();i++){
//		std::cout<<masses[0][i]<<std::endl;
//	}
	SimplePrSMPtrVec2D match;
//	std::cout<< "->" << std::endl;
	for(unsigned int i=0;i<masses[0].size();i++){
		std::vector<std::vector<int>> results =index_->compConvolution(masses[0],masses[1],i,mng_->ptm_fast_filter_result_num_);
		SimplePrSMPtrVec temp_match;
		for(unsigned int j =0;j <results.size();j++){
//			std::cout<<seqs_[results[j][0]]->getName()+"->"+prot::convertToString(results[j][1])<<std::endl;
			temp_match.push_back(SimplePrSMPtr(new SimplePrSM(ms->getHeaderPtr(),seqs_[results[j][0]],results[j][1])));
		}
		match.push_back(temp_match);
	}
	return match;
}
SimplePrSMPtrVec PtmFastFilterHiMem::sort(SimplePrSMPtrVec2D matches){
	SimplePrSMPtrVec sorted_match;

	for(unsigned int i=0;i<matches.size();i++){
		for(unsigned int j =0;j< matches[i].size();j++){
			sorted_match.push_back(matches[i][j]);
		}
	}

	//todo:sort method is not the same with java so that the result is different;
	std::sort(sorted_match.begin(),sorted_match.end(),simple_prsm_down);

	SimplePrSMPtrVec unique_match;
	for(unsigned int i=0;i< sorted_match.size();i++){
		bool found = false;
		std::string seq_name = sorted_match[i]->getSeqName();
		for(unsigned int j=0;j<unique_match.size();j++){
//			std::cout<< seq_name << std::endl;
			if(seq_name.compare(unique_match[j]->getSeqName())==0){
//				std::cout<< unique_match[j]->getSeqName() << std::endl;
				found=true;
				break;
			}
		}
		if(!found){
			unique_match.push_back(sorted_match[i]);
		}
	}

	//system.print?
	return unique_match;
}

} /* namespace prot */
