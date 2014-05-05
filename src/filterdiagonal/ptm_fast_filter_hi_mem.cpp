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

PtmFastFilterHiMem::PtmFastFilterHiMem(ProteoformPtrVec seqs,
                                       PtmFastFilterMngPtr mng){
    mng_ = mng;
    seqs_ = seqs;
    index_ = CompShiftHiMemPtr(new CompShiftHiMem(seqs,mng));
}

SimplePrsmPtrVec PtmFastFilterHiMem::getBestMatch(PrmMsPtr ms){
    SimplePrsmPtrVec2D matches = compute(ms);
    SimplePrsmPtrVec unique_match = sort(matches);
    unsigned int num = mng_->ptm_fast_filter_result_num_;
    if(num > unique_match.size()){
        num = unique_match.size();
    }
    SimplePrsmPtrVec result;
    for(unsigned int i=0;i<num;i++){
        SimplePrsmPtr match = unique_match[i];
        if(match->getScore() > 0.0){
            result.push_back(match);
        }
        else{
            break;
        }
    }
    return result;
}
SimplePrsmPtrVec2D PtmFastFilterHiMem::compute(PrmMsPtr ms){
    std::vector<std::vector<int>> masses 
      = prot::getIntMassErrorList(ms,mng_->ptm_fast_filter_scale_,true,false);
    SimplePrsmPtrVec2D match;
    for(unsigned int i=0;i<masses[0].size();i++){
        std::vector<std::vector<int>> results 
        =index_->compConvolution(masses[0],masses[1],i,
                                 mng_->ptm_fast_filter_result_num_);
        SimplePrsmPtrVec temp_match;
        for(unsigned int j =0;j <results.size();j++){
            temp_match.push_back(
          SimplePrsmPtr(new SimplePrsm(ms->getHeaderPtr(),
                                       seqs_[results[j][0]],
                                       results[j][1])));
        }
        match.push_back(temp_match);
    }
    return match;
}
SimplePrsmPtrVec PtmFastFilterHiMem::sort(SimplePrsmPtrVec2D matches){
    SimplePrsmPtrVec sorted_match;

    for(unsigned int i=0;i<matches.size();i++){
        for(unsigned int j =0;j< matches[i].size();j++){
            sorted_match.push_back(matches[i][j]);
        }
    }

    std::sort(sorted_match.begin(),sorted_match.end(),simplePrsmDown);

    SimplePrsmPtrVec unique_match;
    for(unsigned int i=0;i< sorted_match.size();i++){
        bool found = false;
        std::string seq_name = sorted_match[i]->getSeqName();
        for(unsigned int j=0;j<unique_match.size();j++){
            if(seq_name.compare(unique_match[j]->getSeqName())==0){
                found=true;
                break;
            }
        }
        if(!found){
            unique_match.push_back(sorted_match[i]);
        }
    }
    return unique_match;
}

} /* namespace prot */
