/*
 * diagonal.cpp
 *
 *  Created on: Jan 1, 2014
 *      Author: xunlikun
 */

#include "ptmsearch/diagonal.hpp"
#include "base/logger.hpp"
#include "float.h"
#include "limits.h"

namespace prot {
DiagonalHeaderPtrVec refineHeadersBgnEnd(
        int first_pos,
        ProteoformPtr seq,
        DeconvMsPtr deconv_ms,
        ExtendMsPtr ms_three,
        PtmMngPtr mng,
        DiagonalHeaderPtrVec headers){
    DiagonalHeaderPtrVec result_list;
    for(unsigned int i=0;i<headers.size();i++){
        TheoPeakPtrVec ions = prot::getProteoformTheoPeak(
                seq,
                deconv_ms->getHeaderPtr()->getActivationPtr(),
                mng->sp_para_->getMinMass());
        int bgn = headers[i]->getMatchFirstResPos()-first_pos;
        int end = headers[i]->getMatchLastResPos()+1-first_pos;
        PeakIonPairPtrVec pairs = prot::findPairs(ms_three,ions,bgn,end);

        if(pairs.size()<1){
            int pair_size = pairs.size();
            LOG_WARN("Empty Segment is found "+prot::convertToString(pair_size));
        }
        else{
            int new_bgn = first_pos + getNewBgn(pairs);
            int new_end = first_pos + getNewEnd(pairs);
            headers[i]->setMatchFirstResPos(new_bgn);
            headers[i]->setMatchLastResPos(new_end);
            result_list.push_back(headers[i]);
        }
    }
    return result_list;
}

int getNewBgn(PeakIonPairPtrVec pairs){
    int newBgn = INT_MAX;
    for(unsigned int i=0;i<pairs.size();i++)
    {
        if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos() < newBgn){
            newBgn = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
        }
    }
    return newBgn;
}
int getNewEnd(PeakIonPairPtrVec pairs){
    int newEnd = 0;
    for(unsigned int i=0;i<pairs.size();i++)
    {
        if(pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos() > newEnd){
            newEnd = pairs[i]->getTheoPeakPtr()->getIonPtr()->getPos();
        }
    }
    return newEnd;
}
} /* namespace prot */
