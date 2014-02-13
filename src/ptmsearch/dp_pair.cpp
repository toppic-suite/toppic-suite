/*
 * dp_pair.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: xunlikun
 */

#include "ptmsearch/dp_pair.hpp"

namespace prot {
DPPair::DPPair(int x,int y,double pair_score,double diff,
        int order,int n_shift,DiagonalHeaderPtr header):Pair(x,y){
    diff_ = diff;
    pair_score_ = pair_score;
    order_= order;
    header_ = header;
    for(int i=0;i<n_shift+1;i++){
        prevs_.push_back(nullptr);
        scores_.push_back(-std::numeric_limits<double>::max());
        types_.push_back(PATH_TYPE_NULL);
    }
}
void DPPair::updateTable(int s,double score,int path_type,DPPairPtr prev_pair){
    scores_[s] = score;
    types_[s] = path_type;
    prevs_[s] = prev_pair;
}
} /* namespace prot */
