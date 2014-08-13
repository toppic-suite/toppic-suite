#include "ptmsearch/dp_pair.hpp"

namespace prot {
DPPair::DPPair(int x,int y,double pair_score,double diff,
               int order,int n_shift,DiagonalHeaderPtr header_ptr):Pair(x,y){
  diff_ = diff;
  pair_score_ = pair_score;
  order_= order;
  header_ptr_ = header_ptr;
  for(int i=0;i<n_shift+1;i++){
    prev_pair_ptrs_.push_back(nullptr);
    scores_.push_back(-std::numeric_limits<double>::max());
    types_.push_back(PATH_TYPE_NULL);
  }
}

void DPPair::updateTable(int s,double score,int path_type,DPPairPtr prev_pair_ptr){
    scores_[s] = score;
    types_[s] = path_type;
    prev_pair_ptrs_[s] = prev_pair_ptr;
}

} /* namespace prot */
