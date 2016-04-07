
#include "tag_filter.hpp"

namespace prot{

seq_tag geneSeqTag(const std::string & seq) {


}

TagFilter::TagFilter(const ProteoformPtrVec &proteo_ptrs,
                     TagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    seq_tag_map_[proteo_ptrs[i]->getFastaSeqPtr()->getName()]
        = geneSeqTag(proteo_ptrs[i]->getFastaSeqPtr()->getSeq());
  }
}

SimplePrsmPtrVec TagFilter::getBestMatch(const PrmMsPtrVec &ms_ptr_vec) {

}


}
