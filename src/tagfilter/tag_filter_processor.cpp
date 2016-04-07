
#include "base/proteoform.hpp"
#include "base/base_data.hpp"
#include "base/db_block.hpp"
#include "tag_filter_processor.hpp"

namespace prot{

TagFilterProcessor::TagFilterProcessor(TagFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
}

void TagFilterProcessor::process(){

}

void TagFilterProcessor::processBlock(DbBlockPtr block_ptr, 
                                      int total_block_num) {

}

}
