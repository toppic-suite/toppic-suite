
#ifndef PROT_TAG_FILTER_PROCESSOR
#define PROT_TAG_FILTER_PROCESSOR

#include "base/db_block.hpp"
#include "tagfilter/tag_filter_mng.hpp"

namespace prot {

class TagFilterProcessor {
 public:
  TagFilterProcessor(TagFilterMngPtr mng_ptr);
  void process();

 private:
  TagFilterMngPtr mng_ptr_;
  void processBlock(DbBlockPtr block_ptr, int total_block_num);

};

typedef std::shared_ptr<TagFilterProcessor> TagFilterProcessorPtr;

}

#endif
