#ifndef PROT_BASE_DB_BLOCK_HPP_
#define PROT_BASE_DB_BLOCK_HPP_

#include <memory>
#include <vector>

namespace prot {

class DbBlock;
typedef std::shared_ptr<DbBlock> DbBlockPtr;
typedef std::vector<DbBlockPtr> DbBlockPtrVec;

class DbBlock {
 public:
  DbBlock(int block_index, int seq_index);

  int getBlockIdx() {return block_index_;}
  int getSeqIdx() {return seq_index_;}

  static DbBlockPtrVec readDbBlockIndex(const std::string &index_file_name);

 private:
  int block_index_;
  int seq_index_;
};

}
#endif
