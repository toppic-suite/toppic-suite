// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "filter/index/joint_index.hpp"

#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <thread>

#include "common/thread/simple_thread_pool.hpp"
#include "common/util/logger.hpp"
#include "filter/index/one_ptm_index.hpp"
#include "filter/index/zero_ptm_index.hpp"
#include "filter/massmatch/mass_match_factory.hpp"
#include "seq/db_block.hpp"
#include "seq/proteoform_factory.hpp"
#include "seq/proteoform_util.hpp"

namespace toppic {

namespace joint_index {

// The tasks write their files without the serialization mutex the stage
// processes hold: each task serializes its own MassMatch into its own
// archive, and the writers of an index build (one 60-400 MB file every few
// seconds per thread) would otherwise queue behind one another for a
// growing share of the run once each task writes two files.
namespace {
// the kinds of index, in the order of the managers' file name vectors
const int kKindNum = 4;
const int kDiagKind = 1;  // the PRM diagonal index

// A block's proteoforms and the shift lists the four kinds are built from,
// read once per block and shared by its four tasks (read only).
struct BlockInput {
  ProteoformPtrVec proteo_ptrs;
  std::vector<std::vector<double>> shift_2d;
  std::vector<std::vector<double>> n_term_acet_2d;
  std::vector<std::vector<double>> rev_shift_2d;
};

using BlockInputPtr = std::shared_ptr<BlockInput>;

BlockInputPtr readBlock(int block_idx, const PrsmParaPtr& prsm_para_ptr) {
  std::string db_block_file_name =
      prsm_para_ptr->getSearchDbFileNameWithFolder() + "_" +
      std::to_string(block_idx);
  BlockInputPtr input_ptr = std::make_shared<BlockInput>();
  input_ptr->proteo_ptrs = proteoform_factory::readFastaToProteoformPtrVec(
      db_block_file_name, prsm_para_ptr->getFixModPtrVec());
  input_ptr->shift_2d = proteoform_util::getNTermShift2D(
      input_ptr->proteo_ptrs, prsm_para_ptr->getProtModPtrVec());
  input_ptr->n_term_acet_2d = proteoform_util::getNTermAcet2D(
      input_ptr->proteo_ptrs, prsm_para_ptr->getProtModPtrVec());
  input_ptr->rev_shift_2d.assign(input_ptr->proteo_ptrs.size(),
                                 std::vector<double>(1, 0));
  return input_ptr;
}

// The index of one kind at the given column limits, as the zero-shift and
// one-shift stages build it.
MassMatchPtr buildIndex(int kind, const BlockInputPtr& input_ptr,
                        double max_proteoform_mass, double scale) {
  switch (kind) {
    case 0:
      return mass_match_factory::getPrmTermMassMatchPtr(
          input_ptr->proteo_ptrs, input_ptr->shift_2d, max_proteoform_mass,
          scale);
    case kDiagKind:
      return mass_match_factory::getPrmDiagMassMatchPtr(
          input_ptr->proteo_ptrs, max_proteoform_mass, scale);
    case 2:
      return mass_match_factory::getSrmTermMassMatchPtr(
          input_ptr->proteo_ptrs, input_ptr->rev_shift_2d,
          input_ptr->n_term_acet_2d, max_proteoform_mass, scale);
    default:
      return mass_match_factory::getSrmDiagMassMatchPtr(
          input_ptr->proteo_ptrs, input_ptr->n_term_acet_2d,
          max_proteoform_mass, scale);
  }
}

// The one-shift index files are the zero-shift ones cut at the one-shift
// column limit when the scales agree and that limit is not larger; the file
// name vectors are indexed by kind.
bool canJoinStages(const ZeroPtmFilterMngPtr& zero_mng_ptr,
                   const OnePtmFilterMngPtr& one_mng_ptr) {
  return zero_mng_ptr->filter_scale_ == one_mng_ptr->filter_scale_ &&
         one_mng_ptr->max_proteoform_mass_ <=
             zero_mng_ptr->max_proteoform_mass_ &&
         zero_mng_ptr->getIndexFilePara() == one_mng_ptr->getIndexFilePara() &&
         static_cast<int>(zero_mng_ptr->zero_ptm_file_vec_.size()) ==
             kKindNum &&
         static_cast<int>(one_mng_ptr->one_ptm_file_vec_.size()) == kKindNum;
}

// Writes the zero-shift and one-shift files of one kind for one block.
std::function<void()> geneIndexTask(const BlockInputPtr& input_ptr,
                                    int block_idx, int kind,
                                    const ZeroPtmFilterMngPtr& zero_mng_ptr,
                                    const OnePtmFilterMngPtr& one_mng_ptr) {
  return [input_ptr, block_idx, kind, zero_mng_ptr, one_mng_ptr]() {
    std::string block_str = std::to_string(block_idx);
    std::string dir_name =
        zero_mng_ptr->prsm_para_ptr_->getOriDbName() + "_idx";
    std::string zero_file = zero_mng_ptr->zero_ptm_file_vec_[kind] +
                            zero_mng_ptr->getIndexFilePara() + block_str;
    std::string one_file = one_mng_ptr->one_ptm_file_vec_[kind] +
                           one_mng_ptr->getIndexFilePara() + block_str;

    MassMatchPtr index_ptr =
        buildIndex(kind, input_ptr, zero_mng_ptr->max_proteoform_mass_,
                   zero_mng_ptr->filter_scale_);
    index_ptr->serializeMassMatch(zero_file, dir_name);
    if (index_ptr->canTruncateColumns(one_mng_ptr->max_proteoform_mass_)) {
      index_ptr->truncateColumns(one_mng_ptr->max_proteoform_mass_);
    } else {
      // never for prefix or suffix residue masses, which are nondecreasing
      LOG_WARN("The masses of block " << block_idx
               << " are not nondecreasing; the one-shift index file "
               << one_file << " is built separately.");
      index_ptr = nullptr;
      index_ptr = buildIndex(kind, input_ptr, one_mng_ptr->max_proteoform_mass_,
                             one_mng_ptr->filter_scale_);
    }
    index_ptr->serializeMassMatch(one_file, dir_name);
    index_ptr = nullptr;
  };
}
}  // namespace

void process(const ZeroPtmFilterMngPtr& zero_mng_ptr,
             const OnePtmFilterMngPtr& one_mng_ptr) {
  if (!canJoinStages(zero_mng_ptr, one_mng_ptr)) {
    zero_ptm_index::process(zero_mng_ptr);
    one_ptm_index::process(one_mng_ptr);
    return;
  }

  PrsmParaPtr prsm_para_ptr = zero_mng_ptr->prsm_para_ptr_;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  DbBlockPtrVec db_block_ptr_vec = DbBlock::readDbBlockIndex(db_file_name);
  int block_num = db_block_ptr_vec.size();

  std::cout << "Generating index files --- started" << std::endl;
  SimpleThreadPoolPtr pool_ptr =
      std::make_shared<SimpleThreadPool>(zero_mng_ptr->thread_num_);
  for (int i = 0; i < block_num; i++) {
    std::cout << "Index files - processing " << (i + 1) << " of " << block_num
              << " files." << std::endl;
    int block_idx = db_block_ptr_vec[i]->getBlockIdx();
    BlockInputPtr input_ptr = readBlock(block_idx, prsm_para_ptr);
    for (int kind = 0; kind < kKindNum; kind++) {
      while (pool_ptr->getQueueSize() > 0 ||
             pool_ptr->getIdleThreadNum() == 0) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
      pool_ptr->enqueue(
          geneIndexTask(input_ptr, block_idx, kind, zero_mng_ptr, one_mng_ptr));
    }
  }
  pool_ptr->shutDown();
  std::cout << "Generating index files --- finished" << std::endl;
}

}  // namespace joint_index

}  // namespace toppic
