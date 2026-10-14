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

#include "filter/tag/tag_filter_merge.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "common/base/proteoform_type.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "ms/spec/msalign_util.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/simple_prsm_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace toppic {

namespace tag_filter_merge {

namespace {

const double kPrecMassTole = 1e-4;
const double kShiftTole = 0.01;

// Reads a simple PrSM file one spectrum at a time (the candidates are in
// spectrum order).
class GroupReader {
 public:
  explicit GroupReader(const std::string& file_name) : reader_(file_name) {
    next_ptr_ = reader_.readOnePrsm();
  }

  bool hasNext() const { return next_ptr_ != nullptr; }

  int nextSpecId() const { return next_ptr_->getSpectrumId(); }

  SimplePrsmPtrVec readGroup() {
    SimplePrsmPtrVec group;
    int spec_id = next_ptr_->getSpectrumId();
    while (next_ptr_ != nullptr && next_ptr_->getSpectrumId() == spec_id) {
      group.push_back(next_ptr_);
      next_ptr_ = reader_.readOnePrsm();
    }
    return group;
  }

  void close() { reader_.close(); }

 private:
  SimplePrsmReader reader_;
  SimplePrsmPtr next_ptr_;
};

bool sameShift(double a, double b) { return std::abs(a - b) <= kShiftTole; }

bool sameShiftList(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (size_t i = 0; i < a.size(); i++) {
    if (!sameShift(a[i], b[i])) {
      return false;
    }
  }
  return true;
}

bool sameProtein(const SimplePrsmPtr& a, const SimplePrsmPtr& b) {
  return a->getSeqName() == b->getSeqName() &&
         std::abs(a->getPrecMass() - b->getPrecMass()) <= kPrecMassTole;
}

// adds the shifts of src not within kShiftTole of one of dst's
void addShifts(const std::vector<double>& src, std::vector<double>* dst) {
  for (size_t i = 0; i < src.size(); i++) {
    bool found = false;
    for (size_t j = 0; j < dst->size() && !found; j++) {
      found = sameShift(src[i], (*dst)[j]);
    }
    if (!found) {
      dst->push_back(src[i]);
    }
  }
}

struct Counts {
  size_t filter_num = 0;
  size_t tag_num = 0;
  size_t added_num = 0;
  size_t joined_num = 0;
};

// the candidates of one spectrum: the filter's, then the new tag ones
void combine(const SimplePrsmPtrVec& filter_ptrs, const SimplePrsmPtrVec& tag_ptrs,
             TagFilterMng::Stage stage, SimplePrsmPtrVec* out, Counts* counts) {
  out->insert(out->end(), filter_ptrs.begin(), filter_ptrs.end());
  counts->filter_num += filter_ptrs.size();
  counts->tag_num += tag_ptrs.size();
  for (size_t i = 0; i < tag_ptrs.size(); i++) {
    const SimplePrsmPtr& tag_ptr = tag_ptrs[i];
    bool listed = false;
    for (size_t j = 0; j < filter_ptrs.size() && !listed; j++) {
      const SimplePrsmPtr& filter_ptr = filter_ptrs[j];
      if (!sameProtein(tag_ptr, filter_ptr)) {
        continue;
      }
      switch (stage) {
        case TagFilterMng::Stage::kOne:
          addShifts(tag_ptr->getNTermShifts(), &filter_ptr->getNTermShifts());
          addShifts(tag_ptr->getCTermShifts(), &filter_ptr->getCTermShifts());
          listed = true;
          counts->joined_num++;
          break;
        case TagFilterMng::Stage::kVar:
          listed = tag_ptr->getSeqStartPos() == filter_ptr->getSeqStartPos() &&
                   tag_ptr->getSeqEndPos() == filter_ptr->getSeqEndPos() &&
                   sameShiftList(tag_ptr->getNTermShifts(),
                                 filter_ptr->getNTermShifts()) &&
                   sameShiftList(tag_ptr->getTruncNTermModShifts(),
                                 filter_ptr->getTruncNTermModShifts());
          break;
        case TagFilterMng::Stage::kMulti:
        case TagFilterMng::Stage::kTopmg:
        default:
          listed = true;
          break;
      }
    }
    if (!listed) {
      out->push_back(tag_ptr);
      counts->added_num++;
    }
  }
}

void mergeFiles(const std::string& filter_file, const std::string& tag_file,
                const std::string& out_file,
                const std::unordered_map<int, int>& spec_pos,
                TagFilterMng::Stage stage, Counts* counts) {
  GroupReader filter_reader(filter_file);
  GroupReader tag_reader(tag_file);
  SimplePrsmXmlWriter writer(out_file);
  auto posOf = [&](int spec_id) {
    auto it = spec_pos.find(spec_id);
    if (it == spec_pos.end()) {
      LOG_ERROR("Spectrum id " << spec_id << " is not in the spectrum file!");
      exit(EXIT_FAILURE);
    }
    return it->second;
  };
  while (filter_reader.hasNext() || tag_reader.hasNext()) {
    SimplePrsmPtrVec filter_ptrs;
    SimplePrsmPtrVec tag_ptrs;
    if (!tag_reader.hasNext()) {
      filter_ptrs = filter_reader.readGroup();
    } else if (!filter_reader.hasNext()) {
      tag_ptrs = tag_reader.readGroup();
    } else {
      int filter_pos = posOf(filter_reader.nextSpecId());
      int tag_pos = posOf(tag_reader.nextSpecId());
      if (filter_pos <= tag_pos) {
        filter_ptrs = filter_reader.readGroup();
      }
      if (tag_pos <= filter_pos) {
        tag_ptrs = tag_reader.readGroup();
      }
    }
    SimplePrsmPtrVec out;
    combine(filter_ptrs, tag_ptrs, stage, &out, counts);
    writer.write(out);
  }
  filter_reader.close();
  tag_reader.close();
  writer.close();
}

}  // namespace

void process(const std::string& sp_file_name, const std::string& filter_ext,
             const std::string& tag_ext, const std::string& out_ext,
             TagFilterMng::Stage stage) {
  std::string base_name = file_util::basename(sp_file_name);
  std::vector<int> spec_ids = msalign_util::readSpIds(sp_file_name);
  std::unordered_map<int, int> spec_pos;
  for (size_t i = 0; i < spec_ids.size(); i++) {
    spec_pos.emplace(spec_ids[i], static_cast<int>(i));
  }
  std::vector<std::string> suffixes;
  if (stage == TagFilterMng::Stage::kVar || stage == TagFilterMng::Stage::kOne) {
    const std::vector<ProteoformTypePtr> type_ptrs = {
        ProteoformType::COMPLETE, ProteoformType::PREFIX,
        ProteoformType::SUFFIX, ProteoformType::INTERNAL};
    for (size_t t = 0; t < type_ptrs.size(); t++) {
      suffixes.push_back("_" + type_ptrs[t]->getName());
    }
  } else {
    suffixes.push_back("");
  }
  Counts counts;
  for (size_t i = 0; i < suffixes.size(); i++) {
    mergeFiles(base_name + "." + filter_ext + suffixes[i],
               base_name + "." + tag_ext + suffixes[i],
               base_name + "." + out_ext + suffixes[i], spec_pos, stage,
               &counts);
  }
  std::cout << "Tag filtering (" << TagFilterMng::getStageName(stage)
            << ") - lists joined: " << counts.filter_num
            << " filter candidates, " << counts.tag_num
            << " tag candidates, " << counts.added_num << " added";
  if (stage == TagFilterMng::Stage::kOne) {
    std::cout << ", " << counts.joined_num
              << " joined to listed candidates";
  }
  std::cout << "." << std::endl;
}

}  // namespace tag_filter_merge

}  // namespace toppic
