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

#include "prsm/simple_prsm_mem_merge.hpp"

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "common/base/proteoform_type.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace toppic {

namespace {

// The text record of a candidate as SimplePrsmReader reads it back from a
// block file: the serialized element split into trimmed non-empty lines. The
// merge wrote those lines back, so writing them keeps the file byte for byte
// the same; only the candidates that survive are serialized.
SimplePrsmStrPtr toPrsmStr(const SimplePrsmPtr& prsm_ptr) {
  XmlDOMDocument doc("simple_prsm_list");
  XmlDOMElement root = doc.getDocumentElement();
  XmlDOMElement element = prsm_ptr->toXml(&doc, root);
  std::string str = xml_dom_util::writeToString(element);
  std::vector<std::string> lines;
  std::string line;
  std::istringstream stream(str);
  while (std::getline(stream, line)) {
    str_util::trim(line);
    if (line != "") {
      lines.push_back(line);
    }
  }
  return std::make_shared<SimplePrsmStr>(lines);
}

const int mutex_num = 64;

}  // namespace

// Total order of the candidates of one spectrum: the order they are written
// in (decreasing score, then sequence name, as SimplePrsmStr sorts them) with
// the block, the call inside the block and the position in the call as the
// tie breaker. The tie breaker makes both the trimming and the written order
// independent of the order in which the blocks finish, which the score and
// the sequence name alone do not: a spectrum often has several candidates of
// the same protein and score that differ only in the precursor variant.
bool SimplePrsmMemMerge::cmpEntry(const Entry& a, const Entry& b) {
  double a_score = a.prsm_ptr_->getScore();
  double b_score = b.prsm_ptr_->getScore();
  if (a_score != b_score) {
    return a_score > b_score;
  }
  std::string a_name = a.prsm_ptr_->getSeqName();
  std::string b_name = b.prsm_ptr_->getSeqName();
  if (a_name != b_name) {
    return a_name < b_name;
  }
  if (a.block_idx_ != b.block_idx_) {
    return a.block_idx_ < b.block_idx_;
  }
  if (a.call_idx_ != b.call_idx_) {
    return a.call_idx_ < b.call_idx_;
  }
  return a.pos_idx_ < b.pos_idx_;
}

SimplePrsmMemMerge::SimplePrsmMemMerge(const std::string& sp_file_name,
                                       const std::string& output_file_ext,
                                       const std::vector<int>& top_nums,
                                       const std::vector<int>& spec_ids,
                                       int group_num)
    : sp_file_name_(sp_file_name),
      output_file_ext_(output_file_ext),
      top_nums_(top_nums),
      group_num_(group_num) {
  for (size_t i = 0; i < spec_ids.size(); i++) {
    spec_pos_.emplace(spec_ids[i], static_cast<int>(i));
  }
  results_.resize(top_nums.size() * group_num);
  for (size_t i = 0; i < results_.size(); i++) {
    results_[i].resize(spec_ids.size());
  }
  for (int i = 0; i < mutex_num; i++) {
    mutexes_.push_back(std::unique_ptr<std::mutex>(new std::mutex()));
  }
}

void SimplePrsmMemMerge::add(int type_idx, int block_idx, int call_idx,
                             const SimplePrsmPtrVec& prsm_ptrs, int group) {
  if (prsm_ptrs.size() == 0) {
    return;
  }
  int spec_id = prsm_ptrs[0]->getSpectrumId();
  auto it = spec_pos_.find(spec_id);
  if (it == spec_pos_.end()) {
    LOG_ERROR("Spectrum id " << spec_id << " is not in the spectrum file "
                             << sp_file_name_ << "!");
    exit(EXIT_FAILURE);
  }
  int spec_pos = it->second;
  size_t top_num = static_cast<size_t>(top_nums_[type_idx]);
  std::lock_guard<std::mutex> lock(*mutexes_[spec_pos % mutex_num]);
  std::vector<Entry>& results =
      results_[group * top_nums_.size() + type_idx][spec_pos];
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    Entry entry;
    entry.prsm_ptr_ = prsm_ptrs[i];
    entry.block_idx_ = block_idx;
    entry.call_idx_ = call_idx;
    entry.pos_idx_ = static_cast<int>(i);
    results.push_back(entry);
  }
  if (results.size() > 2 * top_num) {
    // keep the best top_num candidates of the total order, the ones that are
    // written; which of them are kept does not depend on the order the blocks
    // reached this point in
    std::nth_element(results.begin(), results.begin() + top_num, results.end(),
                     cmpEntry);
    results.resize(top_num);
  }
}

void SimplePrsmMemMerge::write() {
  std::string base_name = file_util::basename(sp_file_name_);
  const char* type_names[4] = {
      ProteoformType::COMPLETE->getName().c_str(),
      ProteoformType::PREFIX->getName().c_str(),
      ProteoformType::SUFFIX->getName().c_str(),
      ProteoformType::INTERNAL->getName().c_str()};
  for (size_t type_idx = 0; type_idx < top_nums_.size(); type_idx++) {
    SimplePrsmXmlWriter writer(base_name + "." + output_file_ext_ + "_" +
                               type_names[type_idx]);
    size_t top_num = static_cast<size_t>(top_nums_[type_idx]);
    size_t spec_num = results_[type_idx].size();
    for (size_t spec_pos = 0; spec_pos < spec_num; spec_pos++) {
      // the best top_num of each group, the first group's first; a later
      // group's candidate whose protein and range are listed already is not
      // listed again, its truncation N-terminal modification shifts are
      // added to the listed one's
      std::vector<std::string> keys;
      SimplePrsmPtrVec selected;
      for (int group = 0; group < group_num_; group++) {
        std::vector<Entry>& results =
            results_[group * top_nums_.size() + type_idx][spec_pos];
        if (results.size() == 0) {
          continue;
        }
        std::sort(results.begin(), results.end(), cmpEntry);
        for (size_t i = 0; i < results.size() && i < top_num; i++) {
          const SimplePrsmPtr& prsm_ptr = results[i].prsm_ptr_;
          if (group_num_ > 1) {
            std::string key = prsm_ptr->getSeqName() + " " +
                              std::to_string(prsm_ptr->getSeqStartPos()) +
                              " " + std::to_string(prsm_ptr->getSeqEndPos());
            auto it = std::find(keys.begin(), keys.end(), key);
            if (it != keys.end()) {
              const SimplePrsmPtr& listed_ptr = selected[it - keys.begin()];
              const std::vector<double>& shifts =
                  prsm_ptr->getTruncNTermModShifts();
              for (size_t j = 0; j < shifts.size(); j++) {
                listed_ptr->addTruncNTermModShift(shifts[j]);
              }
              continue;
            }
            keys.push_back(key);
          }
          selected.push_back(prsm_ptr);
        }
        std::vector<Entry>().swap(results);
      }
      for (size_t i = 0; i < selected.size(); i++) {
        writer.write(toPrsmStr(selected[i]));
      }
    }
    writer.close();
  }
}

}  // namespace toppic
