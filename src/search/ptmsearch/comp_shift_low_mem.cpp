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

#include "search/ptmsearch/comp_shift_low_mem.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <mutex>

#include "common/util/logger.hpp"

namespace toppic {

namespace {
// How many (spectrum mass, prefix mass) pairs ahead the counting loop of the
// full-range findBestShift prefetches the count bins of a pair's window.
const int prefetch_pair_num = 32;

// The count buffer of the full-range pass is kept for the next candidate
// (it is reused by the candidates of a thread) unless it is larger than
// this, which happens only for a very long protein: one bin per 0.01 Da of
// the range of the differences.
const size_t max_kept_bin_num = 64 * 1024 * 1024;

// The chunked pass counts the differences in chunks of this many bins
// (128 kB of shorts, which stay in a core's L2 cache), and keeps the count
// events of bins that reach chunked_min_count; when the replay of the
// events does not fill the shift list, the pass is repeated with the
// threshold one lower down to chunked_last_min_count, then the full-range
// loop runs (toppic_test experiment 150).
const int chunk_bins = 1 << 16;
const int chunked_min_count = 4;
const int chunked_last_min_count = 2;

// Research dump (toppic_test experiment 150): with TOPPIC_SHIFT_DUMP=<file>
// every call's inputs are appended as int32 records (a_len, b_len, total,
// min_gap, a, errors, b) for an offline benchmark of the finder.
namespace shift_dump {
std::mutex dump_mutex;
const char* dumpFile() {
  static const char* file = std::getenv("TOPPIC_SHIFT_DUMP");
  return file;
}
bool enabled() { return dumpFile() != nullptr; }
void write(const std::vector<int>& a, const std::vector<int>& e,
           const std::vector<int>& b, int total, int min_gap) {
  std::lock_guard<std::mutex> lock(dump_mutex);
  static std::ofstream out(dumpFile(), std::ios::binary);
  int32_t head[4] = {static_cast<int32_t>(a.size()),
                     static_cast<int32_t>(b.size()), total, min_gap};
  out.write(reinterpret_cast<const char*>(head), sizeof(head));
  out.write(reinterpret_cast<const char*>(a.data()), a.size() * sizeof(int));
  out.write(reinterpret_cast<const char*>(e.data()), e.size() * sizeof(int));
  out.write(reinterpret_cast<const char*>(b.data()), b.size() * sizeof(int));
}
}  // namespace shift_dump
}  // namespace

std::vector<std::vector<int>> CompShiftLowMem::findBestShift(
    const std::vector<int>& a, const std::vector<int>& b) {
  return findBestShift(a, b, 1, 0);
}

std::vector<double> CompShiftLowMem::findBestShift(
    const std::vector<std::pair<int, int>>& a_e, const std::vector<int>& b,
    int total, int min_gap, double scale) {
  std::vector<int> a;
  std::vector<int> e;
  for (size_t i = 0; i < a_e.size(); i++) {
    a.push_back(a_e[i].first);
    e.push_back(a_e[i].second);
  }
  std::vector<std::vector<int>> list = findBestShift(a, e, b, total, min_gap);
  std::vector<double> result;
  for (size_t i = 0; i < list.size(); i++) {
    result.push_back(list[i][0] / scale);
  }
  return result;
}

inline void CompShiftLowMem::resetNumbers(const std::vector<int>& a,
                                          const std::vector<int>& b) {
  int a_length = a.size();
  int b_length = b.size();
  int delta = a[0] - b[b_length - 1];
  int minD = delta - 1;
  for (int i = 0; i < a_length; i++) {
    int a_value = a[i];
    for (int j = 0; j < b_length; j++) {
      int d = a_value - b[j] - minD;
      num_[d - 1] = 0;
      num_[d] = 0;
      num_[d + 1] = 0;
    }
  }
}

inline void CompShiftLowMem::resetNumbers() {
  std::fill(num_.begin(), num_.end(), 0);
}

// The shift list of the full-range loop (findBestShiftFull) from the events
// of the bins that reach min_count, counted chunk by chunk. The full-range
// loop runs over the pairs (i, j) in order and increments the window
// [d - e_i, d + e_i] of every pair, d = a_i - b_j - minD, bin by bin; its
// shift list depends on the order in which the bins reach their counts.
// Only the events (bin, count) with a count of at least the smallest count
// of the final list matter for it: an event below it is rejected by
// cur_min (which never falls once the list is full) or sits in the tail of
// the list until it is pushed out, and it neither rejects nor evicts an
// entry of a larger count. So the events of count >= min_count, replayed
// in the loop's order, give the same list when the replay fills it, and
// the caller tries again with a lower threshold when it does not.
//
// The chunks run from the top of the range down: for a spectrum mass i the
// pairs of a chunk have increasing j (d falls with j), so its events of
// all chunks come in the loop's order, and the replay runs the masses in
// order. A window across a chunk's lower edge is incremented lower bins
// first by the loop, so its events of the upper chunk wait in pending_
// until the lower chunk has emitted its part. The counts of a bin are
// complete within its chunk: every pair whose window meets the chunk is
// visited there, in the loop's order.
bool CompShiftLowMem::findBestShiftChunked(const std::vector<int>& a,
                                           const std::vector<int>& errors,
                                           const std::vector<int>& b,
                                           int total, int minimum_gap,
                                           int min_count) {
  const int a_length = a.size();
  const int b_length = b.size();
  const int minD = a[0] - b[b_length - 1] - 1;
  const int maxD = a[a_length - 1] - b[0] + 1 + errors[a_length - 1];
  const int bin_num = maxD - minD + 1;
  if (chunk_counts_.empty()) {
    chunk_counts_.resize(chunk_bins);
  }
  if (static_cast<int>(events_.size()) < a_length) {
    events_.resize(a_length);
    pending_.resize(a_length);
  }
  for (int i = 0; i < a_length; i++) {
    events_[i].clear();
    pending_[i].clear();
  }
  pair_begins_.assign(a_length, 0);
  pair_ends_.assign(a_length, 0);
  for (int lo = ((bin_num - 1) / chunk_bins) * chunk_bins; lo >= 0;
       lo -= chunk_bins) {
    const int hi = std::min(lo + chunk_bins, bin_num);
    std::fill(chunk_counts_.begin(), chunk_counts_.begin() + (hi - lo), 0);
    short* counts = chunk_counts_.data() - lo;
    for (int i = 0; i < a_length; i++) {
      const int base = a[i] - minD;
      const int error = errors[i];
      // the pairs whose window meets [lo, hi): d in [lo - e, hi - 1 + e],
      // b[j] = base - d; the ranges move up as the chunks move down
      const int b_min = base - (hi - 1 + error);
      const int b_max = base - (lo - error);
      int j_begin = pair_begins_[i];
      int j_end = pair_ends_[i];
      while (j_begin < b_length && b[j_begin] < b_min) {
        j_begin++;
      }
      if (j_end < j_begin) {
        j_end = j_begin;
      }
      while (j_end < b_length && b[j_end] <= b_max) {
        j_end++;
      }
      pair_begins_[i] = j_begin;
      pair_ends_[i] = j_end;
      std::vector<Event>& events = events_[i];
      for (int j = j_begin; j < j_end; j++) {
        const int d = base - b[j];
        const int start = d > error ? d - error : 0;
        const int end = d + error;
        const int chunk_start = start > lo ? start : lo;
        const int chunk_end = end < hi - 1 ? end : hi - 1;
        std::vector<Event>& out = start < lo ? pending_[i] : events;
        short max_count = 0;
        for (int k = chunk_start; k <= chunk_end; k++) {
          short count = ++counts[k];
          max_count = count > max_count ? count : max_count;
        }
        if (max_count >= min_count) {
          for (int k = chunk_start; k <= chunk_end; k++) {
            if (counts[k] >= min_count) {
              out.push_back(Event{k, counts[k]});
            }
          }
        }
        if (end > hi - 1) {
          events.insert(events.end(), pending_[i].begin(), pending_[i].end());
          pending_[i].clear();
        }
      }
    }
  }
  initBestShifts(total);
  int current_minimum = 1;
  for (int i = 0; i < a_length; i++) {
    for (const Event& event : events_[i]) {
      current_minimum =
          checkShift(event.k, event.v, current_minimum, total, minimum_gap);
    }
  }
  return static_cast<int>(best_.size()) == total;
}

std::vector<std::vector<int>> CompShiftLowMem::findBestShift(
    const std::vector<int>& a, const std::vector<int>& errors,
    const std::vector<int>& b, int total, int minimum_gap) {
  if (b.empty()) {
    return std::vector<std::vector<int>>();
  }
  if (shift_dump::enabled()) {
    shift_dump::write(a, errors, b, total, minimum_gap);
  }
  // The full-range loop touches one random cache line per pair in a buffer
  // of tens of MB and zeroes the buffer for every candidate; with 16
  // threads it is bound by the memory bandwidth (50 ms per candidate on
  // 480_F3_02 against 5 ms single-threaded). The chunked pass counts in
  // the L2 cache and replays the events of the frequent bins (5 ms per
  // candidate on 16 threads, the same lists; toppic_test experiment 150).
  for (int min_count = chunked_min_count; min_count >= chunked_last_min_count;
       min_count--) {
    if (findBestShiftChunked(a, errors, b, total, minimum_gap, min_count)) {
      return getBestShifts(a[0] - b[b.size() - 1] - 1);
    }
  }
  return findBestShiftFull(a, errors, b, total, minimum_gap);
}

std::vector<std::vector<int>> CompShiftLowMem::findBestShiftFull(
    const std::vector<int>& a, const std::vector<int>& errors,
    const std::vector<int>& b, int total, int minimum_gap) {
  int a_length = a.size();
  int b_length = b.size();
  int delta = a[0] - b[b_length - 1];
  int minD = delta - 1;
  int maxD = a[a_length - 1] - b[0] + 1 + errors[a_length - 1];
  if (maxD - minD + 1 >= (int)num_.size()) {
    // one allocation instead of one push_back per bin; the buffer is tens of
    // MB and is reused by the candidates of a thread
    num_.resize(maxD - minD + 1, 0);
  }

  // One sequential fill of the bins of this call. On 16 threads it is cheaper
  // than zeroing the window of every pair, one random cache line per pair,
  // although it writes more bytes (toppic_test experiment 78).
  std::fill(num_.begin(), num_.begin() + (maxD - minD + 1), 0);
  initBestShifts(total);

  int current_minimum = 1;
  short* num = num_.data();

  // The pairs touch a random cache line each, so the loop is bound by
  // memory latency: the window of the pair prefetch_pair_num pairs ahead is
  // prefetched, and a pair's window is incremented as a whole (vectorised,
  // no branch per bin) before it is scanned for bins that reach the current
  // minimum, which happens for about 1 in 5,000 increments. checkShift reads
  // only the count of the bin it is given and a window holds each bin once,
  // so the bins, the order and the counts checkShift sees are those of an
  // increment-and-check loop over the bins.
  for (int i = 0; i < a_length; i++) {
    int base = a[i] - minD;
    int error = errors[i];
    for (int j = 0; j < b_length; j++) {
      if (j + prefetch_pair_num < b_length) {
        int ahead = base - b[j + prefetch_pair_num] - error;
        __builtin_prefetch(num + (ahead > 0 ? ahead : 0), 1);
      }
      int d = base - b[j];
      int start = d > error ? d - error : 0;
      int end = d + error;
      short max_count = 0;
      for (int k = start; k <= end; k++) {
        short count = ++num[k];
        max_count = count > max_count ? count : max_count;
      }
      if (max_count >= current_minimum) {
        for (int k = start; k <= end; k++) {
          // checkShift returns cur_min unchanged when num[k] < cur_min
          if (num[k] >= current_minimum) {
            current_minimum =
                checkShift(k, num[k], current_minimum, total, minimum_gap);
          }
        }
      }
    }
  }

  std::vector<std::vector<int>> best_shifts = getBestShifts(minD);
  if (num_.size() > max_kept_bin_num) {
    std::vector<short>().swap(num_);
  }
  return best_shifts;
}

inline std::vector<std::vector<int>> CompShiftLowMem::findBestShift(
    const std::vector<int>& a, const std::vector<int>& b, int total,
    int min_gap) {
  const int a_length = a.size();
  const int b_length = b.size();
  const int delta = a[0] - b[b_length - 1];
  const int minD = delta - 1;
  const int maxD = a[a_length - 1] - b[0] + 1;
  if (maxD - minD + 1 >= (int)num_.size()) {
    num_.resize(maxD - minD + 1, 0);
  }

  resetNumbers(a, b);
  initBestShifts(total);

  int cur_min = 1;

  for (int i = 0; i < a_length; i++) {
    int a_value = a[i];
    for (int j = 0; j < b_length; j++) {
      int d = a_value - b[j] - minD;
      num_[d - 1]++;
      cur_min = checkShift(d - 1, num_[d - 1], cur_min, total, min_gap);
      num_[d]++;
      cur_min = checkShift(d, num_[d], cur_min, total, min_gap);
      num_[d + 1]++;
      cur_min = checkShift(d + 1, num_[d + 1], cur_min, total, min_gap);
    }
  }

  return getBestShifts(minD);
}

inline void CompShiftLowMem::initBestShifts(int total) {
  best_.clear();
  best_.reserve(total + 1);
}

// The kept shifts as (difference, count) pairs.
std::vector<std::vector<int>> CompShiftLowMem::getBestShifts(int min_d) {
  std::vector<std::vector<int>> best_shifts;
  best_shifts.reserve(best_.size());
  for (const Shift& shift : best_) {
    best_shifts.push_back({shift.bin + min_d, shift.count});
  }
  return best_shifts;
}

// Offers bin d with its count to the kept shifts: a kept shift within
// min_gap bins with a smaller count is removed, one with the same or a larger
// count rejects d; the list keeps the best total shifts. Returns the smallest
// count that can still enter the list. best_ is a flat, reused list: this is
// called for every count increment that reaches the threshold, so it must not
// allocate (the former list of 2-int vectors allocated for every entry it
// looked at and every entry it inserted).
inline int CompShiftLowMem::checkShift(int d, short new_value, int cur_min,
                                       int total, int min_gap) {
  if (new_value < cur_min) {
    return cur_min;
  }
  for (size_t i = 0; i < best_.size(); i++) {
    if (std::abs(best_[i].bin - d) <= min_gap) {
      if (best_[i].count < new_value) {
        best_.erase(best_.begin() + i);
        i--;
      } else {
        return cur_min;
      }
    }
  }

  int insert_pos = best_.size() - 1;

  while (insert_pos >= 0 && best_[insert_pos].count < new_value) {
    insert_pos--;
  }
  best_.insert(best_.begin() + insert_pos + 1, Shift{d, new_value});
  if ((int)best_.size() > total) {
    best_.pop_back();
  }

  if ((int)best_.size() == total) {
    return best_.back().count + 1;
  } else {
    return 1;
  }
}

} /* namespace toppic */
