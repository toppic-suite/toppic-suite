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

#ifndef TOPPIC_SEARCH_PTMSEARCH_COMP_SHIFT_LOW_MEM_HPP_
#define TOPPIC_SEARCH_PTMSEARCH_COMP_SHIFT_LOW_MEM_HPP_

#include <memory>
#include <vector>

namespace toppic {

// The diagonal finder of the multiple-shift search: the best shifts between
// a spectrum's mass list a (with tolerances) and a protein's prefix mass
// list b, counted in bins of the scaled masses. One instance is reused by
// the candidates of a thread (its buffers are kept).
class CompShiftLowMem {
 public:
  std::vector<std::vector<int>> findBestShift(const std::vector<int>& a,
                                              const std::vector<int>& b);

  std::vector<double> findBestShift(const std::vector<std::pair<int, int>>& a_e,
                                    const std::vector<int>& b, int total,
                                    int minimum_gap, double scale);

  std::vector<std::vector<int>> findBestShift(const std::vector<int>& a,
                                              const std::vector<int>& errors,
                                              const std::vector<int>& b,
                                              int total, int minimum_gap);

 private:
  // A kept shift: its bin (an index of the difference range until the end
  // of the call) and its count.
  struct Shift {
    int bin;
    int count;
  };

  // A count event of the chunked pass: bin k reached count v.
  struct Event {
    int k;
    short v;
  };

  // the count bins of the full-range pass (findBestShiftFull)
  std::vector<short> num_;
  // The best shifts of the current call, count decreasing, at most total + 1
  // entries; reused by the calls of a thread like num_, so that checkShift
  // allocates nothing.
  std::vector<Shift> best_;

  // the chunked pass: the count bins of one chunk, the events of every
  // spectrum mass in the order of the full-range loop, the events of a
  // window across a chunk edge waiting for the lower chunk, the pair ranges
  std::vector<short> chunk_counts_;
  std::vector<std::vector<Event>> events_;
  std::vector<std::vector<Event>> pending_;
  std::vector<int> pair_begins_;
  std::vector<int> pair_ends_;

  std::vector<std::vector<int>> findBestShift(const std::vector<int>& a,
                                              const std::vector<int>& b,
                                              int total, int minimum_gap);
  std::vector<std::vector<int>> findBestShiftFull(
      const std::vector<int>& a, const std::vector<int>& errors,
      const std::vector<int>& b, int total, int minimum_gap);
  bool findBestShiftChunked(const std::vector<int>& a,
                            const std::vector<int>& errors,
                            const std::vector<int>& b, int total,
                            int minimum_gap, int min_count);
  void initBestShifts(int total);
  std::vector<std::vector<int>> getBestShifts(int min_d);
  int checkShift(int d, short count, int cur_min, int total, int min_gap);

  void resetNumbers(const std::vector<int>& a, const std::vector<int>& b);
  void resetNumbers();
};

using CompShiftLowMemPtr = std::shared_ptr<CompShiftLowMem>;

} /* namespace toppic */

#endif /* COMP_SHIFT_LOW_MEM_HPP_ */
