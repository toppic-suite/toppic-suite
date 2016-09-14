// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <algorithm>

#include "base/logger.hpp"
#include "zeroptmfilter/filter_protein.hpp"

namespace prot {

inline bool cmpScore(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    return a.second > b.second;
}

inline void sortResult(std::vector<std::pair<int,int>> &results) {
  std::pair<int,int> tmp;
  for (size_t i = 0; i < results.size(); i++) {
    for (size_t j = i+1; j < results.size(); j++) {
      if (results[i].second < results[j].second) {
        tmp = results[i];
        results[i] = results[j];
        results[j] = tmp;
      }   
    }   
  }
}


FilterProteinPtrVec FilterProtein::geneResults(std::vector<std::pair<int,int>> &single_type_results, 
                                               int threshold, int single_type_num) {

  FilterProteinPtrVec prot_results;
  //std::sort(single_type_results.begin(), single_type_results.end(), cmpScore);
  sortResult(single_type_results);
  int output_num =0;
  for(int i=0;i< single_type_num;i++){
    if (i >= (int)single_type_results.size()) {
      break;
    }
    //LOG_DEBUG("rank " << i  <<  " score " << single_type_results[i].second);
    if(single_type_results[i].second >= threshold){
      output_num++;
    }
    else {
      break;
    }
  }
  for(int i=0;i<output_num;i++){
    int prot_id = single_type_results[i].first;
    int score = single_type_results[i].second;
    FilterProteinPtr prot_ptr = FilterProteinPtr(new FilterProtein(prot_id, score));
    prot_results.push_back(prot_ptr);
  }
  return prot_results;
}

} /* namespace prot */

