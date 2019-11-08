//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include <cmath>

#include "common/util/logger.hpp"
#include "search/ptmsearch/comp_shift_low_mem.hpp"

namespace toppic {

std::vector<std::vector<int>> CompShiftLowMem::findBestShift(
    const std::vector<int> &a,const std::vector<int> &b){
  return findBestShift(a,b,1,0);
}

std::vector<double> CompShiftLowMem::findBestShift(const std::vector<std::pair<int,int>> &a_e,
                                                   const std::vector<int> &b,
                                                   int total,int min_gap,
                                                   double scale){
  std::vector<int> a;
  std::vector<int> e;
  for (size_t i = 0; i < a_e.size(); i++) {
    a.push_back(a_e[i].first);
    e.push_back(a_e[i].second);
  }
  std::vector<std::vector<int>> list = findBestShift(a,e, b,total,min_gap);
  std::vector<double> result;
  for(size_t i = 0;i<list.size();i++){
    result.push_back(list[i][0]/scale);
  }
  return result;
}

/*
std::vector<double> CompShiftLowMem::findBestShift(const std::vector<int> &a,
                                                   const std::vector<int> &errors,
                                                   const std::vector<int> &b,
                                                   int total,int min_gap,
                                                   double scale){
  std::vector<std::vector<int>> list = findBestShift(a,errors,b,total,min_gap);
  std::vector<double> result;
  for(size_t i = 0;i<list.size();i++){
    result.push_back(list[i][0]/scale);
  }
  return result;
}
*/

inline void CompShiftLowMem::resetNumbers(const std::vector<int> &a,
                                          const std::vector<int> &errors,
                                          const std::vector<int> &b) {
  int a_length = a.size();
  int b_length = b.size();
  int delta = a[0]-b[b_length-1];
  int minD = delta -1;
  for(int i=0;i<a_length;i++){
    int a_value = a[i];
    for(int j =0;j <b_length;j++){
      int d = a_value - b[j] - minD;
      int start = 0;
      if(d > errors[i]){
        start = d- errors[i];
      }
      for(int k = start;k<=d+errors[i];k++){
        num_[k]=0;
      }
    }
  }
}

inline void CompShiftLowMem::resetNumbers(const std::vector<int> &a,
                                          const std::vector<int> &b) {
  int a_length = a.size();
  int b_length = b.size();
  int delta = a[0]-b[b_length-1];
  int minD = delta -1;
  for(int i=0;i<a_length;i++){
    int a_value = a[i];
    for(int j =0;j <b_length;j++){
      int d = a_value - b[j]- minD;
      num_[d-1] =0;
      num_[d] =0;
      num_[d+1] = 0;
    }
  }
}

inline void CompShiftLowMem::resetNumbers() {
  std::fill(num_.begin(), num_.end(), 0);
}

std::vector<std::vector<int>> CompShiftLowMem::findBestShift(
    const std::vector<int> &a, const std::vector<int> &errors,
    const std::vector<int> &b, int total,int minimum_gap){
  int a_length = a.size();
  int b_length = b.size();
  std::vector<std::vector<int>> ans;
  if(b_length == 0){
    return ans;
  }
  int delta = a[0]-b[b_length-1];
  int minD = delta-1;
  int maxD = a[a_length-1]-b[0]+1+errors[a_length-1];
  if(maxD-minD + 1 >= (int)num_.size()){
    size_t required_len = maxD - minD + 1;
    for (size_t i = num_.size(); i < required_len; i++) {
      num_.push_back(0);
    }
  }

  resetNumbers(a, errors, b);
  // resetNumbers();

  int current_minimum =1;

  for(int i=0;i<a_length;i++){
    int a_value = a[i];
    for(int j =0;j <b_length;j++){
      int d = a_value - b[j] - minD;
      int start = 0;
      if(d >errors[i]){
        start = d- errors[i];
      }
      for(int k = start;k<=d+errors[i];k++){
        num_[k]++;
        current_minimum = checkD(ans,k,current_minimum,total,minimum_gap);
      }
    }
  }

  for(size_t i=0;i<ans.size();i++){
    ans[i][0] += minD;
  }
  return ans;
}

inline std::vector<std::vector<int>> CompShiftLowMem::findBestShift(
    const std::vector<int> &a, const std::vector<int> &b,int total,int min_gap) {
  const int a_length = a.size();
  const int b_length = b.size();
  std::vector<std::vector<int>> ans;
  const int delta = a[0]-b[b_length-1];
  const int minD = delta-1;
  const int maxD = a[a_length-1]-b[0]+1;
  if(maxD-minD + 1 >= (int)num_.size()) {
    int required_len = maxD - minD + 1;
    for (int i = num_.size(); i < required_len; i++) {
      num_.push_back(0);
    }
  }

  resetNumbers(a, b);
  //resetNumbers();

  int cur_min =1;

  for(int i=0;i<a_length;i++){
    int a_value = a[i];
    for(int j =0;j <b_length;j++){
      int d = a_value - b[j]- minD;
      num_[d-1] ++;
      cur_min=checkD(ans,d-1,cur_min,total,min_gap);
      num_[d] ++;
      cur_min=checkD(ans,d,cur_min,total,min_gap);
      num_[d+1] ++;
      cur_min=checkD(ans,d+1,cur_min,total,min_gap);
    }
  }

  for(size_t i=0;i<ans.size();i++){
    ans[i][0] +=minD;
  }
  return ans;
}

inline int CompShiftLowMem::checkD(std::vector<std::vector<int>> &ans,int d,
                                   int cur_min,int total,int min_gap){
  short new_value = num_[d];
  if(new_value < cur_min){
    return cur_min;
  }
  for(size_t i =0;i<ans.size();i++){
    std::vector<int>  cur = ans[i];
    if(std::abs(cur[0]-d)<= min_gap){
      if(cur[1]<new_value){
        ans.erase(ans.begin()+i);
        i--;
      }
      else {
        return cur_min;
      }
    }
  }

  int insert_pos = ans.size()-1;

  while(insert_pos >= 0 && ans[insert_pos][1] < new_value){
    insert_pos--;
  }
  std::vector<int> insert_temp = {d,new_value};
  ans.insert(ans.begin()+insert_pos+1, insert_temp);
  if((int)ans.size()>total){
    ans.pop_back();
  }

  if ((int)ans.size() == total) {
    return ans[ans.size()-1][1]+1;
  }
  else {
    return 1;
  }
}

} /* namespace toppic */
