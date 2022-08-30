//
// Created by abbash on 8/24/22.
//

#include "exp_envelope.hpp"
#include <algorithm>

int toppic::ExpEnvelope::get_match_peak_num(){
  int num = 0;
  for (auto p : peak_list_){
    if (!p.isEmpty())
      num = num + 1;
  }
  return num;
}

std::vector<double> toppic::ExpEnvelope::get_inte_list(){
  std::vector<double> inte_list;
  for (auto p : peak_list_) {
    if (!p.isEmpty())
      inte_list.push_back(p.getInte());
    else
      inte_list.push_back(0);
  }
  return inte_list;
}

std::vector<double> toppic::ExpEnvelope::get_non_empty_pos_list(){
  std::vector<double> pos_list;
  for (auto p : peak_list_) {
    if (!p.isEmpty())
      pos_list.push_back(p.getPos());
  }
  return pos_list;
}

void toppic::ExpEnvelope::get_min_max_pos(double* min_pos, double* max_pos){
  std::vector<double> pos_list = get_non_empty_pos_list();
  *min_pos = *std::min_element(pos_list.begin(), pos_list.end());
  *max_pos = *std::max_element(pos_list.begin(), pos_list.end());
}
