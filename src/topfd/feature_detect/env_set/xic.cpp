//
// Created by abbash on 8/24/22.
//

#include "xic.hpp"

toppic::Xic::Xic(int start_spec_id, int base_spec_id, std::vector<double>  inte_list){
  start_spec_id_ = start_spec_id;
  base_spec_id_ = base_spec_id;
  inte_list_ = inte_list;
}

toppic::Xic::Xic(const Xic & x) {
  start_spec_id_ = x.start_spec_id_;
  base_spec_id_ = x.base_spec_id_;
  for (auto inte : x.inte_list_)
    inte_list_.push_back(inte);
  for (auto smoothed_inte: x.smoothed_inte_list_)
    smoothed_inte_list_.push_back(smoothed_inte);
}

void toppic::Xic::moving_avg (int size) {
  std::vector<double> data = inte_list_;
  int padding_len = (size -1) / 2;
  std::vector<double> left_padding (padding_len, data[0]);
  data.insert(data.begin(), left_padding.begin(), left_padding.end());
  std::vector<double> right_padding (padding_len, data[data.size()-1]);
  data.insert(data.end(), right_padding.begin(), right_padding.end());

  double sum = 0.0;
  int cnt = 0;
  for (int i = 0; i < data.size(); i++) {
    sum += data[i];
    cnt++;
    if (cnt >= size) {
      smoothed_inte_list_.push_back((sum / (double) size));
      sum -= data[cnt - size];
    }
  }
}
