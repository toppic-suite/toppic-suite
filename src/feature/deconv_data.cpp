#include <cmath>

#include "feature/raw_ms_util.hpp" 
#include "feature/deconv_data.hpp" 

namespace prot {

DeconvData::DeconvData(PeakPtrVec &peak_list, double max_mass, int max_charge,
                       double win_size): 
    peak_list_(peak_list),
    max_mass_(max_mass),
    max_charge_(max_charge),
    win_size_(win_size) {
      win_num_ = (int) std::ceil(RawMsUtil::findMaxPos(peak_list_) / win_size_) + 2;
      initWinId();
      initWinBgnEnd();
    }

// initialize the window id 
void DeconvData::initWinId() {
  for (size_t i = 0; i < peak_list_.size(); i++) {
    win_ids_.push_back(std::floor(peak_list_[i]->getPosition() / win_size_));
  }
}

// initialize the first and last peak for each window 
void DeconvData::initWinBgnEnd() {
  int peak_num = peak_list_.size();
  win_bgn_peaks_.resize(win_num_, peak_num);
  win_end_peaks_.resize(win_num_, -1);
  win_bgn_peaks_[0] = 0;
  win_end_peaks_[win_num_ - 1] = peak_num - 1;
  for (int i = 0; i < peak_num; i++) {
    int id = win_ids_[i];
    if (i < win_bgn_peaks_[id]) {
      win_bgn_peaks_[id] = i;
    }
    if (i > win_end_peaks_[id]) {
      win_end_peaks_[id] = i;
    }
  }

  int bgn = peak_num;
  for (int i = win_num_ - 1; i >= 0; i--) {
    if (win_bgn_peaks_[i] == peak_num) {
      win_bgn_peaks_[i] = bgn;
    } else {
      bgn = win_bgn_peaks_[i];
    }
  }
  int end = -1;
  for (int i = 0; i < win_num_; i++) {
    if (win_end_peaks_[i] == -1) {
      win_end_peaks_[i] = end;
    } else {
      end = win_end_peaks_[i];
    }
  }
}

}
