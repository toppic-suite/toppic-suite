//
// Created by abbash on 8/24/22.
//

#ifndef TOPPIC_EXP_ENVELOPE_HPP
#define TOPPIC_EXP_ENVELOPE_HPP

#include <vector>
#include "topfd/feature_detect/envelope/simple_peak.hpp"
#include "topfd/feature_detect/spectrum/exp_peak.hpp"

namespace toppic {
    class ExpEnvelope {
    public:
        ExpEnvelope(){ spec_id_ = -1; }

        ExpEnvelope(int spec_id, std::vector<ExpPeak> peak_list){
          spec_id_ = spec_id;
          peak_list_ = peak_list;
        }

        int get_match_peak_num();
        std::vector<double> get_inte_list();
        std::vector<double> get_non_empty_pos_list();
        void get_min_max_pos(double* min_pos, double* max_pos);
        int get_peak_num() { return peak_list_.size(); }
        ExpPeak get_peak(int idx) { return peak_list_[idx]; }

        const int getSpecId() const { return spec_id_; }
        void setSpecId(int spec_id) { spec_id_ = spec_id; }

        std::vector<ExpPeak> getExpEnvList() { return peak_list_; }
        void setExpEnvList(std::vector<ExpPeak> peak_list) { peak_list_ = peak_list; }

        bool isEmpty(){
          if (spec_id_ == -1 && peak_list_.empty())
            return true;
          return false;
        }

    private:
        int spec_id_;
        std::vector<ExpPeak> peak_list_;

    };
}


#endif //TOPPIC_EXP_ENVELOPE_HPP
