//
// Created by abbash on 8/24/22.
//

#ifndef TOPPIC_XIC_HPP
#define TOPPIC_XIC_HPP

#include <vector>
#include <numeric>

namespace toppic {
    class Xic {
    public:
        Xic(){ start_spec_id_ = -1; base_spec_id_ = -1; }
        Xic(int start_spec_id, int base_spec_id, std::vector<double>  inte_list);
        Xic(const Xic & x);

        void moving_avg(int n);
        void refine_boundary() { moving_avg(2); }

        int getStartSpecId() { return start_spec_id_; }
        void setStartSpecId(int startSpecId) { start_spec_id_ = startSpecId; }

        int getBaseSpecId() { return base_spec_id_; }
        void setBaseSpecId(int baseSpecId) { base_spec_id_ = baseSpecId; }

        std::vector<double> getInteList() { return inte_list_; }
        void setInteList(std::vector<double> inteList) { inte_list_ = inteList; }

        double get_inte_list_sum(){ return std::accumulate(inte_list_.begin(), inte_list_.end(), 0.0); }

        std::vector<double> getSmoothedInteList() { return smoothed_inte_list_; }
        void setSmoothedInteList(std::vector<double> smoothedInteList) { smoothed_inte_list_ = smoothedInteList; }

        bool isEmpty(){
          if (start_spec_id_ == -1 && base_spec_id_ == -1 && inte_list_.size() == 0 && smoothed_inte_list_.size() == 0)
            return true;
          return false;
        }

    private:
        int start_spec_id_;
        int base_spec_id_;
        std::vector<double> inte_list_;
        std::vector<double> smoothed_inte_list_;
    };
}


#endif //TOPPIC_XIC_HPP
