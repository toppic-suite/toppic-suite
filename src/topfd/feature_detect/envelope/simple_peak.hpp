//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_SIMPLE_PEAK_HPP
#define TOPPIC_SIMPLE_PEAK_HPP

#include "topfd/feature_detect/spectrum/peak_matrix.hpp"

namespace toppic {
    class SimplePeak {
      public:
        SimplePeak() { pos_ = -1; inte_ = -1; }
        SimplePeak(double pos, double inte);
        SimplePeak(SimplePeak const &p);

        bool isEmpty(){
          if (pos_ == -1 and inte_ == -1) return true;
          return false;
        }
        double getPos() const { return pos_; }
        void setPos(double pos) { pos_ = pos; }

        double getInte() const { return inte_; }
        void setInte(double inte) { inte_ = inte; }

        int getStartIdx() const { return start_idx_; }
        void setStartIdx(int startIdx) { start_idx_ = startIdx; }

        int getEndIdx() const { return end_idx_; }
        void setEndIdx(int endIdx) { end_idx_ = endIdx; }

        static bool cmpInteDec(const SimplePeak &a, const SimplePeak &b) { return a.getInte() > b.getInte();}
        static bool cmpPos(const SimplePeak &a, const SimplePeak &b) { return a.getPos() < b.getPos();}
    private:
        double pos_;
        double inte_;
        int start_idx_ = -1;
        int end_idx_ = -1;

    };
}

#endif //TOPPIC_SIMPLE_PEAK_HPP
