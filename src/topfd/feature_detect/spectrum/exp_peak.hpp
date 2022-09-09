//
// Created by abbash on 8/23/22.
//

#ifndef TOPPIC_EXP_PEAK_HPP
#define TOPPIC_EXP_PEAK_HPP

#include "ms/spec/peak.hpp"

namespace toppic {
    class ExpPeak;
    class ExpPeak {
      public:
        ExpPeak(){ peak_id_ = -1; spec_id_ = -1; pos_ = -1; inte_ = -1; ori_inte_ = -1; start_idx_ = -1; end_idx_ = -1; neighbor_ = false; }
        ExpPeak(int peak_id, int spec_id, PeakPtr peak){ peak_id_ = peak_id; spec_id_ = spec_id; pos_ = peak->getPosition(); inte_ = peak->getIntensity(); ori_inte_ = peak->getIntensity(); start_idx_ = -1; end_idx_ = -1; neighbor_ = false; }
        ExpPeak(const ExpPeak &p) { peak_id_ = p.getPeakId(); spec_id_ = p.getSpecId(); pos_ = p.getPos(); inte_ = p.getInte(); ori_inte_ = p.getOriInte(); start_idx_ = p.getStartIdx(); end_idx_ = p.getEndIdx(); neighbor_ = p.getNeighbor(); }

        bool isEmpty(){
          if (peak_id_ == -1 && spec_id_ == -1 && pos_ == -1 && inte_ == -1) return true;
          return false;
        }

        int getPeakId() const { return peak_id_; }
        void setPeakId(int peak_id) { peak_id_ = peak_id; }

        int getSpecId() const { return spec_id_; }
        void setSpecId(int spec_id) { spec_id_ = spec_id; }

        double getPos() const { return pos_; }
        void setPos(double pos) { pos_ = pos; }

        double getInte() const { return inte_; }
        void setInte(double inte) { inte_ = inte; }

        double getOriInte() const { return ori_inte_; }
        void setOriInte(double ori_inte) { ori_inte_ = ori_inte; }

        int getStartIdx() const { return start_idx_; }
        void setStartIdx(int start_idx) { start_idx_ = start_idx; }

        int getEndIdx() const { return end_idx_; }
        void setEndIdx(int end_idx) { end_idx_ = end_idx; }

        bool getNeighbor() const { return neighbor_; }
        void setNeighbor(bool neighbor) { neighbor_ = neighbor; }

        std::string getString() {
          return "Peak ID: " + std::to_string(peak_id_) + " " +
                 "Spec ID: " + std::to_string(spec_id_) + " " +
                 "Pos: " + std::to_string(pos_) + " " +
                 "Inte: " + std::to_string(inte_) + " " +
                 "Neighbor: " + std::to_string(neighbor_) + "\n";
        }

    private:
        int peak_id_;
        int spec_id_;
        double pos_;
        double inte_;
        double ori_inte_;
        int start_idx_;
        int end_idx_;
        bool neighbor_;
    };

}
#endif //TOPPIC_EXP_PEAK_HPP
