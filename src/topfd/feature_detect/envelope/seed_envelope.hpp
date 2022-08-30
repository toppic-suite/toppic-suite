//
// Created by abbash on 7/29/22.
//

#ifndef TOPPIC_SEED_ENVELOPE_HPP
#define TOPPIC_SEED_ENVELOPE_HPP

#include <vector>
#include <algorithm>
#include "simple_peak.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_base.hpp"

namespace toppic {
    class SeedEnvelope {
    public:
        SeedEnvelope();
        SeedEnvelope(DeconvMsPtr deconv_data);
        SeedEnvelope(int spec_id, int env_id, double pos, double mass, double inte, int charge, std::vector<double> pos_list, std::vector<double> inte_list);
        SeedEnvelope(const SeedEnvelope &env);

        std::vector<double> get_pos_list();
        std::vector<double> get_inte_list();
        double get_inte_sum();
        void rm_peaks(double min_pos, double  max_pos);
        double get_base_pos(); ///////////
        void keep_top_three();
        void change_charge(int new_charge);

        double get_max_pos(){ return peak_list_[peak_list_.size() - 1].getPos(); }
        int get_peak_num() { return peak_list_.size(); }
        SimplePeak get_peak(int idx) { return peak_list_[idx]; }
        double get_proton_mass(){ return 1.007276466879; }
        double get_isotope_mass(){ return 1.00235; }
        SeedEnvelope get_new_charge_env(int new_charge);
        void remove_low_inte_peaks(double ratio, double base_inte);
        void shift(int shift_num);
        SeedEnvelope get_shifted_seed_envelope(EnvBase env_base, int shift_num);

        static bool cmpInteDec(const SeedEnvelope &a, const SeedEnvelope &b) {
          return a.getInte() > b.getInte();}


        int getSpecId() const{ return spec_id_; }
        void setSpecId(int specId) { spec_id_ = specId; }

        int getEnvId() const { return env_id_; }
        void setEnvId(int envId) { env_id_ = envId; }

        double getPos() const { return pos_; }
        void setPos(double pos) { pos_ = pos; }

        double getMass() const { return mass_; }
        void setMass(double mass) { mass_ = mass; }

        double getInte() const { return inte_; }
        void setInte(double inte) { inte_ = inte; }

        int getCharge() const { return charge_; }
        void setCharge(int charge) { charge_ = charge; }

        const std::vector<SimplePeak> &getPeakList() const { return peak_list_; }
        void setPeakList(const std::vector<SimplePeak> &peakList) { peak_list_ = peakList; }

        bool isEmpty(){
          if (spec_id_ == -1 && env_id_ == -1 && pos_ == -1 && mass_ == -1 && inte_ == -1 && charge_ == -1 && peak_list_.size() == 0)
            return true;
          return false;
        }

    private:
        int spec_id_;
        int env_id_;
        double pos_;
        double mass_;
        double inte_;
        int charge_;
        std::vector<SimplePeak> peak_list_;
    };
}

#endif //TOPPIC_SEED_ENVELOPE_HPP
