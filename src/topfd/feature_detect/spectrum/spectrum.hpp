//
// Created by abbash on 8/22/22.
//

#ifndef TOPPIC_SPECTRUM_HPP
#define TOPPIC_SPECTRUM_HPP

#include <vector>

namespace toppic {
    class Spectrum {
    public:
        Spectrum(int spec_id, int scan_num, double rt){
          spec_id_ = spec_id;
          scan_num_ = scan_num;
          rt_ = rt;
        }

        Spectrum() {}

        int getSpecId() const { return spec_id_; }
        void setSpecId(int specId) { spec_id_ = specId; }

        int getScanNum() const { return scan_num_; }
        void setScanNum(int scanNum) { scan_num_ = scanNum; }

        double getRt() const { return rt_; }
        void setRt(double rt) { rt_ = rt; }

    private:
        int spec_id_;
        int scan_num_;
        double rt_;

    };
    typedef std::vector<Spectrum> spec_list;
}

#endif //TOPPIC_SPECTRUM_HPP
