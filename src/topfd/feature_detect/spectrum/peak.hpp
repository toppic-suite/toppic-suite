//
// Created by abbash on 8/23/22.
//

#ifndef TOPPIC_PEAK_HPP
#define TOPPIC_PEAK_HPP

namespace toppic {
    class Peak {
    public:
        int getSpecId() const;
        void setSpecId(int specId);

        double getPos() const;
        void setPos(double pos);

        double getInte() const;
        void setInte(double inte);

        double getOriInte() const;
        void setOriInte(double oriInte);

    private:
        int spec_id_;
        double pos_;
        double inte_;
        double ori_inte_;

    };
}


#endif //TOPPIC_PEAK_HPP
