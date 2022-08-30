//
// Created by abbash on 8/23/22.
//

#include "peak.hpp"

int toppic::Peak::getSpecId() const { return spec_id_; }
void toppic::Peak::setSpecId(int specId) { spec_id_ = specId; }

double toppic::Peak::getPos() const { return pos_; }
void toppic::Peak::setPos(double pos) { pos_ = pos; }

double toppic::Peak::getInte() const { return inte_; }
void toppic::Peak::setInte(double inte) { inte_ = inte; }

double toppic::Peak::getOriInte() const { return ori_inte_; }
void toppic::Peak::setOriInte(double oriInte) { ori_inte_ = oriInte; }
