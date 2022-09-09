//
// Created by abbash on 7/29/22.
//

#include "seed_envelope.hpp"

std::shared_ptr<toppic::Envelope> get_env_(double mono_mass, int charge, double mono_mz){
  toppic::EnvelopePtr ref_env_ptr = toppic::EnvBase::getStaticEnvByMonoMass(mono_mass);
  if (ref_env_ptr == nullptr) {
    return nullptr;
  }
  toppic::EnvelopePtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mono_mz, charge);
  return theo_env_ptr;
}

toppic::SeedEnvelope::SeedEnvelope(DeconvPeakPtr& p){
  spec_id_ = p->getSpId();
  env_id_ = p->getId();
  mass_ = p->getMonoMass();
  pos_ = p->getMonoMz();
  inte_ = p->getIntensity();
  charge_ = p->getCharge();
  std::shared_ptr<toppic::Envelope> theo_env =  get_env_(mass_, charge_, pos_);
  for (int j = 0; j < theo_env->getPeakNum(); j++)
    peak_list_.push_back(SimplePeak(theo_env->getMz(j), theo_env->getIntensity(j)));
}

toppic::SeedEnvelope::SeedEnvelope() {
  spec_id_ = -1;
  env_id_ = -1;
  pos_ = -1;
  mass_ = -1;
  inte_ = -1;
  charge_ = -1;
}

toppic::SeedEnvelope::SeedEnvelope(int spec_id, int env_id, double pos, double mass, double inte, int charge, std::vector<double> pos_list, std::vector<double> inte_list) {
  spec_id_ = spec_id;
  env_id_ = env_id;
  pos_ = pos;
  mass_ = mass;
  inte_ = inte;
  charge_ = charge;
  for (int i = 0; i < pos_list.size(); i++)
    peak_list_.push_back(SimplePeak(pos_list[i], inte_list[i]));
}

toppic::SeedEnvelope::SeedEnvelope(const toppic::SeedEnvelope &env) {
  spec_id_ = env.spec_id_;
  env_id_ = env.env_id_;
  pos_ = env.pos_;
  mass_ = env.mass_;
  inte_ = env.inte_;
  charge_ = env.charge_;
  for (const auto & i : env.peak_list_)
    peak_list_.push_back(SimplePeak(i.getPos(), i.getInte()));
}

std::vector<double> toppic::SeedEnvelope::get_pos_list() {
  std::vector<double> pos_list;
  for (auto p : peak_list_)
    pos_list.push_back(p.getPos());
  return pos_list;
}

std::vector<double> toppic::SeedEnvelope::get_inte_list() {
  std::vector<double> inte_list;
  for (auto p : peak_list_)
    inte_list.push_back(p.getInte());
  return inte_list;
}

double toppic::SeedEnvelope::get_inte_sum(){
  double inte_sum = 0;
  for (auto p : peak_list_)
    inte_sum = inte_sum + p.getInte();
  return inte_sum;
}

void toppic::SeedEnvelope::rm_peaks(double min_pos, double max_pos){
  std::vector<SimplePeak> peak_list;
  for (auto p : peak_list_) {
    if (p.getPos() >= min_pos and p.getPos() <= max_pos)
      peak_list.push_back(p);
    peak_list_ = peak_list;
  }
}

void toppic::SeedEnvelope::change_charge(int new_charge){
  for (auto &p : peak_list_) {
    double new_pos = (((p.getPos() - get_proton_mass()) * charge_)/new_charge) + get_proton_mass();
    p.setPos(new_pos);
  }
  pos_ = (((pos_ - get_proton_mass()) * charge_)/new_charge) + get_proton_mass();
  charge_ = new_charge;
}

toppic::SeedEnvelope toppic::SeedEnvelope::get_new_charge_env(int new_charge){
  SeedEnvelope new_env = *this;
  new_env.change_charge(new_charge);
  return new_env;
}

void toppic::SeedEnvelope::remove_low_inte_peaks(double ratio, double base_inte, double snr) {
  std::vector<SimplePeak> peak_list;
  for (auto &p : peak_list_) {
    if (p.getInte() * ratio >= (base_inte * snr))
      peak_list.push_back(p);
    peak_list_ = peak_list;
  }
}

void toppic::SeedEnvelope::shift(double shift_num){
  double shift_mass = shift_num * get_isotope_mass();
  double shift_mz = shift_mass/charge_;
  pos_ = pos_ + shift_mz;
  mass_ = mass_ + shift_mass;
  for (auto p : peak_list_)
    p.setPos(p.getPos()+shift_mz);
}

void toppic::SeedEnvelope::keep_top_three(){
  std::sort(peak_list_.begin(), peak_list_.end(), SimplePeak::cmpInteDec);
  std::vector<SimplePeak> vec (peak_list_.begin(), peak_list_.begin() + 3);
  peak_list_ = vec;
  std::sort(peak_list_.begin(), peak_list_.end(), SimplePeak::cmpPos);
}

double toppic::SeedEnvelope::get_base_pos() {
  std::sort(peak_list_.begin(), peak_list_.end(), SimplePeak::cmpInteDec);
  double base_pos = peak_list_[0].getPos();
  std::sort(peak_list_.begin(), peak_list_.end(), SimplePeak::cmpPos);
  return base_pos;
}

toppic::SeedEnvelope toppic::SeedEnvelope::get_shifted_seed_envelope(EnvBase env_base, int shift_num){
  double mass = mass_ + (shift_num * get_isotope_mass());
  double mz = (mass + (charge_ * get_proton_mass())) / charge_;
  std::shared_ptr<toppic::Envelope> env = get_env_(mass, charge_, mz);
  std::vector<double> env_peaks_mz;
  std::vector<double> env_peaks_inte;

  for (int i = 0; i < env->getPeakNum(); i++){
    env_peaks_mz.push_back(env->getMz(i));
    env_peaks_inte.push_back(env->getIntensity(i));
  }
  SeedEnvelope sp_peak = SeedEnvelope(spec_id_, env_id_, mz, mass, inte_, charge_, env_peaks_mz, env_peaks_inte);
  return sp_peak;
}