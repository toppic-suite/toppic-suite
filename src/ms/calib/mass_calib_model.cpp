// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "ms/calib/mass_calib_model.hpp"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"

namespace toppic {

namespace {

// %.17g: every double is written so that it is read back unchanged
std::string fmt(double v) {
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

std::string levelJson(const MassCalibLevel& level) {
  return "{\"c0\": " + fmt(level.c0) + ", \"c_mz\": " + fmt(level.c_mz) +
         ", \"c_rt_per_hour\": " + fmt(level.c_rt_per_hour) + "}";
}

double getNum(const boost::property_tree::ptree& tree, const std::string& key,
              const std::string& file_name) {
  auto v = tree.get_optional<double>(key);
  if (!v) {
    LOG_ERROR("The calibration model " << file_name << " has no " << key
                                       << "!");
    std::exit(EXIT_FAILURE);
  }
  return *v;
}

}  // namespace

double MassCalibModel::getErrorPpm(int ms_level, double mz,
                                   double rt_min) const {
  const MassCalibLevel& c = ms_level <= 1 ? ms1 : ms2;
  return c.c0 + c.c_mz * (mz - mz0) / 1000.0 +
         c.c_rt_per_hour * (rt_min - rt0_min) / 60.0;
}

double MassCalibModel::correctMz(int ms_level, double mz,
                                 double rt_min) const {
  return mz * (1.0 - getErrorPpm(ms_level, mz, rt_min) * 1e-6);
}

double MassCalibModel::correctMass(int ms_level, double mass, int charge,
                                   double rt_min) const {
  double proton = mass_constant::getProtonMass();
  double mz = (mass + charge * proton) / charge;
  double e = getErrorPpm(ms_level, mz, rt_min) * 1e-6;
  return mass - e * (mass + charge * proton);
}

std::string MassCalibModel::toJson() const {
  std::ostringstream out;
  out << "{\n"
      << " \"mz0\": " << fmt(mz0) << ",\n"
      << " \"rt0_min\": " << fmt(rt0_min) << ",\n"
      << " \"ms1\": " << levelJson(ms1) << ",\n"
      << " \"ms2\": " << levelJson(ms2) << ",\n"
      << " \"n_prec\": " << prec_num << ",\n"
      << " \"n_frag\": " << frag_num << ",\n"
      << " \"n_frag_prsms\": " << frag_prsm_num << ",\n"
      << " \"robust_sd_prec\": " << fmt(prec_robust_sd) << ",\n"
      << " \"robust_sd_frag\": " << fmt(frag_robust_sd) << ",\n"
      << " \"note\": \"e (ppm) = c0 + c_mz (m/z - mz0)/1000 + c_rt_per_hour "
         "(rt_min - rt0_min)/60; corrected m/z = m/z (1 - e 1e-6)\"\n"
      << "}\n";
  return out.str();
}

void MassCalibModel::write(const std::string& file_name) const {
  std::ofstream out(file_name);
  if (!out) {
    LOG_ERROR("Cannot write the calibration model " << file_name << "!");
    std::exit(EXIT_FAILURE);
  }
  out << toJson();
}

MassCalibModel MassCalibModel::read(const std::string& file_name) {
  boost::property_tree::ptree tree;
  try {
    boost::property_tree::read_json(file_name, tree);
  } catch (const std::exception& e) {
    LOG_ERROR("Cannot read the calibration model " << file_name << ": "
                                                   << e.what());
    std::exit(EXIT_FAILURE);
  }
  MassCalibModel model;
  model.mz0 = getNum(tree, "mz0", file_name);
  model.rt0_min = getNum(tree, "rt0_min", file_name);
  MassCalibLevel* levels[2] = {&model.ms1, &model.ms2};
  const char* names[2] = {"ms1", "ms2"};
  for (int i = 0; i < 2; i++) {
    std::string n = names[i];
    levels[i]->c0 = getNum(tree, n + ".c0", file_name);
    levels[i]->c_mz = getNum(tree, n + ".c_mz", file_name);
    levels[i]->c_rt_per_hour = getNum(tree, n + ".c_rt_per_hour", file_name);
  }
  model.prec_num = tree.get<int>("n_prec", 0);
  model.frag_num = tree.get<int>("n_frag", 0);
  model.frag_prsm_num = tree.get<int>("n_frag_prsms", 0);
  model.prec_robust_sd = tree.get<double>("robust_sd_prec", 0.0);
  model.frag_robust_sd = tree.get<double>("robust_sd_frag", 0.0);
  return model;
}

}  // namespace toppic
