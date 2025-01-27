// Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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


#include <algorithm>
#include <numeric>

#include "common/util/logger.hpp"
#include "topdia/pseudo_spec/mzrt_feature.hpp"
#include "topdia/pseudo_spec/generate_pseudo_spectrum.hpp"

namespace toppic {
MzrtFeature::MzrtFeature(int id, int fraction_id, int env_num, double mass,
                         double mono_mz, int charge, double intensity,
                         int mz_begin, int mz_end, double time_begin,
                         double time_end, int spec_id_begin, int spec_id_end,
                         double time_apex, double ec_score,
                         std::vector<double> xic,
                         std::vector<double> normalized_xic,
                         std::vector<double> envelope_mz,
                         std::vector<double> envelope_inte, int apex_cycle) {
  id_ = id;
  fraction_id_ = fraction_id;
  env_num_ = env_num;
  mass_ = mass;
  mono_mz_ = mono_mz;
  charge_ = charge;
  intensity_ = intensity;
  mz_begin_ = mz_begin;
  mz_end_ = mz_end;
  time_begin_ = time_begin;
  time_end_ = time_end;
  time_apex_ = time_apex;
  spec_id_begin_ = spec_id_begin;
  spec_id_end_ = spec_id_end;
  ec_score_ = ec_score;
  xic_ = xic;
  apex_cycle_ = apex_cycle;
  envelope_mz_ = envelope_mz;
  envelope_inte_ = envelope_inte;
  normalized_xic_ = normalized_xic;
  used_ = false;
  cycle_span_ = countNonZero(xic);
}

MzrtFeaturePtrVec MzrtFeature::read_record(std::string filename) {
  MzrtFeaturePtrVec data;  // Store read data
  std::ifstream file(filename);
  if (!file.is_open()) {
    LOG_ERROR("Error opening file: " << filename); 
    return data;  // Return empty vector if file cannot be opened
  }

  std::string line;
  std::getline(file, line);
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> tokens;
    while (std::getline(ss, cell, ',')) {
      tokens.push_back(cell);
    }
    if (tokens.size() >= 0) {
      int id = std::stoi(tokens[0]);
      int fraction_id = std::stoi(tokens[1]);
      int env_num = std::stoi(tokens[2]);
      double mass = std::stod(tokens[3]);
      double mono_mz = std::stod(tokens[4]);
      int charge = std::stoi(tokens[5]);
      double intensity = std::stod(tokens[6]);
      int mz_begin = std::stod(tokens[7]);
      int mz_end = std::stod(tokens[8]);
      double time_begin = std::stod(tokens[9]);
      double time_end = std::stod(tokens[10]);
      int spec_id_begin = std::stoi(tokens[11]);
      int spec_id_end = std::stoi(tokens[12]);
      double time_apex = std::stod(tokens[13]);
      double ec_score = std::stod(tokens[14]);
      std::vector<double> xic = parseXIC(tokens[15]);
      std::vector<double> normalized_xic = normalizeXIC(xic);
      std::vector<double> envelope_mz, envelope_inte;
      parseEnvelope(tokens[16], envelope_mz, envelope_inte);
      int apex_cycle =
          std::distance(xic.begin(), std::max_element(xic.begin(), xic.end()));
      MzrtFeaturePtr feature = std::make_shared<MzrtFeature>(
          id, fraction_id, env_num, mass, mono_mz, charge, intensity, mz_begin,
          mz_end, time_begin, time_end, spec_id_begin, spec_id_end, time_apex,
          ec_score, xic, normalized_xic, envelope_mz, envelope_inte,
          apex_cycle);
      if (std::accumulate(xic.begin(), xic.end(), 0) == 0) continue;
      data.push_back(feature);
    }
  }
  file.close();
  std::sort(data.begin(), data.end(),
            GeneratePseudoSpectrum::compareFeaturesInte);
  return data;
}

std::vector<double> MzrtFeature::parseXIC(const std::string &line) {
  std::vector<double> result;
  std::stringstream ss(line);
  std::string cell;
  while (std::getline(ss, cell, ';')) {
    result.push_back(std::stod(cell));
  }
  return result;
}

std::vector<double> MzrtFeature::normalizeXIC(const std::vector<double> &xic) {
  // Calculate the sum of all elements in xic
  double sum = 0.0;
  for (const auto &value : xic) {
    sum += value;
  }
  // Normalize the vector by dividing each element by the sum
  std::vector<double> normalized_xic;
  for (auto &value : xic) {
    normalized_xic.push_back(value / sum);
  }
  return normalized_xic;
}

void MzrtFeature::parseEnvelope(const std::string &input,
                                std::vector<double> &array1,
                                std::vector<double> &array2) {
  std::stringstream ss(input);
  std::string pair;
  while (std::getline(ss, pair, ';')) {
    std::stringstream pairStream(pair);
    std::string value1, value2;
    std::getline(pairStream, value1, '&');
    std::getline(pairStream, value2, '&');
    array1.push_back(std::stod(value1));
    array2.push_back(std::stod(value2));
  }
}

int MzrtFeature::countNonZero(const std::vector<double> &xic) {
  int nonZeroCount = 0;
  for (int num : xic)
    if (num != 0) nonZeroCount++;
  return nonZeroCount;
}

}  // namespace toppic
