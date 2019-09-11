//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <fstream>

#include "xml2json/rapidjson/document.h"
#include "xml2json/rapidjson/prettywriter.h"
#include "xml2json/rapidjson/stringbuffer.h"

#include "common/util/logger.hpp"
#include "deconv/msreader/raw_ms_writer.hpp" 

namespace toppic {

namespace raw_ms_writer {

void write(std::string &file_name, RawMsPtr ms_ptr, MatchEnvPtrVec &envs) {

  rapidjson::Document doc;

  // define the document as an object rather than an array
  doc.SetObject();
  
  // must pass an allocator when the object may need to allocate memory
  rapidjson::Document::AllocatorType& allocator = doc.GetAllocator();
  
  // create a rapidjson array type with similar syntax to std::vector
  rapidjson::Value peaks(rapidjson::kArrayType);

  PeakPtrVec raw_peaks = ms_ptr->getPeakPtrVec();
  for (size_t i = 0; i < raw_peaks.size(); i++) {
    rapidjson::Value peak(rapidjson::kObjectType);
    peak.AddMember("peak_id", i, allocator);
    peak.AddMember("mz", raw_peaks[i]->getPosition(), allocator);
    peak.AddMember("intensity", raw_peaks[i]->getIntensity(), allocator);
    peaks.PushBack(peak, allocator);
  }

  doc.AddMember("peaks", peaks, allocator);

  rapidjson::Value envelopes(rapidjson::kArrayType);
  for (size_t i = 0; i < envs.size(); i++) {
    rapidjson::Value env(rapidjson::kObjectType);
    EnvelopePtr theo_env = envs[i]->getTheoEnvPtr();
    env.AddMember("mono_mass", theo_env->getMonoNeutralMass(), allocator);
    env.AddMember("charge", theo_env->getCharge(), allocator);

    rapidjson::Value env_peaks(rapidjson::kArrayType);
    for (int k = 0; k < theo_env->getPeakNum(); k++) {
      rapidjson::Value peak(rapidjson::kObjectType);
      peak.AddMember("mz", theo_env->getMz(k), allocator);
      peak.AddMember("intensity", theo_env->getIntensity(k), allocator);
      env_peaks.PushBack(peak, allocator);
    }
    env.AddMember("env_peaks", env_peaks, allocator);
    envelopes.PushBack(env, allocator);
  }
  doc.AddMember("envelopes", envelopes, allocator);

  // Convert JSON document to string
  rapidjson::StringBuffer buffer;
  rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
  doc.Accept(writer);
  std::ofstream output;
  output.open(file_name.c_str());
  output << "spectrum_data =" << std::endl;
  output << buffer.GetString();
  output.close();
}

}

}
