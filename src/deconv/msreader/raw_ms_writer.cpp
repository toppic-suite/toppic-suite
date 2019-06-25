//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

void write(std::string &file_name, RawMsPtr ms_ptr) {

  rapidjson::Document doc;

  // define the document as an object rather than an array
  doc.SetObject();
  
  // must pass an allocator when the object may need to allocate memory
  rapidjson::Document::AllocatorType& allocator = doc.GetAllocator();
  
  // create a rapidjson array type with similar syntax to std::vector
  rapidjson::Value peaks(rapidjson::kArrayType);

  PeakPtrVec raw_peaks = ms_ptr->getPeakPtrVec();
  for (size_t i = 0; i < raw_peaks.size(); i++) {
    rapidjson::Value object(rapidjson::kObjectType);
    object.AddMember("peak_id", i, allocator);
    object.AddMember("mz", raw_peaks[i]->getPosition(), allocator);
    object.AddMember("intensity", raw_peaks[i]->getIntensity(), allocator);
    peaks.PushBack(object, allocator);
  }

  doc.AddMember("peaks", peaks, allocator);

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
