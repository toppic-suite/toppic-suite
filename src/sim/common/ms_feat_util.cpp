//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include "ms/msmap/ms_map.hpp"
#include "sim/common/ms_feat_util.hpp"
#include "sim/png/pngwriter.hpp"

namespace toppic {

namespace ms_feat_util {
  
void writepng(MsMapPtr map_ptr, const std::string& png_file_name) {
  int scale =  100;
  int width = 512;
  int height = 512;
  int center = width/2;
  double mz_win_size = width/scale;

  pngwriter png(width, height, 1.0, (png_file_name +"_data.png").c_str());

  double left_mz = 0;
  double right_mz = 2000;
  double ref_mz = (left_mz + right_mz) / 2.0;
  // get max intensity in m/z range
  double max_inte = 0;
  MsMapPeakPtr2D peak_2d_list = map_ptr->get2DPeaks();
  for (size_t i = 0; i < peak_2d_list.size(); i++) {
    for (size_t j = 0; j < peak_2d_list[i].size(); j++) {
      MsMapPeakPtr peak = peak_2d_list[i][j];
      double peak_pos = peak->getPosition();
      if (peak_pos >= left_mz && peak_pos < right_mz) {
        if (peak->getIntensity() > max_inte) {
          max_inte = peak->getIntensity();
        }
      }
    }
  }

  // get the background row ranges
  int spec_range = map_ptr->getRowNum();
  int png_start_spec = 0;
  int png_end_spec = spec_range - 1;
  // add background
  for (size_t i = png_start_spec; i <= png_end_spec; i++) {
    for (size_t j = 0; j < peak_2d_list[i].size(); j++) {
      MsMapPeakPtr peak = peak_2d_list[i][j];
      double peak_pos = peak->getPosition();
      if (peak_pos >= left_mz && peak_pos < right_mz) {
        int x = (peak_pos - ref_mz)*scale + center;
        int y = i - png_start_spec;
        double inte = peak->getIntensity() / max_inte;
        png.plot(x, y, inte, 0.0, 0.0);
      }
    }
  }

  png.close();
}

}

}