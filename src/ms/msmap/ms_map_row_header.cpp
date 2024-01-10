//
// Created by abbash on 8/22/22.
//

#include "ms/msmap/ms_map_row_header.hpp"

namespace toppic {

MsMapRowHeader::MsMapRowHeader(int row_id, int raw_spec_id, int scan_num, double rt):
  row_id_(row_id),
  raw_spec_id_(raw_spec_id),
  scan_num_(scan_num),
  rt_(rt) {}

}

