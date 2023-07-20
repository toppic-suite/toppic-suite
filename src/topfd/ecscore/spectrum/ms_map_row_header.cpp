//
// Created by abbash on 8/22/22.
//

#include "topfd/ecscore/spectrum/ms_map_row_header.hpp"

namespace toppic {

MsMapRowHeader::MsMapRowHeader(int spec_id, int scan_num, double rt):
  spec_id_(spec_id),
  scan_num_(scan_num),
  rt_(rt) {}

}

