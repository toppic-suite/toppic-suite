//
// Created by abbash on 8/22/22.
//

#include "topfd/ecscore/spectrum/matrix_spectrum.hpp"

namespace toppic {

MatrixSpectrum::MatrixSpectrum(int spec_id, int scan_num, double rt): 
  spec_id_(spec_id),
  scan_num_(scan_num),
  rt_(rt) {}

}

