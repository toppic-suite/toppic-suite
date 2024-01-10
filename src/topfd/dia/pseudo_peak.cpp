//
// Created by abbash on 12/30/23.
//

#include "topfd/dia/pseudo_peak.hpp"

namespace toppic {

PseudoPeak::PseudoPeak(double mass, double mono_mz, int charge, 
                       double intensity, double score, double corr, 
                       double coverage, double apex_diff, double apex_diff_scan, 
                       double rt_low, double rt_high): 
  mass_(mass), mono_mz_(mono_mz),
  charge_(charge),
  intensity_(intensity),
  score_(score), corr_(corr),
  coverage_(coverage),
  apex_diff_(apex_diff),
  apex_diff_scan_(apex_diff_scan),
  rt_low_(rt_low), 
  rt_high_(rt_high) {}
}
