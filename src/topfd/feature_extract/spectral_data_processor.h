//
// Created by abbash on 9/17/2021.
//

#ifndef TOPPIC_SPECTRAL_DATA_PROCESSOR_H
#define TOPPIC_SPECTRAL_DATA_PROCESSOR_H

#include <algorithm>
#include "src/ms/spec/deconv_ms.hpp"
#include "src/topfd/common/topfd_para.hpp"
#include "ms/spec/deconv_peak.hpp"
#include "ms/spec/ms.hpp"
#include "ms/spec/simple_msalign_reader.hpp"
#include "topfd/msreader/raw_ms_reader.hpp"
#include "ms/spec/baseline_util.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_base.hpp"
#include "common/base/mass_constant.hpp"
#include "topfd/msreader/raw_ms_group_reader.hpp"
#include "ms/spec/raw_ms_util.hpp"
#include "common/util/logger.hpp"
#include "ms/env/env_filter.hpp"
#include "ms/env/env_assign.hpp"
#include "topfd/spec/deconv_data_util.hpp"
#include "topfd/dp/dp_a.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/env/match_env_refine.hpp"
#include "topfd/feature_extract/mzml_data_processor.h"

namespace toppic {

class SpectralProcessor {
 public:
  SpectralProcessor(std::string ms1_file_name, TopfdParaPtr topfd_para_ptr, mzmlProcessorPtr mzml_processor);

  DeconvPeakPtrVec shortlistPeaks(DeconvPeakPtrVec matched_peaks, int rt_begin_sp_id, int rt_end_sp_id);

  DeconvPeakPtrVec getMatchedPeaks(DeconvPeakPtr spectrum_peak, double tol, TopfdParaPtr topfd_para_ptr);

  DeconvPeakPtr get_spec_peak(int idx) { return spectrum_peaks_[idx]; }

  void removePeaks(DeconvPeakPtrVec shortlistMatchedPeaks);

  bool exists(DeconvPeakPtr peak);

  int getNumPeaks() { return spectrum_peaks_.size(); }

  void print_n(int x, DeconvPeakPtr peak) {LOG_ERROR("IDX: " << x << " and ID: " << peak->getId()); }

  void print(int x) {LOG_ERROR(x << " not cleared."); }

 private:
  int _binary_search(double d);

  std::vector<double> _getExtendMasses(double d);

  void get_envelope(DeconvPeakPtrVec peaks, size_t base, PeakPtrVec peak_list, EnvParaPtr env_para_ptr);

  bool _check_overlap(DeconvPeakPtr spectrum_peak, DeconvPeakPtr peak, TopfdParaPtr topfd_para_ptr);

  DeconvMsPtrVec spec_list_;
  DeconvPeakPtrVec spectrum_peaks_;
  DeconvPeakPtrVec spectrum_peaks_mass_sorted_;

};

typedef std::shared_ptr<SpectralProcessor> SpectralProcessorPtr;

} //namespace
#endif //TOPPIC_SPECTRAL_DATA_PROCESSOR_H
