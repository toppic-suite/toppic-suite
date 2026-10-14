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

#include "console/topdia_argument.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/mem_check.hpp"
#include "common/util/version.hpp"

namespace toppic {

namespace {

// Parse `str` as a double, rejecting empty strings and trailing characters
// (e.g. "1.5x" or "3abc"). Returns true and sets `val` only on a full parse.
bool toDouble(const std::string& str, double& val) {
  try {
    size_t pos = 0;
    val = std::stod(str, &pos);
    return pos == str.size();
  } catch (const std::exception&) {
    return false;
  }
}

// Parse `str` as an int, rejecting empty strings and trailing characters.
// Returns true and sets `val` only on a full parse.
bool toInt(const std::string& str, int& val) {
  try {
    size_t pos = 0;
    val = std::stoi(str, &pos);
    return pos == str.size();
  } catch (const std::exception&) {
    return false;
  }
}

// The deconvolution defaults of version 1.9.0.0, which topdia keeps; see
// getTopfdParaPtrForTopdia().
const double TOPDIA_MS2_RANK_PEAK_DISTANCE = 12;
const double TOPDIA_MASS_GROUP_BOUNDARY = 1500;

}  // namespace

Argument::Argument() {
  topfd_para_ptr_ = getTopfdParaPtrForTopdia();
  topdia_para_ptr_ = std::make_shared<TopdiaPara>();
}

TopfdParaPtr Argument::getTopfdParaPtrForTopdia() {
  TopfdParaPtr topfd_para_ptr = std::make_shared<TopfdPara>();
  topfd_para_ptr->setMs1EcscoreCutoff(0);
  topfd_para_ptr->setMs1MinScanNum(2);
  topfd_para_ptr->setPrecWindowWidth(4.0);
  topfd_para_ptr->setOutputCsvFeatureFile(true);
  // TopDIA keeps the deconvolution defaults of version 1.9.0.0. The three
  // defaults TopFD changed in 2.0.0.0 - the 2 Dalton MS/MS envelope ranking
  // window, the MS/MS harmonic check and the 4250 Dalton mass group boundary -
  // were measured on DDA data. A DIA MS/MS scan holds the fragments of every
  // co-isolated precursor of a wide isolation window, so they roughly double
  // the reported fragment masses and with them the fragment features that
  // pseudo spectrum generation holds in memory; on the tutorial data set that
  // was enough to make topdia run out of memory (experiment 61 of
  // toppic_test). The hidden --ms2-rank-peak-distance and
  // --disable-ms2-harmonic-check options still override these.
  topfd_para_ptr->setMsTwoRankPeakDistance(TOPDIA_MS2_RANK_PEAK_DISTANCE);
  topfd_para_ptr->setMsTwoHarmonicCheck(false);
  topfd_para_ptr->setMassGroupBoundary(TOPDIA_MASS_GROUP_BOUNDARY);
  // TopDIA also keeps MS-Deconv's MS1 and MS/MS deconvolution by default:
  // TopFD's default greedy searches were measured on DDA data only
  // (toppic_test experiments 104, 106 and 109). --ms1-deconv and
  // --ms2-deconv select them.
  topfd_para_ptr->setMs1DeconvMethod("msdeconv");
  topfd_para_ptr->setMs2DeconvMethod("msdeconv");
  // nor does it correct precursor isotopes with fragment pairs (TopFD's
  // default, measured on DDA data in experiment 116): the fragments of a DIA
  // MS/MS scan come from every co-isolated precursor.
  topfd_para_ptr->setMs2PrecPairRefine(false);
  return topfd_para_ptr;
}

void Argument::showUsage(
    const boost::program_options::options_description& desc) {
  std::cout << "Usage: topdia [options] spectrum-file-names" << std::endl;
  std::cout << desc << std::endl;
  std::cout << "Version: " << Version::getVersion() << std::endl;
}

bool Argument::parse(int argc, char* argv[]) {
  std::string max_charge = "";
  std::string max_mass = "";
  std::string mz_error = "";
  std::string ms1_rank_peak_distance = "";
  std::string ms2_rank_peak_distance = "";
  std::string ms_two_sn_ratio = "";
  std::string ms_one_sn_ratio = "";
  std::string prec_window = "";
  std::string thread_number = "";
  std::string activation = "";

  std::string ms1_ecscore_cutoff = "";
  std::string ms2_ecscore_cutoff = "";
  std::string ms1_min_scan_num = "";
  std::string ms2_min_scan_num = "";
  std::string pseudo_score_cutoff = "";
  std::string pseudo_min_peaks = "";
  std::string ms1_seed_env_inte_corr_tole_cutoff = "";
  std::string ms2_seed_env_inte_corr_tole_cutoff = "";
  // deconvolution methods
  std::string deconv_ms1 = "";
  std::string deconv_ms2 = "";
  // the score that ranks the envelopes
  std::string env_rank = "";
  // the EnvCNN score cutoffs of the MS1 and MS/MS envelopes
  std::string ms1_env_cnn_score_cutoff = "";
  std::string ms2_env_cnn_score_cutoff = "";

  // Define and parse the program options
  try {
    namespace po = boost::program_options;
    // Documented options shown in the help message.
    po::options_description visible_desc("Options");
    visible_desc.add_options()("help,h", "Print this help message.")(
        "activation,a", po::value<std::string>(&activation),
        "<CID|ETD|ETCID|ETHCD|HCD|MPD|UVPD|FILE>. Fragmentation method of "
        "MS/MS spectra. When FILE is used, the fragmentation methods of "
        "spectra are given in the input spectral data file. ETCID (ETciD) and "
        "ETHCD (EThcD) have the ion types of ETD. The default value is FILE.")(
        "max-charge,c", po::value<std::string>(&max_charge),
        "<a positive integer>. Set the maximum charge state of precursor and "
        "fragment ions. The default value is 30.")(
        "max-mass,m", po::value<std::string>(&max_mass),
        "<a positive number>. Set the maximum monoisotopic mass of precursor "
        "and fragment ions. The default value is 50,000 Daltons.")(
        "mz-error,e", po::value<std::string>(&mz_error),
        "<a positive number>. Set the error tolerance of m/z values of "
        "spectral peaks. The default value is 0.02 m/z.")(
        "ms1-sn-ratio,r", po::value<std::string>(&ms_one_sn_ratio),
        "<a non-negative number>. Set the signal-to-noise ratio for MS1 "
        "spectra. The default value is 3.")(
        "ms2-sn-ratio,s", po::value<std::string>(&ms_two_sn_ratio),
        "<a non-negative number>. Set the signal-to-noise ratio for MS/MS "
        "spectra. The default value is 1.")(
        "missing-ms1,o", "MS1 spectra are missing in the input file.")(
        "precursor-window,w", po::value<std::string>(&prec_window),
        "<a positive number>. Set the default precursor window size. The "
        "default value is 4.0 m/z. When the input file contains precursor "
        "window information, the parameter is ignored.")(
        "ms1-deconv,p", po::value<std::string>(&deconv_ms1),
        "<greedy|msdeconv>. Set the deconvolution method of MS1 spectra: "
        "msdeconv, the MS-Deconv dynamic programming, or greedy, TopFD's "
        "default greedy search, which decides the charge state of the most "
        "intense peak with all the charge states of its mass and subtracts "
        "the mass's envelopes. The default value is msdeconv.")(
        "ms2-deconv,q", po::value<std::string>(&deconv_ms2),
        "<greedy|msdeconv|flashdeconv>. Set the deconvolution method of "
        "MS/MS spectra: msdeconv, the MS-Deconv dynamic programming; greedy, "
        "TopFD's default greedy search adapted to fragment ions; or "
        "flashdeconv, the FLASHDeconv algorithm (experimental). The default "
        "value is msdeconv.")(
        "ms1-envcnn-cutoff,f",
        po::value<std::string>(&ms1_env_cnn_score_cutoff),
        "<a number in [0,1]>. Set the EnvCNN score cutoff for isotopic "
        "envelopes in MS1 spectra: envelopes with a lower EnvCNN score are "
        "not reported, except those added for precursor isolation windows. "
        "Used with the greedy and msdeconv MS1 deconvolution methods. The "
        "default value is 0.2.")(
        "ms2-envcnn-cutoff,g",
        po::value<std::string>(&ms2_env_cnn_score_cutoff),
        "<a number in [0,1]>. Set the EnvCNN score cutoff for isotopic "
        "envelopes in MS/MS spectra: envelopes with a lower EnvCNN score are "
        "not reported. Used with the greedy and msdeconv MS/MS deconvolution "
        "methods, not with flashdeconv. The default value is 0.")(
        "env-rank,n", po::value<std::string>(&env_rank),
        "<envcnn|msdeconv>. Set the score used to rank isotopic envelopes: "
        "envcnn, the EnvCNN score, or msdeconv, the MS-Deconv score. The "
        "ranking decides which envelopes -j keeps; the EnvCNN cutoffs (-f "
        "and -g) use the EnvCNN score with either value. The default value "
        "is envcnn.")(
        "frag-num-filtering,j",
        "Filter fragment ion envelopes in MS/MS scans based on the estimated "
        "number of fragment ions: the envelopes ranked highest by the score "
        "of -n are kept. Used with the greedy and msdeconv MS/MS "
        "deconvolution methods, not with flashdeconv.")(
        "ms1-ecscore-cutoff,t", po::value<std::string>(&ms1_ecscore_cutoff),
        "<a number in [0,1]>. Set the MS1 ECScore cutoff value for proteoform "
        "features. The default value is 0.")(
        "ms2-ecscore-cutoff,T", po::value<std::string>(&ms2_ecscore_cutoff),
        "<a number in [0,1]>. Set the MS2 ECScore cutoff value for fragment "
        "features. The default value is 0.")(
        "ms1-min-scan-number,b", po::value<std::string>(&ms1_min_scan_num),
        "<1|2|3>. The minimum number of MS1 scans in which a proteoform "
        "feature is detected. The default value is 2.")(
        "ms2-min-scan-number,B", po::value<std::string>(&ms2_min_scan_num),
        "<1|2|3>. The minimum number of MS2 scans in which a fragment "
        "feature is detected. The default value is 1.")(
        "single-scan-noise,i",
        "Use the peak intensity noise levels in single scans (MS1 scans in "
        "proteoform feature detection, MS/MS scans in fragment feature "
        "detection) to filter out low intensity peaks. The default method is "
        "to use the peak intensity noise level of the whole LC-MS map (for "
        "fragment features, of the MS/MS scans of an isolation window) to "
        "filter out low intensity peaks.")(
        "pseudo-cutoff,v", po::value<std::string>(&pseudo_score_cutoff),
        "<a number in [0,1]>. Set the Pseudo Score cutoff value for "
        "generating pseudo-MS/MS spectra: once a pseudo-MS/MS spectrum has "
        "the number of peaks of -V, fragment features with a lower Pseudo "
        "Score are not added. The default value is 0.55.")(
        "pseudo-peak-number,V", po::value<std::string>(&pseudo_min_peaks),
        "<an integer of at least 10>. The minimum number of peaks in a "
        "pseudo-MS/MS spectrum: until a spectrum has this number of peaks, "
        "fragment features are added regardless of the Pseudo Score cutoff "
        "(-v), when available. The default value is 25.")(
        "ms1-correlation-cutoff,y",
        po::value<std::string>(&ms1_seed_env_inte_corr_tole_cutoff),
        "<a number in [0,1]>. Set the MS1 seed envelope intensity correlation "
        "cutoff value for extracting features. The default value is 0.5.")(
        "ms2-correlation-cutoff,Y",
        po::value<std::string>(&ms2_seed_env_inte_corr_tole_cutoff),
        "<a number in [0,1]>. Set the MS2 seed envelope intensity correlation "
        "cutoff value for extracting features. The default value is 0.5.")(
        "thread-number,u", po::value<std::string>(&thread_number),
        "<a positive integer>. Number of threads used in spectral "
        "deconvolution. The default value is 1.")(
        "no-sql,z", "Do not write the deconvoluted spectra to an SQLite database.");

    // Advanced options accepted on the command line but hidden from the help
    // message; the positional spectrum file argument is also hidden here.
    po::options_description hidden_desc("Hidden options");
    hidden_desc.add_options()(
        "disable-zstd",
        "Write the SQLite database uncompressed. By default its peak and "
        "envelope tables are compressed with zstd when it is complete "
        "(the format of the sqlite-zstd extension; topzstd -d decompresses "
        "it).")(
        "ms1-rank-peak-distance",
        po::value<std::string>(&ms1_rank_peak_distance),
        "<a positive number>. The window of the filtering-by-mz step for MS1 "
        "spectra, in Daltons: only the best scoring envelope of one charge "
        "within the window is kept. The default value is 12.")(
        "ms2-rank-peak-distance",
        po::value<std::string>(&ms2_rank_peak_distance),
        "<a positive number>. The same window for MS/MS spectra. The default "
        "value is 12.")(
        "disable-ms2-harmonic-check",
        "Skip the harmonic check of MS/MS spectra, which replaces an envelope "
        "by a better scoring harmonic parent mass. TopDIA skips it by "
        "default.")(
        "keep,k",
        "Report monoisotopic masses extracted from low quality isotopic "
        "envelopes.")(
        "spectrum-file-name",
        po::value<std::vector<std::string> >()->multitoken()->required(),
        "Spectrum file name with its path.");

    // All options (visible + hidden) are used for parsing; only visible_desc
    // is shown to the user, so the help text and the parser cannot drift apart.
    po::options_description desc("All options");
    desc.add(visible_desc).add(hidden_desc);

    po::positional_options_description positional_options;
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv)
                    .options(desc)
                    .positional(positional_options)
                    .run(),
                vm);
      if (vm.count("help")) {
        showUsage(visible_desc);
        return false;
      }
      po::notify(vm);
      // throws on error, so do after help in case there are any problems
    } catch (boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(visible_desc);
      return false;
    } catch (boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(visible_desc);
      return false;
    }

    // get the execution directory
    std::string argv_0(argv[0]);
    std::string exec_dir = file_util::getExecutiveDir(argv_0);
    if (file_util::checkSpace(exec_dir)) {
      LOG_ERROR("Current directory "
                << exec_dir
                << " contains space and will cause errors in the program!");
      exit(EXIT_FAILURE);
    }

    topfd_para_ptr_->setResourceDir(file_util::getResourceDir(exec_dir));

    if (vm.count("activation")) {
      topfd_para_ptr_->setActivation(activation);
    }

    if (vm.count("max-charge")) {
      int charge = 0;
      if (!toInt(max_charge, charge) || charge <= 0) {
        LOG_ERROR("Max charge " << max_charge
                                << " should be a positive integer.");
        return false;
      }
      topfd_para_ptr_->setMaxCharge(charge);
    }

    if (vm.count("keep")) {
      topfd_para_ptr_->setKeepUnusedPeaks(true);
    }

    if (vm.count("ms1-rank-peak-distance")) {
      double distance = 0;
      if (!toDouble(ms1_rank_peak_distance, distance) ||
          !std::isfinite(distance) || distance <= 0) {
        LOG_ERROR("MS1 rank peak distance "
                  << ms1_rank_peak_distance << " should be a positive number.");
        return false;
      }
      topfd_para_ptr_->setMsOneRankPeakDistance(distance);
    }

    if (vm.count("ms2-rank-peak-distance")) {
      double distance = 0;
      if (!toDouble(ms2_rank_peak_distance, distance) ||
          !std::isfinite(distance) || distance <= 0) {
        LOG_ERROR("MS/MS rank peak distance "
                  << ms2_rank_peak_distance << " should be a positive number.");
        return false;
      }
      topfd_para_ptr_->setMsTwoRankPeakDistance(distance);
    }

    if (vm.count("disable-ms2-harmonic-check")) {
      topfd_para_ptr_->setMsTwoHarmonicCheck(false);
    }

    if (vm.count("env-rank") && !topfd_para_ptr_->setDeconvScore(env_rank)) {
      LOG_ERROR("Envelope ranking score " << env_rank
                << " should be envcnn or msdeconv.");
      return false;
    }

    if (vm.count("ms1-envcnn-cutoff")) {
      double cutoff = 0;
      if (!toDouble(ms1_env_cnn_score_cutoff, cutoff) || cutoff < 0 ||
          cutoff > 1) {
        LOG_ERROR("MS1 EnvCNN cutoff " << ms1_env_cnn_score_cutoff
                                       << " should be a number in [0,1].");
        return false;
      }
      topfd_para_ptr_->setMs1EnvCnnScoreCutoff(cutoff);
    }

    if (vm.count("ms2-envcnn-cutoff")) {
      double cutoff = 0;
      if (!toDouble(ms2_env_cnn_score_cutoff, cutoff) || cutoff < 0 ||
          cutoff > 1) {
        LOG_ERROR("MS/MS EnvCNN cutoff " << ms2_env_cnn_score_cutoff
                                         << " should be a number in [0,1].");
        return false;
      }
      topfd_para_ptr_->setMs2EnvCnnScoreCutoff(cutoff);
    }

    if (vm.count("ms1-deconv") &&
        !topfd_para_ptr_->setMs1DeconvMethod(deconv_ms1)) {
      LOG_ERROR("MS1 deconvolution method " << deconv_ms1
                << " should be greedy or msdeconv.");
      return false;
    }

    if (vm.count("ms2-deconv") &&
        !topfd_para_ptr_->setMs2DeconvMethod(deconv_ms2)) {
      LOG_ERROR("MS/MS deconvolution method " << deconv_ms2
                << " should be greedy, msdeconv or flashdeconv.");
      return false;
    }

    if (vm.count("max-mass")) {
      double mass = 0;
      if (!toDouble(max_mass, mass) || mass <= 0) {
        LOG_ERROR("Max mass " << max_mass << " should be a positive number.");
        return false;
      }
      topfd_para_ptr_->setMaxMass(mass);
    }

    if (vm.count("mz-error")) {
      double error = 0;
      if (!toDouble(mz_error, error) || error <= 0) {
        LOG_ERROR("M/z error " << mz_error << " should be a positive number.");
        return false;
      }
      topfd_para_ptr_->setMzError(error);
    }

    if (vm.count("ms2-sn-ratio")) {
      double sn_ratio = 0;
      if (!toDouble(ms_two_sn_ratio, sn_ratio) || sn_ratio < 0) {
        LOG_ERROR("MS/MS S/N ratio " << ms_two_sn_ratio
                                     << " should be a non-negative number.");
        return false;
      }
      topfd_para_ptr_->setMsTwoSnRatio(sn_ratio);
    }

    if (vm.count("ms1-sn-ratio")) {
      double sn_ratio = 0;
      if (!toDouble(ms_one_sn_ratio, sn_ratio) || sn_ratio < 0) {
        LOG_ERROR("MS1 S/N ratio " << ms_one_sn_ratio
                                   << " should be a non-negative number.");
        return false;
      }
      topfd_para_ptr_->setMsOneSnRatio(sn_ratio);
    }

    if (vm.count("missing-ms1")) {
      topfd_para_ptr_->setMissingLevelOne(true);
    }

    if (vm.count("precursor-window")) {
      double window = 0;
      if (!toDouble(prec_window, window) || window <= 0) {
        LOG_ERROR("Precursor window " << prec_window
                                      << " should be a positive number.");
        return false;
      }
      topfd_para_ptr_->setPrecWindowWidth(window);
    }

    if (vm.count("ms1-ecscore-cutoff")) {
      double cutoff = 0;
      if (!toDouble(ms1_ecscore_cutoff, cutoff) || cutoff < 0 || cutoff > 1) {
        LOG_ERROR("MS1 ECScore cutoff " << ms1_ecscore_cutoff
                                        << " should be a number in [0,1].");
        return false;
      }
      topfd_para_ptr_->setMs1EcscoreCutoff(cutoff);
    }

    if (vm.count("ms2-ecscore-cutoff")) {
      double cutoff = 0;
      if (!toDouble(ms2_ecscore_cutoff, cutoff) || cutoff < 0 || cutoff > 1) {
        LOG_ERROR("MS2 ECScore cutoff " << ms2_ecscore_cutoff
                                        << " should be a number in [0,1].");
        return false;
      }
      topfd_para_ptr_->setMs2EcscoreCutoff(cutoff);
    }

    if (vm.count("pseudo-cutoff")) {
      double cutoff = 0;
      if (!toDouble(pseudo_score_cutoff, cutoff) || cutoff < 0 || cutoff > 1) {
        LOG_ERROR("Pseudo score cutoff " << pseudo_score_cutoff
                                         << " should be a number in [0,1].");
        return false;
      }
      topdia_para_ptr_->setPseudoScoreCutoff(cutoff);
    }

    if (vm.count("single-scan-noise")) {
      topfd_para_ptr_->setUseSingleScanNoiseLevel(true);
    }

    if (vm.count("ms1-min-scan-number")) {
      int n = 0;
      if (!toInt(ms1_min_scan_num, n) || n < 1 || n > 3) {
        LOG_ERROR("MS1 min scan number " << ms1_min_scan_num
                                         << " should be 1, 2, or 3.");
        return false;
      }
      topfd_para_ptr_->setMs1MinScanNum(n);
    }

    if (vm.count("ms2-min-scan-number")) {
      int n = 0;
      if (!toInt(ms2_min_scan_num, n) || n < 1 || n > 3) {
        LOG_ERROR("MS2 min scan number " << ms2_min_scan_num
                                         << " should be 1, 2, or 3.");
        return false;
      }
      topfd_para_ptr_->setMs2MinScanNum(n);
    }

    if (vm.count("pseudo-peak-number")) {
      int n = 0;
      if (!toInt(pseudo_min_peaks, n) || n < 10) {
        LOG_ERROR("Pseudo peak number " << pseudo_min_peaks
                                        << " should be an integer >= 10.");
        return false;
      }
      topdia_para_ptr_->setPseudoMinPeaks(n);
    }

    if (vm.count("ms1-correlation-cutoff")) {
      double cutoff = 0;
      if (!toDouble(ms1_seed_env_inte_corr_tole_cutoff, cutoff) || cutoff < 0 ||
          cutoff > 1) {
        LOG_ERROR("MS1 intensity correlation cutoff "
                  << ms1_seed_env_inte_corr_tole_cutoff
                  << " should be a number in [0,1].");
        return false;
      }
      topdia_para_ptr_->setMs1SeedEnvInteCorrToleCutoff(cutoff);
    }

    if (vm.count("ms2-correlation-cutoff")) {
      double cutoff = 0;
      if (!toDouble(ms2_seed_env_inte_corr_tole_cutoff, cutoff) || cutoff < 0 ||
          cutoff > 1) {
        LOG_ERROR("MS2 intensity correlation cutoff "
                  << ms2_seed_env_inte_corr_tole_cutoff
                  << " should be a number in [0,1].");
        return false;
      }
      topdia_para_ptr_->setMs2SeedEnvInteCorrToleCutoff(cutoff);
    }

    if (vm.count("spectrum-file-name")) {
      spec_file_list_ =
          vm["spectrum-file-name"].as<std::vector<std::string> >();
    }

    if (vm.count("thread-number")) {
      int num = 0;
      if (!toInt(thread_number, num) || num <= 0) {
        LOG_ERROR("Thread number " << thread_number
                                   << " should be a positive integer.");
        return false;
      }
      topfd_para_ptr_->setThreadNum(num);
    }

    if (vm.count("frag-num-filtering")) {
      topfd_para_ptr_->setAANumBasedFilter(true);
    }

    if (vm.count("no-sql")) {
      topfd_para_ptr_->setGeneSql(false);
    }

    if (vm.count("disable-zstd")) {
      topfd_para_ptr_->setZstd(false);
    }
  } catch (std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line " << e.what()
              << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool Argument::validateArguments() {
  if (!std::filesystem::exists(topfd_para_ptr_->getResourceDir())) {
    LOG_ERROR("The directory "
              << topfd_para_ptr_->getResourceDir() << " does not exist!\n"
              << "Please check if the file directory or name contains special "
                 "characters such as spaces or quotation marks.");
    return false;
  }

  for (size_t k = 0; k < spec_file_list_.size(); k++) {
    if (!std::filesystem::exists(spec_file_list_[k])) {
      LOG_ERROR(spec_file_list_[k]
                << " does not exist!\n"
                << "Please check if file directory or name contains special "
                   "characters such as spaces or quotation marks, or the file "
                   "has been deleted.");
      return false;
    }
  }
  int thread_number = topfd_para_ptr_->getThreadNum();
  int valid = mem_check::checkThreadNum(thread_number, "topdia");
  if (!valid) {
    return false;
  }

  // validate activation method
  std::string activation = topfd_para_ptr_->getActivation();
  if (activation != "FILE" && activation != "CID" && activation != "ETD" &&
      activation != "ETCID" && activation != "ETHCD" && activation != "MPD" &&
      activation != "HCD" && activation != "UVPD") {
    LOG_ERROR(
        "Activation method should be one out of "
        "|FILE|CID|ETD|ETCID|ETHCD|HCD|MPD|UVPD.");
    return false;
  }

  return true;
}

}  // namespace toppic
