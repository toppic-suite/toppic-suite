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

// topcalib: the two-pass m/z recalibration of toppic_test experiment 152
// (ms/calib). "fit" fits the MS1 and MS/MS error models on the confident
// PrSMs of a zero-shift TopPIC search; "apply" corrects the precursor and
// fragment masses of a TopFD MS/MS msalign file with a model (no TopFD
// rerun). The C++ form of the experiment's calib.py and recal_msalign.py.

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "common/base/base_data.hpp"
#include "common/util/file_util.hpp"
#include "common/util/version.hpp"
#include "ms/calib/mass_calib_fit.hpp"
#include "ms/calib/mass_calib_model.hpp"
#include "ms/calib/msalign_mass_calib.hpp"

namespace {

namespace po = boost::program_options;

void printUsage(std::ostream& out) {
  out << "Usage: topcalib fit [options] <PrSM table> <MS/MS msalign file> "
         "<model file>"
      << std::endl
      << "       topcalib apply <model file> <MS/MS msalign file> "
         "<output msalign file>"
      << std::endl
      << std::endl
      << "fit:   fits the m/z error models (ppm) on the PrSMs of a zero-shift "
         "TopPIC"
      << std::endl
      << "       search (a *_toppic_prsm_single.tsv file, e.g. of toppic -s 0 "
         "-D) and"
      << std::endl
      << "       the msalign file it searched, and writes them as JSON:"
      << std::endl
      << "         MS1:   e = c0 + c_rt (rt - rt0) / 60" << std::endl
      << "         MS/MS: e = c0 + c_mz (m/z - 800) / 1000 + c_rt (rt - rt0) "
         "/ 60"
      << std::endl
      << "       (rt in minutes; MS1 from the precursor errors, MS/MS from "
         "the b / y"
      << std::endl
      << "       fragment errors)." << std::endl
      << "apply: corrects the precursor and fragment masses of a TopFD MS/MS "
         "msalign"
      << std::endl
      << "       file with a model: m/z' = m/z (1 - e 1e-6)." << std::endl;
}

bool parseOptions(int argc, char* argv[], const po::options_description& desc,
                  int pos_num, po::variables_map* vm) {
  po::options_description all;
  all.add(desc);
  all.add_options()("files", po::value<std::vector<std::string>>(), "");
  po::positional_options_description pos;
  pos.add("files", -1);
  try {
    po::store(po::command_line_parser(argc, argv)
                  .options(all)
                  .positional(pos)
                  .run(),
              *vm);
    po::notify(*vm);
  } catch (const po::error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return false;
  }
  if (vm->count("help")) return true;
  int n = vm->count("files")
              ? static_cast<int>((*vm)["files"].as<std::vector<std::string>>().size())
              : 0;
  if (n != pos_num) {
    std::cerr << "Error: " << pos_num << " file names expected, " << n
              << " given." << std::endl;
    return false;
  }
  return true;
}

void checkInput(const std::string& file_name) {
  if (!std::filesystem::exists(file_name)) {
    std::cerr << "Error: " << file_name << " does not exist." << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void printLevel(const char* name, const toppic::MassCalibLevel& c) {
  std::cout << "  " << name << ": c0 " << c.c0 << " ppm, c_mz " << c.c_mz
            << " ppm per 1000 m/z, c_rt " << c.c_rt_per_hour << " ppm per hour"
            << std::endl;
}

int runFit(int argc, char* argv[], const std::string& argv_0) {
  double max_evalue = 1e-5;
  double tol_ppm = 20.0;
  po::options_description desc("Options of topcalib fit");
  desc.add_options()("help,h", "Print the help message.")(
      "max-evalue,e", po::value<double>(&max_evalue),
      "<positive number>. The largest E-value of a calibrant PrSM (no "
      "unexpected shift and no variable PTM are required too). Default: "
      "1e-5.")(
      "tolerance,t", po::value<double>(&tol_ppm),
      "<positive number>. The tolerance (ppm) of the fragment matching. "
      "Default: 20.");
  po::variables_map vm;
  if (!parseOptions(argc, argv, desc, 3, &vm)) {
    printUsage(std::cerr);
    return EXIT_FAILURE;
  }
  if (vm.count("help")) {
    printUsage(std::cout);
    std::cout << std::endl << desc << std::endl;
    return EXIT_SUCCESS;
  }
  if (!(max_evalue > 0) || !(tol_ppm > 0)) {
    std::cerr << "Error: --max-evalue and --tolerance must be positive."
              << std::endl;
    return EXIT_FAILURE;
  }
  std::vector<std::string> files = vm["files"].as<std::vector<std::string>>();
  checkInput(files[0]);
  checkInput(files[1]);
  std::string resource_dir =
      toppic::file_util::getResourceDir(toppic::file_util::getExecutiveDir(argv_0));
  toppic::base_data::init(resource_dir);

  std::cout << "Reading the calibrants - started." << std::endl;
  toppic::MassCalibrantSet set =
      toppic::mass_calib_fit::readCalibrants(files[0], files[1], max_evalue, tol_ppm);
  std::cout << "Reading the calibrants - finished: " << set.used_prsm_num
            << " of " << set.prsm_num << " PrSMs with E-value <= "
            << max_evalue << ", no unexpected shift and no variable PTM; "
            << set.prec.size() << " precursor and " << set.frag.size()
            << " fragment calibrants." << std::endl;
  std::cout << "Fitting the models - started." << std::endl;
  toppic::MassCalibModel model = toppic::mass_calib_fit::fitModel(set);
  std::cout << "Fitting the models - finished." << std::endl;
  std::cout << std::setprecision(4);
  printLevel("MS1", model.ms1);
  printLevel("MS/MS", model.ms2);
  std::cout << "  rt0 " << model.rt0_min << " min; robust SD of the residuals: "
            << "precursors " << model.prec_robust_sd << " ppm, fragments "
            << model.frag_robust_sd << " ppm (" << model.frag_prsm_num
            << " PrSMs)" << std::endl;
  model.write(files[2]);
  std::cout << "Model written to " << files[2] << "." << std::endl;
  toppic::base_data::release();
  return EXIT_SUCCESS;
}

int runApply(int argc, char* argv[]) {
  po::options_description desc("Options of topcalib apply");
  desc.add_options()("help,h", "Print the help message.");
  po::variables_map vm;
  if (!parseOptions(argc, argv, desc, 3, &vm)) {
    printUsage(std::cerr);
    return EXIT_FAILURE;
  }
  if (vm.count("help")) {
    printUsage(std::cout);
    std::cout << std::endl << desc << std::endl;
    return EXIT_SUCCESS;
  }
  std::vector<std::string> files = vm["files"].as<std::vector<std::string>>();
  checkInput(files[0]);
  checkInput(files[1]);
  if (std::filesystem::exists(files[2]) &&
      std::filesystem::equivalent(files[1], files[2])) {
    std::cerr << "Error: the output file is the input file." << std::endl;
    return EXIT_FAILURE;
  }
  toppic::MassCalibModel model = toppic::MassCalibModel::read(files[0]);
  std::cout << "Correcting " << files[1] << " - started." << std::endl;
  int sp_num = toppic::msalign_mass_calib::process(model, files[1], files[2]);
  std::cout << "Correcting - finished: " << sp_num << " spectra written to "
            << files[2] << "." << std::endl;
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char* argv[]) {
  std::cout << "topcalib " << toppic::Version::getVersion() << std::endl;
  std::string cmd = argc > 1 ? argv[1] : "";
  if (cmd == "fit") return runFit(argc - 1, argv + 1, argv[0]);
  if (cmd == "apply") return runApply(argc - 1, argv + 1);
  if (cmd == "-h" || cmd == "--help") {
    printUsage(std::cout);
    return EXIT_SUCCESS;
  }
  printUsage(std::cerr);
  return EXIT_FAILURE;
}
