//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <iostream>
#include <iomanip>

#include "boost/thread/thread.hpp"
#include "common/util/version.hpp"
#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/time_util.hpp"
#include "console/topfd_argument.hpp"

namespace toppic {

Argument::Argument() {
  topfd_para_ptr_ = std::make_shared<TopfdPara>();
}

void Argument::showUsage(boost::program_options::options_description &desc) {
  std::cout << "Usage: toppfd [options] spectrum-file-name" << std::endl;
  std::cout << desc << std::endl;
}

bool Argument::parse(int argc, char* argv[]) {
  std::string max_charge = "";
  std::string max_mass = "";
  std::string mz_error = "";
  std::string ms_two_sn_ratio = "";
  std::string ms_one_sn_ratio = "";
  std::string prec_window = "";
  std::string merged_file_name = "";
  std::string thread_number = "";

  // Define and parse the program options
  try {
    namespace po = boost::program_options;
    po::options_description display_desc("Options");

    display_desc.add_options()
        ("help,h", "Print this help message.")
        ("max-charge,c", po::value<std::string> (&max_charge),
         "<a positive integer>. Set the maximum charge state of precursor and fragment ions. The default value is 30.")
        ("max-mass,m", po::value<std::string> (&max_mass),
         "<a positive number>. Set the maximum monoisotopic mass of precursor and fragment ions. The default value is 100000 Dalton.")
        ("mz-error,e", po::value<std::string> (&mz_error),
         "<a positive number>. Set the error tolerance of m/z values of spectral peaks. The default value is 0.02 m/z.")
        ("ms-one-sn-ratio,r", po::value<std::string> (&ms_one_sn_ratio),
         "<a positive number>. Set the signal/noise ratio for MS1 spectra. The default value is 3.")
        ("ms-two-sn-ratio,t", po::value<std::string> (&ms_two_sn_ratio),
         "<a positive number>. Set the signal/noise ratio for MS/MS spectra. The default value is 1.")
        ("precursor-window,w", po::value<std::string> (&prec_window),
         "<a positive number>. Set the precursor window size. The default value is 3.0 m/z.")
        ("missing-level-one,o","The input spectrum file does not contain MS1 spectra.")
        ("thread-number,u", po::value<std::string> (&thread_number), "<a positive integer>. Number of threads used in the computation. Default value: 1.")
        ("generate-html-folder,g","Generate an html folder containing TopView and spectrum data for visualization.")
        ;

    po::options_description desc("Options");
    desc.add_options() 
        ("help,h", "Print this help message.") 
        ("max-charge,c", po::value<std::string> (&max_charge), "")
        ("max-mass,m", po::value<std::string> (&max_mass), "")
        ("mz-error,e", po::value<std::string> (&mz_error), "")
        ("ms-one-sn-ratio,r", po::value<std::string> (&ms_one_sn_ratio), "")
        ("ms-two-sn-ratio,t", po::value<std::string> (&ms_two_sn_ratio), "")
        ("precursor-window,w", po::value<std::string> (&prec_window), "")
        ("missing-level-one,o", "")
        //("multiple-mass,u", "Output multiple masses for one envelope.")
        ("thread-number,u", po::value<std::string> (&thread_number), "")
        ("generate-html-folder,g","")
        ("keep,k", "Report monoisotopic masses extracted from low quality isotopic envelopes.")
        ("merged-file-name,f", po::value<std::string> (&merged_file_name), 
         "Merge deconvoluted files and specify the name of the merged file.")
        ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(), 
         "Spectrum file name with its path.")
        ;

    po::positional_options_description positional_options;
    positional_options.add("spectrum-file-name", -1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), vm);
      if (vm.count("help")) {
        showUsage(display_desc);
        return false;
      }
      po::notify(vm);
      // throws on error, so do after help in case there are any problems
    }
    catch(boost::program_options::required_option& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }
    catch(boost::program_options::error& e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      showUsage(display_desc);
      return false;
    }

    // get the execution directory
    std::string argv_0(argv[0]);
    std::string exec_dir = file_util::getExecutiveDir(argv_0);

    topfd_para_ptr_->resource_dir_ = file_util::getResourceDir(exec_dir);

    if (vm.count("max-charge")) {
      topfd_para_ptr_->max_charge_ = std::stoi(max_charge);
    }

    if (vm.count("keep")) {
      topfd_para_ptr_->keep_unused_peaks_ = true;
    }

    if (vm.count("max-mass")) {
      topfd_para_ptr_->max_mass_ = std::stod(max_mass);
    }

    if (vm.count("mz-error")) {
      topfd_para_ptr_->mz_error_ = std::stod(mz_error);
    }

    if (vm.count("ms-two-sn-ratio")) {
      topfd_para_ptr_->ms_two_sn_ratio_ = std::stod(ms_two_sn_ratio);
    }

    if (vm.count("ms-one-sn-ratio")) {
      topfd_para_ptr_->ms_one_sn_ratio_ = std::stod(ms_one_sn_ratio);
    }

    if (vm.count("missing-level-one")) {
      topfd_para_ptr_->missing_level_one_ = true;
    }

    if (vm.count("multiple-mass")) {
      topfd_para_ptr_->output_multiple_mass_ = true;
    }

    if (vm.count("precursor-window")) {
      topfd_para_ptr_->prec_window_ = std::stod(prec_window);
    }

    if (vm.count("spectrum-file-name")) {
      spec_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string> >(); 
    }

    if (vm.count("merged-file-name")) {
      if (spec_file_list_.size() > 1) {
        topfd_para_ptr_->merge_files_ = true;
        topfd_para_ptr_->merged_file_name_ = merged_file_name;
      }
    }
    if (vm.count("thread-number")) {
      topfd_para_ptr_->thread_number_ = thread_number;
    }
    if (vm.count("generate-html-folder")) {
      topfd_para_ptr_->gene_html_folder_ = true;
    }
  }
  catch(std::exception& e) {
    std::cerr << "Unhandled Exception in parsing command line "
        << e.what() << ", application will now exit" << std::endl;
    return false;
  }

  return validateArguments();
}

bool Argument::validateArguments() {
  if (!file_util::exists(topfd_para_ptr_->resource_dir_)) {
    LOG_ERROR("Resource direcotry " << topfd_para_ptr_->resource_dir_ << " does not exist!");
    return false;
  }

  for (size_t k = 0; k < spec_file_list_.size(); k++) {
    if (!file_util::exists(spec_file_list_[k])) {
      LOG_ERROR(spec_file_list_[k] << " does not exist!");
      return false;
    }
  }
  std::string thread_number = topfd_para_ptr_->thread_number_;
  try {
    int num = std::stoi(thread_number.c_str());
    if (num <= 0) {
      LOG_ERROR("Thread number " << thread_number << " error! The value should be positive.");
      return false;
    }
    int n = static_cast<int>(boost::thread::hardware_concurrency());
    if(num > n){
      LOG_ERROR("Thread number " << thread_number << " error! The value is too large. Only " << n << " threads are supported.");
      return false;
    }
  } catch (int e) {
    LOG_ERROR("Thread number " << thread_number << " should be a number.");
    return false;
  }
  return true;
}

}  // namespace toppic
