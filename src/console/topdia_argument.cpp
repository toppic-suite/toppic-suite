//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
#include "common/util/mem_check.hpp"
#include "console/topdia_argument.hpp"

namespace toppic {

    Argument::Argument() {
        topdia_para_ptr_ = std::make_shared<TopdiaPara>();
    }

    void Argument::showUsage(boost::program_options::options_description &desc) {
        std::cout << "Usage: topdia [options] spectrum-file-name" << std::endl;
        std::cout << desc << std::endl;
        std::cout << "Version: " << Version::getVersion() << std::endl;
    }

    bool Argument::parse(int argc, char* argv[]) {
        std::string max_charge = "";
        std::string max_mass = "";
        std::string mz_error = "";
        std::string ms_two_sn_ratio = "";
        std::string ms_one_sn_ratio = "";
        std::string prec_window = "";
        std::string thread_number = "";
        std::string activation = "";
        std::string ecscore_cutoff = "";
        std::string min_scan_num = "";

        // Define and parse the program options
        try {
            namespace po = boost::program_options;
            po::options_description display_desc("Options");

            display_desc.add_options()
                    ("help,h", "Print this help message.")
                    ("activation,a", po::value<std::string> (&activation),
                     "<CID|ETD|HCD|MPD|UVPD|FILE>. Fragmentation method of MS/MS spectra. When FILE is used, the fragmentation methods of spectra are given in the input spectral data file. Default value: FILE.")
                    ("max-charge,c", po::value<std::string> (&max_charge),
                     "<a positive integer>. Set the maximum charge state of precursor and fragment ions. The default value is 30.")
                    ("max-mass,m", po::value<std::string> (&max_mass),
                     "<a positive number>. Set the maximum monoisotopic mass of precursor and fragment ions. The default value is 70,000 Dalton.")
                    ("mz-error,e", po::value<std::string> (&mz_error),
                     "<a positive number>. Set the error tolerance of m/z values of spectral peaks. The default value is 0.02 m/z.")
                    ("ms-one-sn-ratio,r", po::value<std::string> (&ms_one_sn_ratio),
                     "<a positive number>. Set the signal-to-noise ratio for MS1 spectra. The default value is 3.")
                    ("ms-two-sn-ratio,s", po::value<std::string> (&ms_two_sn_ratio),
                     "<a positive number>. Set the signal-to-noise ratio for MS/MS spectra. The default value is 1.")
                    ("missing-level-one,o","MS1 spectra are missing in the input file.")
                    ("msdeconv,n", "Use the MS-Deconv score to rank isotopic envelopes.")
                    ("precursor-window,w", po::value<std::string> (&prec_window),
                     "<a positive number>. Set the default precursor window size. The default value is 3.0 m/z. When the input file contains the information of precursor windows, the parameter will be ignored.")
                    ("ecscore-cutoff,t", po::value<std::string> (&ecscore_cutoff),
                     "<a positive number in [0,1]>. Set the ECScore cutoff value for proteoform features. The default value is 0.5.")
                    ("min-scan-number,b",po::value<std::string> (&min_scan_num),
                     "<1|2|3>. The minimum number of MS1 scans in which a proteoform feature is detected. The default value is 3.")
                    ("single-scan-noise,i","Use the peak intensity noise levels in single MS1 scans to filter out low intensity peaks in proteoform feature detection. The default method is to use the peak intensity noise level of the whole LC-MS map to filter out low intensity peaks.")
                    ("additional-feature-search,f","Perform additional feature search for MS/MS scans that do not have detected proteoform features in their precursor isolation windows. In additional search, the signal noise ratio is set to 0, the min scan number is set to 1, and the ecscore cutoff is set to 0.")
                    ("disable-final-filtering,d","Skip the final filtering of envelopes in MS/MS scans.")
                    ("thread-number,u", po::value<std::string> (&thread_number), "<a positive integer>. Number of threads used in spectral deconvolution. Default value: 1.")
                    ("skip-html-folder,g","Skip the generation of HTML files for visualization.")
                    ;

            po::options_description desc("Options");
            desc.add_options()
                    ("help,h", "Print this help message.")
                    ("activation,a", po::value<std::string> (&activation), "")
                    ("max-charge,c", po::value<std::string> (&max_charge), "")
                    ("max-mass,m", po::value<std::string> (&max_mass), "")
                    ("mz-error,e", po::value<std::string> (&mz_error), "")
                    ("ms-one-sn-ratio,r", po::value<std::string> (&ms_one_sn_ratio), "")
                    ("ms-two-sn-ratio,s", po::value<std::string> (&ms_two_sn_ratio), "")
                    ("precursor-window,w", po::value<std::string> (&prec_window), "")
                    ("ecscore-cutoff,t", po::value<std::string> (&ecscore_cutoff), "")
                    ("missing-level-one,o", "")
                    ("single-scan-noise,i","")
                    ("additional-feature-search,f","")
                    ("min-scan-number,b",po::value<std::string> (&min_scan_num),"")
                    ("thread-number,u", po::value<std::string> (&thread_number), "")
                    ("skip-html-folder,g","")
                    ("msdeconv,n", "")
                    ("disable-final-filtering,d", "")
                    ("keep,k", "Report monoisotopic masses extracted from low quality isotopic envelopes.")
                    ("spectrum-file-name", po::value<std::vector<std::string> >()->multitoken()->required(),
                     "Spectrum file name with its path.")
                    ;

            po::positional_options_description positional_options;
            positional_options.add("spectrum-file-name", -1);

            po::variables_map vm;
            try {
                po::store(
                        po::command_line_parser(argc, argv).options(desc).positional(positional_options).run(), vm);
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
            if (file_util::checkSpace(exec_dir)) {
                LOG_ERROR("Current directory " << exec_dir << " contains space and will cause errors in the program!")
                exit(EXIT_FAILURE);
            }

            topdia_para_ptr_->setResourceDir(file_util::getResourceDir(exec_dir));

            if (vm.count("activation")) {
                topdia_para_ptr_->setActivation(activation);
            }

            if (vm.count("max-charge")) {
                topdia_para_ptr_->setMaxCharge(std::stoi(max_charge));
            }

            if (vm.count("keep")) {
                topdia_para_ptr_->setKeepUnusedPeaks(true);
            }

            if (vm.count("msdeconv")) {
                topdia_para_ptr_->setUseMsDeconv(true);
            }

            if (vm.count("max-mass")) {
                topdia_para_ptr_->setMaxMass(std::stod(max_mass));
            }

            if (vm.count("mz-error")) {
                topdia_para_ptr_->setMzError(std::stod(mz_error));
            }

            if (vm.count("ms-two-sn-ratio")) {
                topdia_para_ptr_->setMsTwoSnRatio(std::stod(ms_two_sn_ratio));
            }

            if (vm.count("ms-one-sn-ratio")) {
                topdia_para_ptr_->setMsOneSnRatio(std::stod(ms_one_sn_ratio));
            }

            if (vm.count("missing-level-one")) {
                topdia_para_ptr_->setMissingLevelOne(true);
            }

            if (vm.count("multiple-mass")) {
                topdia_para_ptr_->setOutputMultipleMass(true);
            }

            if (vm.count("precursor-window")) {
                topdia_para_ptr_->setPrecWindowWidth(std::stod(prec_window));
            }

            if (vm.count("ecscore-cutoff")) {
                topdia_para_ptr_->setEcscoreCutoff(std::stod(ecscore_cutoff));
            }

            if (vm.count("single-scan-noise")) {
                topdia_para_ptr_->setUseSingleScanNoiseLevel(true);
            }

            if (vm.count("additional-feature-search")) {
                topdia_para_ptr_->setSearchPrecWindow(true);
            }

            if (vm.count("min-scan-number")) {
                try {
                    int n = std::stoi(min_scan_num);
                    if (n < 1 || n > 3) {
                        LOG_ERROR("Min scan number " << min_scan_num << " should be 1, 2, or 3.");
                        return false;
                    }
                    topdia_para_ptr_->setMinScanNum(n);
                } catch (std::exception& e) {
                    LOG_ERROR("Min scan number " << min_scan_num << " should be 1, 2, or 3.");
                    return false;
                }
            }

            if (vm.count("spectrum-file-name")) {
                spec_file_list_ = vm["spectrum-file-name"].as<std::vector<std::string> >();
            }

            if (vm.count("thread-number")) {
                try {
                    topdia_para_ptr_->setThreadNum(std::stoi(thread_number));
                } catch (std::exception& e) {
                    LOG_ERROR("Thread number " << thread_number << " should be a number.");
                    return false;
                }
            }
            if (vm.count("skip-html-folder")) {
                topdia_para_ptr_->setGeneHtmlFolder(false);
            }
            if (vm.count("disable-final-filtering")) {
                topdia_para_ptr_->setDoFinalFiltering(false);
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
        if (!file_util::exists(topdia_para_ptr_->getResourceDir())) {
            LOG_ERROR("The directory " << topdia_para_ptr_->getResourceDir() << " does not exist!\n"
                                       << "Please check if the file directory or name contains special characters such as spaces or quotation marks.");
            return false;
        }

        for (size_t k = 0; k < spec_file_list_.size(); k++) {
            if (!file_util::exists(spec_file_list_[k])) {
                LOG_ERROR(spec_file_list_[k] << " does not exist!\n"
                                             << "Please check if file directory or name contains special characters such as spaces or quotation marks, or the file has been deleted.");
                return false;
            }
        }
        int thread_number = topdia_para_ptr_->getThreadNum();
        int valid = mem_check::checkThreadNum(thread_number, "topdia");
        if (!valid) {
            return false;
        }

        //validate activation method
        std::string activation = topdia_para_ptr_->getActivation();
        if (activation != "FILE" && activation != "CID" && activation != "ETD"
            && activation != "MPD" && activation != "HCD" && activation != "UVPD"){
            //throw InvalidActivation();
            LOG_ERROR("Activation method should be one out of |FILE|CID|ETD|HCD|MPD|UVPD.");
            return false;
        }

        return true;
    }

}  // namespace toppic
