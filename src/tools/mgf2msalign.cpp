//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include <pwiz/data/msdata/DefaultReaderList.hpp>
#include <pwiz/data/msdata/MSDataFile.hpp>
#include <pwiz/data/msdata/SpectrumInfo.hpp>

#include "base/activation_base.hpp"
#include "spec/ms_header.hpp"

namespace pm = pwiz::msdata;

int main(int argc, const char *argv[]) {
  if (argc < 2) {
    std::cout << "Please provide the mgf file." << std::endl;
    std::cout << argv[0] << " <MGF FILE>" << std::endl;
    return EXIT_SUCCESS;
  }

  std::string file_name(argv[1]);

  std::ofstream ms1_out(file_name.substr(0, file_name.length() - 4) + "_ms1.msalign");

  std::ofstream ms2_out(file_name.substr(0, file_name.length() - 4) + "_ms2.msalign");

  std::shared_ptr<pm::MSDataFile> mds_ptr
      = std::make_shared<pm::MSDataFile>(argv[1]);

  pm::SpectrumListPtr spec_list_ptr = mds_ptr->run.spectrumListPtr;

  size_t n = spec_list_ptr->size();

  std::cout << "A total of " << n << " spectra in " << argv[1] << std::endl;

  int ms1_cnt = 0, ms2_cnt = 0;

  for (size_t i = 0; i < n; i++) {
    pm::SpectrumPtr spec = spec_list_ptr->spectrum(i, true);
    std::vector<pm::MZIntensityPair> pairs;
    spec->getMZIntensityPairs(pairs);
    pm::SpectrumInfo spec_info(*spec);
    int ms_level = spec_info.msLevel;
    if (ms_level == 2) {
      double prec_mz;
      if (spec_info.precursors.size() == 0) {
        prec_mz = 0;
      } else {
        prec_mz = spec_info.precursors[0].mz;
      }

      if (prec_mz < 0) {
        prec_mz = 0;
      }

      int prec_charge;
      if (spec_info.precursors.size() == 0) {
        prec_charge = 1;
      } else {
        prec_charge = static_cast<int>(spec_info.precursors[0].charge);
      }

      if (prec_charge  < 0) {
        prec_charge = 1;
      }

      prot::MsHeaderPtr header_ptr = std::make_shared<prot::MsHeader>();
      header_ptr->setId(ms2_cnt);
      ms2_cnt++;
      header_ptr->setScan(spec_info.scanNumber);
      header_ptr->setMsLevel(ms_level);
      header_ptr->setPrecCharge(prec_charge);
      header_ptr->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
      header_ptr->setPrecMonoMz(prec_mz);
      header_ptr->setRetentionTime(spec_info.retentionTime);

      std::string ac_name;
      if (spec->precursors.size() != 0) {
        std::vector<pwiz::data::CVParam> cv_list = spec->precursors[0].activation.cvParams;
        for (size_t i = 0; i < cv_list.size(); i++) {
          if (cv_list[i].cvid == pwiz::cv::MS_CID) {
            ac_name = "CID";
            break;
          } else if (cv_list[i].cvid == pwiz::cv::MS_HCD) {
            ac_name = "HCD";
            break;
          } else if (cv_list[i].cvid == pwiz::cv::MS_ETD) {
            ac_name = "ETD";
            break;
          }
        }
      }
      if (ac_name == "") {
        std::cerr << "No activation information is available. Using HCD." << std::endl;
      }
      prot::ActivationPtr activation_ptr = prot::ActivationBase::getActivationPtrByName(ac_name);
      header_ptr->setActivationPtr(activation_ptr);

      ms2_out << "BEGIN IONS" << std::endl;
      ms2_out << "ID=" << header_ptr->getId() << std::endl;
      ms2_out << "SCANS=" << header_ptr->getScansString() << std::endl;
      ms2_out << "RETENTION_TIME=" << std::setprecision(2)
          << header_ptr->getRetentionTime() << std::endl;
      if (header_ptr->getActivationPtr() != nullptr) {
        ms2_out << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
      }
      ms2_out << "MS_ONE_ID=" << header_ptr->getMsOneId() << std::endl;
      ms2_out << "MS_ONE_SCAN=" << header_ptr->getMsOneScan() << std::endl;
      ms2_out << "PRECURSOR_MZ=" << std::setprecision(5)
          << header_ptr->getPrecMonoMz() << std::endl;
      ms2_out << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
      ms2_out << "PRECURSOR_MASS=" << std::setprecision(5)
          << header_ptr->getPrecMonoMass() << std::endl;

      for (size_t k = 0; k < pairs.size(); k++) {
        if (pairs[k].intensity > 0)
          ms2_out << pairs[k].mz << "\t" << pairs[k].intensity << "\t1" << std::endl;
      }

      ms2_out << "END IONS" << std::endl << std::endl;
    } else {
      prot::MsHeaderPtr header_ptr = std::make_shared<prot::MsHeader>();
      header_ptr->setId(ms1_cnt);
      ms1_cnt++;
      header_ptr->setScan(spec_info.scanNumber);
      header_ptr->setMsLevel(ms_level);
      header_ptr->setPrecCharge(0);
      header_ptr->setTitle("Scan_" + std::to_string(spec_info.scanNumber));
      header_ptr->setRetentionTime(spec_info.retentionTime);

      ms1_out << "BEGIN IONS" << std::endl;
      ms1_out << "ID=" << header_ptr->getId() << std::endl;
      ms1_out << "SCANS=" << header_ptr->getScansString() << std::endl;
      ms1_out << "RETENTION_TIME=" << std::setprecision(2)
          << header_ptr->getRetentionTime() << std::endl;
      if (header_ptr->getActivationPtr() != nullptr) {
        ms1_out << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
      }

      for (size_t k = 0; k < pairs.size(); k++) {
        if (pairs[k].intensity > 0)
          ms1_out << pairs[k].mz << "\t" << pairs[k].intensity << "\t1" << std::endl;
      }

      ms1_out << "END IONS" << std::endl << std::endl;
    }
  }

  ms1_out.close();

  ms2_out.close();
}
