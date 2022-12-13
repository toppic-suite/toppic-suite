
#include "src/topfd/feature_detect/feature/feature_map.hpp"

namespace toppic {
  namespace feature_map {
    void assign_features2() {
      std::cout << "Yapooooooozzzzzaaaa !!!!!!!!!!!!" << std::endl;
    }
  }
}

//    inline double get_mz(double mass, int charge) {
//      double proton = 1.00727;
//      return (mass + (charge * proton)) / charge;
//    }

//    inline void readHeaders(const std::string &file_name, MsHeaderPtrVec &header_ptr_vec) {
//      SimpleMsAlignReader sp_reader(file_name);
//      DeconvMsPtr ms_ptr;
//      while ((ms_ptr = sp_reader.getNextMsPtr()) != nullptr) {
//        header_ptr_vec.push_back(ms_ptr->getMsHeaderPtr());
//      }
//    }

//    inline bool isMatch(FracFeaturePtr feature_ptr, MsHeaderPtr header, FeatureParaPtr para_ptr) {
//      int ms1_scan = header->getMsOneScan();
//      if (ms1_scan < feature_ptr->getScanBegin()) {
//        return false;
//      }
//      if (ms1_scan > feature_ptr->getScanEnd()) {
//        return false;
//      }
//      double prec_mass = header->getPrecMonoMass();
//      std::vector<double> search_masses = para_ptr->getSearchMasses(prec_mass);
//
//      double feature_mass = feature_ptr->getMonoMass();
//
//      double min_diff = std::numeric_limits<double>::max();
//      for (size_t j = 0; j < search_masses.size(); j++) {
//        double mass_diff = std::abs(search_masses[j] - feature_mass);
//        if (mass_diff < min_diff) {
//          min_diff = mass_diff;
//        }
//      }
//      double error_tole = para_ptr->peak_tolerance_ptr_->compStrictErrorTole(prec_mass);
//      if (min_diff <= error_tole) {
//        return true;
//      }
//      return false;
//    }

//    inline FracFeaturePtr getMatchedFeaturePtr(FracFeaturePtrVec &features, MsHeaderPtr header, FeatureParaPtr para_ptr) {
//      for (size_t i = 0; i < features.size(); i++) {
//        FracFeaturePtr ft_ptr = features[i];
//        if (isMatch(features[i], header, para_ptr)) {
//          return features[i];
//        }
//      }
//      return nullptr;
//    }

//    void assign_features() {
//      std::cout << "Yapooooooozzzzzaaaa !!!!!!!!!!!!" << std::endl;
//    }
//      MsHeaderPtrVec header_ptr_vec;
//      readHeaders(ms2_file_name, header_ptr_vec);
//      int isolation_windows_mz = 10;
//      int isolation_windows_num = int(header_ptr_vec.size() / ms1_ptr_vec.size()) + 1;
//      std::cout << "# of MS1 scans: " << ms1_ptr_vec.size() << ", # of MS2 scans: " << header_ptr_vec.size()
//                << " and num of isolation windows: " << isolation_windows_num << std::endl;
//
//      std::vector<std::vector<double>> feature_mapping(ms1_ptr_vec.size(),
//                                                       std::vector<double>(isolation_windows_num, -1));
//      std::vector<std::vector<double>> base_mzs(ms1_ptr_vec.size(), std::vector<double>(isolation_windows_num, 0));
//      for (int spec_id = 0; spec_id < ms1_ptr_vec.size(); spec_id++) {
//        for (int iso_win = 0; iso_win < isolation_windows_num; iso_win++) {
//          double base_mz = precMzs[(spec_id * isolation_windows_num) + iso_win];
//          double base_min_mz = base_mz - (isolation_windows_mz / 2);
//          double base_max_mz = base_mz + (isolation_windows_mz / 2);
//          base_mzs[spec_id][iso_win] = base_mz;
//          for (int feat_id = 0; feat_id < features.size(); feat_id++) {
//            FracFeaturePtr f = features[feat_id];
//            if (spec_id >= f->getMinMs1Id() and spec_id <= f->getMaxMs1Id()) {
//              SingleChargeFeaturePtrVec sfs = f->getSingleFeatures();
//              for (int sf_id = 0; sf_id < sfs.size(); sf_id++) {
//                SingleChargeFeaturePtr sf = sfs[sf_id];
//                double mz = get_mz(f->getMonoMass(), sf->getCharge());
//                if (mz >= base_min_mz and mz < base_max_mz) {
//                  feature_mapping[spec_id][iso_win] = sf->getIntensity();
//                  std::cout << "(" << spec_id << ", " << iso_win << ") ---- ";
//                  std::cout << f->getId() << ", " << f->getMinMs1Id() << ", " << f->getMaxMs1Id() << ", "
//                            << f->getMonoMass() << ", " << mz << std::endl;
//                  break;
//                }
//              }
//            }
//          }
//        }
//      }
//      write_out_files::write_feature_map(feature_mapping);
//      write_out_files::write_base_msz(base_mzs);
//
//      MsHeaderPtr hh = header_ptr_vec[0];
//      std::cout << "Scan # " << hh->getMsOneScan() << ", Spec # " << hh->getMsOneId() << ", Mono Mass "
//                << hh->getPrecMonoMass() << ", Mono Mz "
//                << hh->getPrecMonoMz() << ", Charge " << hh->getPrecCharge() << ", Spec MZ " << hh->getPrecSpMz()
//                << ", Target MZ " << hh->getPrecTargetMz()
//                << ", Win Begin " << hh->getPrecWinBegin() << ", Win End " << hh->getPrecWinEnd() << ", RT "
//                << hh->getRetentionTime() << std::endl;
//
//      MsHeaderPtr hh1 = ms1_ptr_vec[0]->getMsHeaderPtr();
//      std::cout << "Scan # " << hh1->getFirstScanNum() << ", Spec # " << hh1->getId() << ", Mono Mass "
//                << hh1->getPrecMonoMass() << ", Mono Mz " << hh1->getPrecMonoMz() << ", Charge " << hh1->getPrecCharge()
//                << ", Spec MZ " << hh1->getPrecSpMz() << ", Target MZ " << hh1->getPrecTargetMz() << ", Win Begin " <<
//                hh1->getPrecWinBegin() << ", Win End " << hh1->getPrecWinEnd() << ", RT " << hh1->getRetentionTime()
//                << std::endl;
//    }

//  }
//}
