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


#include "generate_pseudo_spectrum.hpp"

namespace toppic {
    GeneratePseudoSpectrum::GeneratePseudoSpectrum(const TopdiaParaPtr &topdia_para_ptr) {
        std::string output_base_name = topdia_para_ptr->getOutputBaseName();

        DeconvMsPtrVec deconv_ms1_ptr_vec;
        std::string ms1_file_name = output_base_name + "_ms1.msalign";
        msalign_reader_util::readAllSpectra(ms1_file_name, deconv_ms1_ptr_vec);
        for (const auto &ms1_data: deconv_ms1_ptr_vec)
            rt_ms1_.push_back(ms1_data->getMsHeaderPtr()->getRetentionTime() / 60);

        DeconvMsPtrVec deconv_ms2_ptr_vec;
        std::string ms2_file_name = output_base_name + "_ms2.msalign";
        msalign_reader_util::readAllSpectra(ms2_file_name, deconv_ms2_ptr_vec);

        // get isolation window base mz
        std::set<double> base_mz_set;
        for (auto &ms2_data: deconv_ms2_ptr_vec) {
            if (ms2_data->getMsHeaderPtr()->getMsLevel() == 1)
                continue;
            double base_mz = (ms2_data->getMsHeaderPtr()->getPrecWinBegin() + ms2_data->getMsHeaderPtr()->getPrecWinEnd())/2;
            base_mz_set.insert(base_mz);
        }
        isolation_window_base_mz_ = std::vector<double>(base_mz_set.begin(), base_mz_set.end());

        // get retention times
        for (auto base_mz: isolation_window_base_mz_) {
            DeconvMsPtrVec deconv_ms2_ptr_shortlisted_vec;
            std::vector<double> rt_ms2_window;
            for (auto &ms2_data: deconv_ms2_ptr_vec) {
                double win_mz = (ms2_data->getMsHeaderPtr()->getPrecWinBegin() + ms2_data->getMsHeaderPtr()->getPrecWinEnd())/2;
                if (win_mz != base_mz)
                    continue;
                deconv_ms2_ptr_shortlisted_vec.push_back(ms2_data);
                rt_ms2_window.push_back(ms2_data->getMsHeaderPtr()->getRetentionTime() / 60);
            }
            if (rt_ms2_window.size() < rt_ms1_.size())
                rt_ms2_window.push_back(rt_ms2_window.back());
            rt_ms2_.push_back(rt_ms2_window);
        }

        // read feature files
        std::string filename = output_base_name + "_frac_ms1.mzrt.csv";
        ms1_features_ = mzrtFeature::read_record(filename);
        for (auto base_mz: isolation_window_base_mz_) {
            filename = output_base_name + "_" + std::to_string(int(base_mz - topdia_para_ptr->getPrecWindowWidth() / 2)) +
                       "_" + std::to_string(int(base_mz + topdia_para_ptr->getPrecWindowWidth() / 2)) + "_frac_ms2.mzrt.csv";
            mzrtFeaturePtrVec feature_list = mzrtFeature::read_record(filename);
            ms2_features_.push_back(feature_list);
        }

        // get interpolated xic
        double max_rt = get_max_rt();
        double interval = 0.01;
        std::vector<double> rt_target;
        for (double i = 0; i <= max_rt; i = i + interval) {
            rt_target.push_back(i);
        }

        for (auto &ms1_feature: ms1_features_)
            ms1_feature->setInterpolatedXic(interp(rt_target, rt_ms1_, ms1_feature->getXic()));

        for (int iso_win_idx = 0; iso_win_idx < static_cast<int>(isolation_window_base_mz_.size()); iso_win_idx++)
            for (auto &ms2_feature: ms2_features_[iso_win_idx])
                ms2_feature->setInterpolatedXic(interp(rt_target, rt_ms2_[iso_win_idx], ms2_feature->getXic()));
    }

    void GeneratePseudoSpectrum::process(const TopdiaParaPtr &topdia_para_ptr) {
        int feature_id = 0;
        EnvParaPtr env_para_ptr = std::make_shared<EnvPara>(topdia_para_ptr->getMzError());
        for (int iso_win_idx = 0; iso_win_idx < static_cast<int>(isolation_window_base_mz_.size()); iso_win_idx++) {
            PseudoSpectrumPtrVec pseudo_spectra;
            double base_mz = isolation_window_base_mz_[iso_win_idx];
            mzrtFeaturePtrVec selected_ms1_features = get_iso_win_ms1_features(iso_win_idx, topdia_para_ptr->getPrecWindowWidth());
            std::sort(selected_ms1_features.begin(), selected_ms1_features.end(), compareFeaturesInte);
            std::cout << "Processing Isolation window " << iso_win_idx << " with " << selected_ms1_features.size()
                      << " features." << std::endl;
            for (auto &ms1_feature: selected_ms1_features) {
                int apex_cycle_distance_tole = std::min(3, ms1_feature->getCycleSpan() / 2);
                ms1_feature->setBaseMz(base_mz);
                std::vector<PseudoPeaks> pseudo_peak_list;
                for (int ms2_feature_idx = 0; ms2_feature_idx < static_cast<int>(ms2_features_[iso_win_idx].size()); ms2_feature_idx++) {
                    if (ms2_features_[iso_win_idx][ms2_feature_idx]->getUsedStatus()) continue;
                    mzrtFeaturePtr ms2_feature = ms2_features_[iso_win_idx][ms2_feature_idx];
                    if (ms2_feature->getMass() < ms1_feature->getMass()) {

                        int apex_cycle_distance = get_apex_cycle_distance(ms1_feature, ms2_feature);
                        if (apex_cycle_distance > apex_cycle_distance_tole) continue;

                        /// Get shared intensity
                        std::vector<double> feature_xic_ms1 = ms1_feature->getInterpolatedXic();
                        std::vector<double> feature_xic_ms2 = ms2_feature->getInterpolatedXic();
                        double shared_area = computeSharedArea(mzrtFeature::normalizeXIC(feature_xic_ms1),
                                                               mzrtFeature::normalizeXIC(feature_xic_ms2));

                        PseudoPeaks peak = PseudoPeaks(ms2_feature->getMass(), ms2_feature->getMonoMz(), ms2_feature->getCharge(),
                                                       ms2_feature->getIntensity(), 0, 0, shared_area, ms2_feature->getCycleSpan(),
                                                       apex_cycle_distance, ms2_feature->getTimeBegin(), ms2_feature->getTimeEnd(),
                                                       ms2_feature->getApexCycle());
                        peak.setMs2FeatureIdx(ms2_feature_idx);
                        pseudo_peak_list.push_back(peak);
                    }
                }
                /// filter peaks based on low-high mass dividers
                score_pseudo_peaks(pseudo_peak_list, ms1_feature);
                std::vector<PseudoPeaks> filtered_pseudo_peak_list = filterPseudoPeaks(env_para_ptr, ms1_feature,
                                                                                       ms2_features_[iso_win_idx],
                                                                                       pseudo_peak_list, topdia_para_ptr->getPseudoScoreCutoff(), topdia_para_ptr->getPseudoMinPeaks());
                ms1_feature->setPseudoPeakNum(static_cast<int>(filtered_pseudo_peak_list.size()));
                PseudoSpectrumPtr pseudo_spec_ptr = std::make_shared<PseudoSpectrum>(ms1_feature, filtered_pseudo_peak_list);
                pseudo_spectra.push_back(pseudo_spec_ptr);
            }

            for (const auto &spec: pseudo_spectra) {
                write_pseudo_spectrum(topdia_para_ptr, feature_id, spec->getMs1Feature(), spec->getFragmentFeatures());
                feature_id++;
            }
        }
    }

    double GeneratePseudoSpectrum::get_max_rt() {
        double max_rt = rt_ms1_.back();
        for (int iso_win_idx = 0; iso_win_idx < static_cast<int>(isolation_window_base_mz_.size()); iso_win_idx++) {
            if (max_rt > rt_ms2_[iso_win_idx].back())
                max_rt = rt_ms2_[iso_win_idx].back();
        }
        return max_rt;
    }

    std::vector<double> GeneratePseudoSpectrum::interp(const std::vector<double> &x,
                                                       const std::vector<double> &xp, const std::vector<double> &fp) {
        std::vector<double> interpolatedValues;
        interpolatedValues.reserve(x.size());

        for (double xi: x) {
            auto it = std::lower_bound(xp.begin(), xp.end(), xi);
            if (it == xp.begin() || it == xp.end()) {
                interpolatedValues.push_back(0.0);
            } else {
                size_t i = it - xp.begin();
                double x1 = xp[i - 1];
                double x2 = xp[i];
                double y1 = fp[i - 1];
                double y2 = fp[i];
                double y = y1 + ((y2 - y1) / (x2 - x1)) * (xi - x1);
                interpolatedValues.push_back(y);
            }
        }
        return interpolatedValues;
    }

    mzrtFeaturePtrVec GeneratePseudoSpectrum::get_iso_win_ms1_features(int isolation_window_base_index, double window_size) {
        mzrtFeaturePtrVec selected_features;
        double target_window = isolation_window_base_mz_[isolation_window_base_index];
        for (const auto &ms1_feature: ms1_features_) {
            std::vector<double> envelope_intensity = ms1_feature->getEnvelopeInte();
            std::vector<double> envelope_mz = ms1_feature->getEnvelopeMz();
            std::vector<double> xic = ms1_feature->getXic();
            double envelope_intensity_sum = std::accumulate(envelope_intensity.begin(), envelope_intensity.end(), 0.0);
            if (envelope_intensity_sum > 0) {
                for (double base_mz: isolation_window_base_mz_) {
                    if (base_mz == target_window and std::accumulate(xic.begin(), xic.end(), 0.0) > 0) {
                        double window_start_mz = base_mz - window_size / 2;
                        double window_end_mz = base_mz + window_size / 2;
                        // get envelope peaks intensity present within isolation window
                        double env_inte_sum_iso_win = 0;
                        for (size_t idx = 0; idx < envelope_intensity.size(); idx++)
                            if (window_start_mz <= envelope_mz[idx] && envelope_mz[idx] < window_end_mz)
                                env_inte_sum_iso_win += envelope_intensity[idx];
                        double coverage = env_inte_sum_iso_win / envelope_intensity_sum;
                        if (coverage > 0.5) {
                            selected_features.push_back(ms1_feature);
                        }
                    }
                }
            }
        }
        std::sort(selected_features.begin(), selected_features.end(), compareFeaturesInte);
        return selected_features;
    }

    int GeneratePseudoSpectrum::get_apex_cycle_distance(const mzrtFeaturePtr &ms1_feature, const mzrtFeaturePtr &ms2_feature) {
        int apex_cycle_distance_tole = std::min(3, ms1_feature->getCycleSpan() / 2);
        int apex_cycle_distance = std::abs(ms1_feature->getApexCycle() - ms2_feature->getApexCycle());
        if (apex_cycle_distance > apex_cycle_distance_tole) {
            std::vector<double> ms1_xic = ms1_feature->getXic();
            std::vector<double> ms2_xic = ms2_feature->getXic();
            if (ms2_feature->getCycleSpan() > 15)
                ms2_xic = moving_avg(ms2_feature->getXic(), 3);

            int ms1_apex_scan = static_cast<int>(std::distance(ms1_xic.begin(), std::max_element(ms1_xic.begin(), ms1_xic.end())));
            int ms2_apex_scan = static_cast<int>(std::distance(ms2_xic.begin(), std::max_element(ms2_xic.begin(), ms2_xic.end())));
            apex_cycle_distance = std::abs(ms1_apex_scan - ms2_apex_scan);
            return apex_cycle_distance;
        }
        return apex_cycle_distance;
    }

    std::vector<double> GeneratePseudoSpectrum::moving_avg(std::vector<double> xic, int size) {
        std::vector<double> smoothed_inte_list;
        std::vector<double> left_padding(1, 0);
        xic.insert(xic.begin(), left_padding.begin(), left_padding.end());
        int num_spec = static_cast<int>(xic.size());
        double sum = 0.0;
        int cnt = 0;
        for (int i = 0; i < num_spec; i++) {
            sum += xic[i];
            cnt++;
            if (cnt >= size) {
                smoothed_inte_list.push_back((sum / (double) size));
                sum -= xic[cnt - size];
            }
        }
        return smoothed_inte_list;
    }

    double GeneratePseudoSpectrum::computeSharedArea(const std::vector<double> &xic1, const std::vector<double> &xic2) {
        double sharedArea = 0.0;
        for (size_t i = 0; i < xic1.size(); ++i) {
            sharedArea += std::min(xic1[i], xic2[i]);
        }
        return sharedArea;
    }

    void GeneratePseudoSpectrum::score_pseudo_peaks(std::vector<PseudoPeaks> &pseudo_peak_list,
                                                    const mzrtFeaturePtr &ms1_feature) {
        std::sort(pseudo_peak_list.begin(), pseudo_peak_list.end(), comparePseudoPeaksInte);
        int total_peaks = static_cast<int>(pseudo_peak_list.size());
        double rank = total_peaks;
        for (int peak_idx = 0; peak_idx < total_peaks; peak_idx++) {
            pseudo_peak_list[peak_idx].setRank(rank / total_peaks);
            double score = get_pred(rank / total_peaks, pseudo_peak_list[peak_idx].getSharedInte(),
                                    (pseudo_peak_list[peak_idx].getMS2CycleSpan() * 1.0) / (ms1_feature->getCycleSpan() * 1.0));
            pseudo_peak_list[peak_idx].setScore(score);
            rank--;
        }
    }

    double GeneratePseudoSpectrum::get_pred(double intensity_ratio, double shared_area, double length_ratio) {
        double B0 = -3.349924626238689;
        double B1 = 1.8679961011204878;
        double B2 = 0.27006086659334383;
        double B3 = 3.983766800414337;
        double y = B0 + B1*intensity_ratio + B2*(length_ratio) + B3*shared_area;

        double py = 1 / (1 + std::exp(-y));
        return py;
    }

    std::vector<PseudoPeaks>
    GeneratePseudoSpectrum::filterPseudoPeaks(const EnvParaPtr &env_para_ptr, const mzrtFeaturePtr &ms1_feature,
                                              mzrtFeaturePtrVec &ms2_features_window,
                                              std::vector<PseudoPeaks> &pseudo_peak_list, double cutoff,
                                              int min_peak_num) {
        std::sort(pseudo_peak_list.begin(), pseudo_peak_list.end(), comparePseudoPeaksScore);
        std::vector<PseudoPeaks> low_mass_features;
        std::vector<PseudoPeaks> high_mass_features;

        int low_mass_num = env_para_ptr->compLowMassNum();
        int high_mass_num = env_para_ptr->compHighMassNum(ms1_feature->getMass());

        // add shorter features with correlation
        int counter = 0;
        for (const auto &i: pseudo_peak_list) {
            if (i.getScore() < cutoff and counter >= min_peak_num) continue;
            if (i.getMass() <= env_para_ptr->low_high_dividor_) {
                if ((int) low_mass_features.size() < low_mass_num) {
                    low_mass_features.push_back(i);
                    ms2_features_window[i.getMs2FeatureIdx()]->setUsedStatus(true);
                    counter++;
                }
            } else if ((int) high_mass_features.size() < high_mass_num) {
                high_mass_features.push_back(i);
                ms2_features_window[i.getMs2FeatureIdx()]->setUsedStatus(true);
                counter++;
            }
        }

        std::vector<PseudoPeaks> result;
        result.insert(std::end(result), std::begin(low_mass_features), std::end(low_mass_features));
        result.insert(std::end(result), std::begin(high_mass_features), std::end(high_mass_features));
        std::sort(result.begin(), result.end(), comparePseudoPeaksScore);
        return result;
    }

    void GeneratePseudoSpectrum::write_pseudo_spectrum(const TopdiaParaPtr &topdia_para_ptr, int ms1_feature_idx,
                                                       const mzrtFeaturePtr &ms1_feature,
                                                       std::vector<PseudoPeaks> &assigned_ms2_features) {
        std::string output_base_name = topdia_para_ptr->getOutputBaseName();
        std::string ms2_msalign_name = output_base_name + "_pseudo_ms2.msalign";
        int num_ms2_specs = 0;

        std::ofstream output;
        output.open(ms2_msalign_name, std::ios_base::app);
        output << std::fixed;

        output << "BEGIN IONS" << std::endl;
        output << "FRACTION_ID=" << ms1_feature_idx << std::endl;
        output << "FILE_NAME=" << ms2_msalign_name << std::endl;
        output << "SPECTRUM_ID=" << num_ms2_specs + ms1_feature_idx << std::endl; ///
        output << "TITLE=" << "Pseudo_Scan_" + std::to_string(ms1_feature_idx) << std::endl;
        output << "SCANS=" << num_ms2_specs + ms1_feature_idx << std::endl; ///
        output << "RETENTION_TIME=" << std::fixed << std::setprecision(2) << ms1_feature->getTimeApex() * 60.0
               << std::endl;
        output << "LEVEL=" << 2 << std::endl;
        output << "MS_ONE_ID=" << ms1_feature->getApexCycle() << std::endl; ///
        output << "MS_ONE_SCAN=" << num_ms2_specs + ms1_feature_idx << std::endl; ///
        output << "PRECURSOR_WINDOW_BEGIN="
               << ms1_feature->getBaseMz() - topdia_para_ptr->getPrecWindowWidth() / 2 << std::endl;
        output << "PRECURSOR_WINDOW_END="
               << ms1_feature->getBaseMz() + topdia_para_ptr->getPrecWindowWidth() / 2 << std::endl;
        output << "ACTIVATION=HCD" << std::endl;
        output << "PRECURSOR_MZ=" << ms1_feature->getMonoMz() << std::endl;
        output << "PRECURSOR_CHARGE=" << ms1_feature->getCharge() << std::endl;
        output << "PRECURSOR_MASS=" << ms1_feature->getMass() << std::endl;
        output << "PRECURSOR_INTENSITY=" << ms1_feature->getIntensity() << std::endl;
        output << "PRECURSOR_LENGTH=" << ms1_feature->getCycleSpan() << std::endl;
        output << "PRECURSOR_FEATURE_ID=" << ms1_feature_idx << std::endl; ///

        for (const auto &peak: assigned_ms2_features) {
            output << std::fixed << std::setprecision(5) << peak.getMass();
            output << "\t" << std::fixed << std::setprecision(2) << peak.getIntensity();
            output << "\t" << peak.getCharge();
            output << "\t" << std::fixed << std::setprecision(2) << peak.getScore();
            output << "\t" << std::fixed << std::setprecision(2) << peak.getApexDiffScan();
            output << "\t" << std::fixed << std::setprecision(2) << peak.getRank();
            output << "\t" << std::fixed << std::setprecision(2) << peak.getMS2CycleSpan();
            output << "\t" << std::fixed << std::setprecision(2) << peak.getSharedInte();
            output << "\t" << peak.getMs2ApexCycle();
            output << std::endl;
        }
        output << "END IONS" << std::endl;
        output << std::endl;
        output.close();
    }
}