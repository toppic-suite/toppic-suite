//
// Created by abbash on 8/30/22.
//

#ifndef TOPPIC_FEATURE_HPP
#define TOPPIC_FEATURE_HPP

#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/score/get_env_cnn_score.hpp"
#include "topfd/feature_detect/score/get_component_score.hpp"

namespace toppic {
class Feature {
public:
    Feature(EnvCollection &env_coll, PeakMatrix &peak_matrix, fdeep::model &model, int feature_id, double snr);

    int getFeatureId() const { return feature_id_; }
    void setFeatureId(int featureId) { feature_id_ = featureId; }

    int getMinScan() const { return min_scan_; }
    void setMinScan(int minScan) { min_scan_ = minScan; }

    int getMaxScan() const { return max_scan_; }
    void setMaxScan(int maxScan) { max_scan_ = maxScan; }

    int getMinCharge() const { return min_charge_; }
    void setMinCharge(int minCharge) { min_charge_ = minCharge; }

    int getMaxCharge() const { return max_charge_; }
    void setMaxCharge(int maxCharge) { max_charge_ = maxCharge; }

    double getMonoMass() const { return mono_mass_; }
    void setMonoMass(double monoMass) { mono_mass_ = monoMass; }

    int getRepCharge() const { return rep_charge_; }
    void setRepCharge(int repCharge) { rep_charge_ = repCharge; }

    double getRepMz() const { return rep_mz_; }
    void setRepMz(double repMz) { rep_mz_ = repMz; }

    double getAbundance() const { return abundance_; }
    void setAbundance(double abundance) { abundance_ = abundance; }

    double getMinElutionTime() const { return min_elution_time_; }
    void setMinElutionTime(double minElutionTime) { min_elution_time_ = minElutionTime; }

    double getMaxElutionTime() const { return max_elution_time_; }
    void setMaxElutionTime(double maxElutionTime) { max_elution_time_ = maxElutionTime; }

    double getApexElutionTime() const { return apex_elution_time_; }
    void setApexElutionTime(double apexElutionTime) { apex_elution_time_ = apexElutionTime; }

    double getElutionLength() const { return elution_length_; }
    void setElutionLength(double elutionLength) { elution_length_ = elutionLength; }

    double getEnvcnnScore() const { return envcnn_score_; }
    void setEnvcnnScore(double envcnnScore) { envcnn_score_ = envcnnScore; }

    double getPercentMatchedPeaks() const { return percent_matched_peaks_; }
    void setPercentMatchedPeaks(double percentMatchedPeaks) { percent_matched_peaks_ = percentMatchedPeaks; }

    double getIntensityCorrelation() const { return intensity_correlation_; }
    void setIntensityCorrelation(double intensityCorrelation) { intensity_correlation_ = intensityCorrelation; }

    double getTop3Correlation() const { return top3_correlation_; }
    void setTop3Correlation(double top3Correlation) { top3_correlation_ = top3Correlation; }

    double getEvenOddPeakRatios() const { return even_odd_peak_ratios_; }
    void setEvenOddPeakRatios(double evenOddPeakRatios) { even_odd_peak_ratios_ = evenOddPeakRatios; }

    double getPercentConsecPeaks() const { return percent_consec_peaks_; }
    void setPercentConsecPeaks(double percentConsecPeaks) { percent_consec_peaks_ = percentConsecPeaks; }

    int getNumTheoPeaks() const { return num_theo_peaks_; }
    void setNumTheoPeaks(int numTheoPeaks) { num_theo_peaks_ = numTheoPeaks; }

    double getMzErrorSum() const { return mz_error_sum_; }
    void setMzErrorSum(double mzErrorSum) { mz_error_sum_ = mzErrorSum; }

    double getScore() const { return score_; }
    void setScore(double score) { score_ = score;}

    int getLabel() const { return label_; }
    void setLabel(int label) { label_ = label; }

private:
    int feature_id_ = 0;
    int min_scan_ = 0;
    int max_scan_ = 0;
    int min_charge_ = 0;
    int max_charge_ = 0;
    double mono_mass_ = 0;
    int rep_charge_ = 0;
    double rep_mz_ = 0;
    double abundance_ = 0;
    double min_elution_time_ = 0;
    double max_elution_time_ = 0;
    double apex_elution_time_ = 0;
    double elution_length_ = 0;
    double envcnn_score_ = 0;
    double percent_matched_peaks_ = 0;
    double intensity_correlation_ = 0;
    double top3_correlation_ = 0;
    double even_odd_peak_ratios_ = 0;
    double percent_consec_peaks_ = 0;
    int num_theo_peaks_ = 0;
    double mz_error_sum_ = 0;
    double score_ = 0;
    int label_ = 0;
};
}


#endif //TOPPIC_FEATURE_HPP
