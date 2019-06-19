//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#include "feature/peak_cluster.hpp"

namespace toppic {

PeakCluster::PeakCluster(MatchEnvPtr match_env) {
  theo_env_ = match_env->getTheoEnvPtr();
  rep_mass_ = theo_env_->getMonoMass();
  rep_charge_ = theo_env_->getCharge();

  int peak_num = theo_env_->getPeakNum();
  rep_summed_peaks_.resize(peak_num, 0.0);

  clearScores();

  flag_ = 0;
  init_score_ = false;
}

void PeakCluster::addEnvelopes(int min_charge, int max_charge, 
                               int min_ms1_id, int max_ms1_id, 
                               int scan_begin, int scan_end,
                               RealEnvPtrVec envs) {
  int row_num = max_charge - min_charge + 1;
  int col_num = max_ms1_id - min_ms1_id + 1;
  
  min_charge_ = min_charge;
  max_charge_ = max_charge;

  min_ms1_id_ = min_ms1_id;
  max_ms1_id_ = max_ms1_id;

  scan_begin_ = scan_begin;
  scan_end_ = scan_end;

  real_envs_.resize(row_num);

  for (int i = 0; i < row_num; i++) {
    real_envs_[i].resize(col_num);
  }

  for (size_t i = 0; i < envs.size(); i++) {
    int row = envs[i]->getCharge() - min_charge_;
    int col = envs[i]->getSpId() - min_ms1_id;
    if (row >= 0 && row < row_num && col >= 0 && col < col_num) {
      real_envs_[row][col] = envs[i];
    }
  }
}

void PeakCluster::clearScores() {
  inte_distr_.resize(2, 0.0);
  best_corr_scores_.resize(2, 0.0);
  best_inte_scores_.resize(2, 0.0);
  best_dist_scores_.resize(2, 1.0);

  best_charges_.resize(2, 0);
  env_dist_scores_.resize(2, 1.0);
  env_corr_scores_.resize(2, 0.0);
  env_inte_scores_.resize(2, 0.0);
  xic_corr_between_best_charges_.resize(2, 0.0);
}

void PeakCluster::updateScore(bool check_pvalue) {
  int row_num = max_charge_ - min_charge_ + 1;
  int col_num = max_ms1_id_ - min_ms1_id_ + 1;
  int base_peak_idx = theo_env_->getReferIdx(); 

  clearScores();

  std::vector<double> best_charge_dist{10.0, 10.0};

  // sum up peak intensities
  int peak_num = theo_env_->getPeakNum();
  std::vector<double> summed_inte(peak_num, 0);

  int xic_len = col_num + 18;
  int xic_start_idx = 9;

  std::vector<std::vector<double>> xic2(2);
  xic2[0].resize(xic_len, 0.0);
  xic2[1].resize(xic_len, 0.0);

  std::vector<std::vector<double>> charge_xic(row_num);

  double temp_best_bc_dist = 10.0;
  double rep_env_bc_dist = 10.0;
  RealEnvPtr rep_env(nullptr);

  double rep_env_bc_dist_2 = 10.0;
  RealEnvPtr rep_env_2(nullptr);

  std::vector<double> tmp_best_dist_score{10.0, 10.0};
  std::vector<double> tmp_best_inte_score(2, 0.0);
  std::vector<double> tmp_best_corr_score(2, 0.0);

  for (int i = 0; i < row_num; i++) {
    int charge = i + min_charge_;

  }

}

/*  
public void UpdateScore(List<Ms1Spectrum> ms1Spectra, bool pValueCheck = true)
{
  for (var i = 0; i < nRows; i++)
  {
    var charge = i + MinCharge;
    var mostAbuMz = TheoreticalEnvelope.GetIsotopeMz(charge, mostAbuIdx);
    Array.Clear(summedIntensity, 0, summedIntensity.Length);

    chargeXic[i] = new double[xicLen];

    var chargeIdx = (charge % 2 == 0) ? EvenCharge : OddCharge;
    var summedMostAbuIsotopeIntensity = 0d;
    var summedReferenceIntensity = 0d;

    for (var j = 0; j < nCols; j++)
    {
      var envelope = Envelopes[i][j];
      var col = minCol + j;

      var localWin = ms1Spectra[col].GetLocalMzWindow(mostAbuMz);

      if (envelope == null) continue;

      envelope.Peaks.SumEnvelopeTo(summedIntensity);
      var mostAbuPeak = envelope.Peaks[mostAbuIdx];

      if (mostAbuPeak != null && mostAbuPeak.Active)
      {
        summedMostAbuIsotopeIntensity += mostAbuPeak.Intensity;
        summedReferenceIntensity += localWin.HighestIntensity;
      }
      AbundanceDistributionAcrossCharge[chargeIdx] += envelope.Abundance;

      var newBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(envelope.Peaks);
      var newCorr = TheoreticalEnvelope.GetPearsonCorrelation(envelope.Peaks);

      var goodEnvelope = (newBcDist < 0.07 || newCorr > 0.7);

      if (goodEnvelope)
      {
        xic2[chargeIdx][xicStartIdx + j] += envelope.Abundance;
        chargeXic[i][xicStartIdx + j] = envelope.Abundance;
      }

      var levelOneEnvelope = true;
      var levelTwoEnvelope = true;

      if (pValueCheck)
      {
        var poissonPValue = localWin.GetPoissonTestPValue(envelope.Peaks, TheoreticalEnvelope.Size);
        var rankSumPValue = localWin.GetRankSumTestPValue(envelope.Peaks, TheoreticalEnvelope.Size);
        levelOneEnvelope = (rankSumPValue < 0.01 && poissonPValue < 0.01);
        //levelTwoEnvelope = (rankSumPValue < 0.05 || poissonPValue < 0.05);
      }

      if (levelOneEnvelope)
      {
        if (newBcDist < BestDistanceScoreAcrossCharge[chargeIdx])
        {
          BestDistanceScoreAcrossCharge[chargeIdx] = newBcDist;
          if (localWin.MedianIntensity > 0)
            BestIntensityScoreAcrossCharge[chargeIdx] = envelope.HighestIntensity / localWin.HighestIntensity;
          else BestIntensityScoreAcrossCharge[chargeIdx] = 1.0d;
        }

        BestCorrelationScoreAcrossCharge[chargeIdx] = Math.Max(BestCorrelationScoreAcrossCharge[chargeIdx], newCorr);

        if (newBcDist < repEnvelopeBcDist)
        {
          repEnvelopeBcDist = newBcDist;
          repEnvelope = envelope;
        }

        // in the initial scoring, classify major and minor envelopes
        if (!_initScore && goodEnvelope) envelope.GoodEnough = true;
      }

      if (levelTwoEnvelope)
      {
        if (newBcDist < tempBestDistanceScoreAcrossCharge[chargeIdx])
        {
          tempBestDistanceScoreAcrossCharge[chargeIdx] = newBcDist;
          if (localWin.MedianIntensity > 0)
            tempBestIntensityScoreAcrossCharge[chargeIdx] = envelope.HighestIntensity / localWin.HighestIntensity;
          else tempBestIntensityScoreAcrossCharge[chargeIdx] = 1.0d;
        }
        tempBestCorrelationScoreAcrossCharge[chargeIdx] = Math.Max(tempBestCorrelationScoreAcrossCharge[chargeIdx], newCorr);

        if (newBcDist < repEnvelopeBcDist2)
        {
          repEnvelopeBcDist2 = newBcDist;
          repEnvelope2 = envelope;
        }
      }
    }

    var bcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(summedIntensity);
    EnvelopeDistanceScoreAcrossCharge[chargeIdx] = Math.Min(bcDist, EnvelopeDistanceScoreAcrossCharge[chargeIdx]);
    EnvelopeCorrelationScoreAcrossCharge[chargeIdx] = Math.Max(TheoreticalEnvelope.GetPearsonCorrelation(summedIntensity), EnvelopeCorrelationScoreAcrossCharge[chargeIdx]);

    if (BestCharge[chargeIdx] < 1 || bcDist < bestChargeDist[chargeIdx])
    {
      BestCharge[chargeIdx] = charge;
      bestChargeDist[chargeIdx] = bcDist;
      if (summedReferenceIntensity > 0)
        EnvelopeIntensityScoreAcrossCharge[chargeIdx] = summedMostAbuIsotopeIntensity/summedReferenceIntensity;
      //if (summedMedianIntensity > 0) EnvelopeIntensityScoreAcrossCharge[chargeIdx] = Math.Min(1.0, 0.1*(summedMostAbuIsotopeIntensity / summedMedianIntensity));
    }

    if (bcDist < tempBestBcDist)
    {
      tempBestBcDist = bcDist;
      Array.Copy(summedIntensity, RepresentativeSummedEnvelop, RepresentativeSummedEnvelop.Length);
    }
  }

  // when good envelope is observed at only either even or odd charge...
  if (BestCorrelationScoreAcrossCharge[0] > 0.7 && BestCorrelationScoreAcrossCharge[1] < 0.5)
  {
    const int i = 1;
    BestCorrelationScoreAcrossCharge[i] = tempBestCorrelationScoreAcrossCharge[i];
    BestIntensityScoreAcrossCharge[i] = tempBestIntensityScoreAcrossCharge[i];
    BestDistanceScoreAcrossCharge[i] = tempBestDistanceScoreAcrossCharge[i];
  }

  if (BestCorrelationScoreAcrossCharge[1] > 0.7 && BestCorrelationScoreAcrossCharge[0] < 0.5)
  {
    const int i = 0;
    BestCorrelationScoreAcrossCharge[i] = tempBestCorrelationScoreAcrossCharge[i];
    BestIntensityScoreAcrossCharge[i] = tempBestIntensityScoreAcrossCharge[i];
    BestDistanceScoreAcrossCharge[i] = tempBestDistanceScoreAcrossCharge[i];
  }

  // normalize abundance across charges
  var s = AbundanceDistributionAcrossCharge[0] + AbundanceDistributionAcrossCharge[1];
  if (s > 0)
  {
    for (var chargeIdx = 0; chargeIdx < 2; chargeIdx++)
    {
      AbundanceDistributionAcrossCharge[chargeIdx] = AbundanceDistributionAcrossCharge[chargeIdx] / s;
    }
  }

  if (nCols > 1)
  {
    var evenChargeIdx = BestCharge[EvenCharge] - MinCharge;
    var oddChargeIdx = BestCharge[OddCharge] - MinCharge;
    XicCorrelationBetweenBestCharges[0] = FitScoreCalculator.GetPearsonCorrelation(Smoother.Smooth(chargeXic[evenChargeIdx]), Smoother.Smooth(chargeXic[oddChargeIdx]));
    XicCorrelationBetweenBestCharges[1] = FitScoreCalculator.GetPearsonCorrelation(Smoother.Smooth(xic2[EvenCharge]), Smoother.Smooth(xic2[OddCharge]));
  }

  if (repEnvelope == null && repEnvelope2 != null) repEnvelope = repEnvelope2;

  if (repEnvelope != null)
  {
    // set representative charge, mz and scanNum
    RepresentativeCharge = repEnvelope.Charge;
    RepresentativeMz = repEnvelope.RepresentativePeak.Mz;
    RepresentativeScanNum = repEnvelope.ScanNum;
  }

  _initScore = true;
}
*/

}

