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

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <src/envcnn/fdeep/tensor5.hpp>
#include "generate_matrix.hpp"
#include "ms/spec/baseline_util.hpp"

namespace toppic {

    std::vector<fdeep::tensor5s> GenerateMatrix::getTensor(MatchEnvPtrVec &ori_envs, PeakPtrVec &peak_list) {
        std::vector<fdeep::tensor5s> tensorsL;
        double baseline_intensity;
        getBaseLineUsingPeaklist(peak_list, baseline_intensity);
        std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpScoreDec);
        for (auto & ori_env : ori_envs) {
            // extract theoretical mass and intensities into separate vectors
            std::vector<double> theo_mass;
            std::vector<double> theoIntes;
            double max_inte;
            double theo_mono_mz;
            getTheoEnvData(ori_env, theo_mass, theoIntes, max_inte, theo_mono_mz);

            // extract intervaled experimental mass and intensities into separate vectors
            std::vector<double> peak_mass;
            std::vector<double> peak_intes;
            getExpEnvData(peak_list, ori_env, theo_mass, peak_mass, peak_intes);

            // Normalize Inte
            // Compute max theo inte
            double tolerance;

            std::vector<std::vector<double>> matrix  = initializeMatrix(tolerance);

            double min_mz = theo_mono_mz - 0.1;
            for (size_t k = 0; k < theo_mass.size(); k++) {
                double t_peak_mass;
                double t_peak_inte;
                double exp_inte;
                double inte_diff;
                double md;
                extractFeature(theo_mass, theoIntes, peak_mass, peak_intes, max_inte, tolerance, k, t_peak_mass,
                               t_peak_inte, exp_inte, inte_diff, md);
                populateMatrix(baseline_intensity, max_inte, matrix, min_mz, t_peak_mass, t_peak_inte,
                               exp_inte, inte_diff, md);
            }

            /// Note: Case where we have already added an experimental peak with the noise!
            addNoisePeaksInMatrix(peak_mass, peak_intes, max_inte, matrix, min_mz);
            generateTensors(tensorsL, matrix);
        }
        return tensorsL;
    }

    void GenerateMatrix::getExpEnvData(const PeakPtrVec &peak_list, MatchEnvPtr &ori_env, std::vector<double> &theo_mass,
                                  std::vector<double> &peak_mass, std::vector<double> &peak_intes) {
        RealEnvPtr real_env = ori_env->getRealEnvPtr();
        double real_mono_mz = real_env->getMonoMz();
        getExpIntervaledPeakData(peak_list, real_mono_mz, theo_mass, peak_mass, peak_intes);
    }

    void GenerateMatrix::getTheoEnvData(MatchEnvPtr &ori_env, std::vector<double> &theo_mass, std::vector<double> &theoIntes,
                                   double &max_inte, double &theo_mono_mz) {
        EnvelopePtr theo_env = ori_env->getTheoEnvPtr();
        theo_mono_mz= theo_env->getMonoMz();
        extractTheoPeakData(theo_env, theo_mass, theoIntes);
        normalizeTheoIntens(theoIntes, max_inte);
    }

    void GenerateMatrix::generateTensors(std::vector<fdeep::tensor5s> &tensorsL,
                                    const std::vector<std::vector<double>> &matrix) {
        fdeep::shape5 tensor5_shape(1, 1, 1, 300, 5);
        fdeep::tensor5 t(tensor5_shape, 0.0f);
        for (int y = 0; y < 300; ++y)
            for (int x = 0; x < 5; ++x)
                t.set(0, 0, 0, y, x, matrix[y][x]);

        std::vector<fdeep::tensor5> tensors;
        tensors.push_back(t);
        tensorsL.push_back(tensors);
    }

    void GenerateMatrix::addNoisePeaksInMatrix(const std::vector<double> &peak_mass,
                                                const std::vector<double> &peak_intes, double max_inte,
                                                std::vector<std::vector<double>> &matrix, double min_mz) {
        for (size_t exp_id = 0; exp_id < peak_mass.size(); exp_id++) {
            double exp_peak = peak_mass[exp_id];
            int bin_index = int((exp_peak - min_mz) * 100);
            // Evaluate Peak Condition
            bool peak_condition = false;
            for (int i = 0; i < 3; i++) {
                // to accomodate +2 and -2 bins - reason tolerance of 0.02
                if ((bin_index + i < 300) && (bin_index - i > -1) && (matrix[bin_index][1] == 0))
                    if (matrix[bin_index - i][0] == 0 && matrix[bin_index + i][0] == 0)
                        peak_condition = true;
                    else {
                        peak_condition = false;
                        break;
                    }
            }
            if (peak_condition == 1)
                matrix[bin_index][1] = (peak_intes[exp_id] / max_inte);
        }
    }

    void GenerateMatrix::extractFeature(const std::vector<double> &theo_mass, const std::vector<double> &theo_intes,
                                        const std::vector<double> &peak_mass, const std::vector<double> &peak_intes,
                                        double max_inte, double tolerance, size_t k, double &t_peak_mass,
                                        double &t_peak_inte, double &exp_inte, double &inte_diff, double &md) {
        t_peak_mass= theo_mass[k];
        t_peak_inte= theo_intes[k];
        exp_inte= 0;
        inte_diff= t_peak_inte;
        md= 0;
        int exp_peak_idx = -1;
        double mass_diff;
        double old_mass_diff = std::numeric_limits<double>::max();
        for (size_t exp_id = 0; exp_id < peak_mass.size(); exp_id++) {
            double exp_peak = peak_mass[exp_id];
            mass_diff = t_peak_mass - exp_peak;
            if (mass_diff <= tolerance && mass_diff >= -tolerance) {
                if (abs(mass_diff) < abs(old_mass_diff)){
                    exp_peak_idx = exp_id;
                    old_mass_diff = mass_diff;
                }
            }
        }
        mass_diff = -t_peak_mass;
        if (exp_peak_idx > -1) {
            double exp_peak = peak_mass[exp_peak_idx];
            exp_inte = peak_intes[exp_peak_idx] / max_inte;
            mass_diff = t_peak_mass - exp_peak;
            inte_diff = t_peak_inte - exp_inte;
        }

        mass_diff = abs(mass_diff);
        if (mass_diff == 0)
            md = 1;
        else if (mass_diff > 0.02)
            md = 0;
        else
            md = (0.02 - mass_diff)/0.02;
    }

    void GenerateMatrix::populateMatrix(double baseline_intensity, double max_inte,
                                        std::vector<std::vector<double>> &matrix, double min_mz,
                                        double t_peak_mass, double t_peak_inte, double exp_inte,
                                        double inte_diff, double md) {
        int bin_index = int((t_peak_mass - min_mz) * 100);
        if (bin_index < 300 && bin_index >= 0 ){
            matrix[bin_index][0] = t_peak_inte;
            matrix[bin_index][1] = exp_inte;
            matrix[bin_index][2] = md;
            matrix[bin_index][3] = inte_diff;
            matrix[bin_index][4] = log10(max_inte/baseline_intensity);
        }
    }

    void GenerateMatrix::normalizeTheoIntens(std::vector<double> &theo_intes, double &max_inte) {
        max_inte= *std::max_element(std::begin(theo_intes), std::end(theo_intes));
        for (double & theo_inte : theo_intes) {
            theo_inte = theo_inte / max_inte;
        }
    }

    void GenerateMatrix::getExpIntervaledPeakData(const PeakPtrVec &peak_list, double real_mono_mz,
                                                  std::vector<double> &theo_mass, std::vector<double> &peak_mass,
                                                  std::vector<double> &peak_intes) {
        double max_theo_mass = *std::max_element(std::begin(theo_mass), std::end(theo_mass));
        double min_theo_mass = *std::min_element(std::begin(theo_mass), std::end(theo_mass));
        for (const auto & k : peak_list) {
            double temp_mz = k->getPosition();
            if (temp_mz >= real_mono_mz - 5 && temp_mz <= real_mono_mz + 5) {
                if ((temp_mz >= min_theo_mass - 0.1) && (temp_mz <= max_theo_mass + 0.1)){
                    peak_mass.push_back(temp_mz);
                    peak_intes.push_back(k->getIntensity());
                }
            }
        }
    }

    void GenerateMatrix::extractTheoPeakData(EnvelopePtr &theo_env, std::vector<double> &theo_mass,
                                             std::vector<double> &theo_intes) {
        for (int i = 0; i < theo_env->getPeakNum(); i++) {
            theo_mass.push_back(theo_env->getMz(i));
            theo_intes.push_back(theo_env->getIntensity(i));
        }
    }

    std::vector<std::vector<double>> GenerateMatrix::initializeMatrix(double &tolerance) {
        tolerance= 0.02;
        int row = 300;
        int col = 5;
        std::vector<std::vector<double>> matrix(row);
        for (int i = 0; i < row; i++) {
            matrix[i] = std::vector<double>(col);
            for (int j = 0; j < col; j++)
                matrix[i][j] = 0;
        }
        return matrix;
    }

    void GenerateMatrix::getBaseLineUsingPeaklist(PeakPtrVec &peak_list, double &baseline_intensity) {
        std::vector<double> intes;
        for (auto & i : peak_list) {
            intes.push_back(i->getIntensity());
        }
        baseline_intensity= baseline_util::getBaseLine(intes);
    }
}

