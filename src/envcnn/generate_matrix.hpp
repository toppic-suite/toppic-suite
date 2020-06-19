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

#ifndef TOPPIC_TOPFD_ENVCNN_GENERATE_MATRIX_HPP_
#define TOPPIC_TOPFD_ENVCNN_GENERATE_MATRIX_HPP_

#include "ms/env/env_para.hpp"
#include "envcnn/fdeep/fdeep.hpp"
#include "ms/env/match_env.hpp"

namespace toppic {

    class GenerateMatrix {
    public:
        static std::vector<fdeep::tensor5s> getTensor(MatchEnvPtrVec &ori_envs, PeakPtrVec &peak_list);
        //MatchEnvPtrVec &ori_envs, double prec_mass, EnvParaPtr env_para_ptr);

        static void getBaseLineUsingPeaklist(PeakPtrVec &peak_list, double &baseline_intensity);

        static std::vector<std::vector<double>> initializeMatrix(double &tolerance) ;

        static void extractTheoPeakData(EnvelopePtr &theo_env, std::vector<double> &theo_mass,
                                        std::vector<double> &theo_intes);

        static void getExpIntervaledPeakData(const PeakPtrVec &peak_list, double real_mono_mz, std::vector<double> &theo_mass,
                                      std::vector<double> &peak_mass, std::vector<double> &peak_intes) ;

        static void normalizeTheoIntens(std::vector<double> &theo_intes, double &max_inte) ;

        static void populateMatrix(double baseline_intensity, double max_inte, std::vector<std::vector<double>> &matrix,
                       double min_mz, double t_peak_mass, double t_peak_inte, double exp_inte, double inte_diff,
                       double md) ;

        static void extractFeature(const std::vector<double> &theo_mass, const std::vector<double> &theo_intes,
                            const std::vector<double> &peak_mass, const std::vector<double> &peak_intes,
                            double max_inte,
                            double tolerance, size_t k, double &t_peak_mass, double &t_peak_inte, double &exp_inte,
                            double &inte_diff, double &md) ;

        static void addNoisePeaksInMatrix(const std::vector<double> &peak_mass, const std::vector<double> &peak_intes,
                              double max_inte, std::vector<std::vector<double>> &matrix, double min_mz) ;

        static void
        generateTensors(std::vector<fdeep::tensor5s> &tensorsL, const std::vector<std::vector<double>> &matrix) ;

        static void
        getTheoEnvData(MatchEnvPtr &ori_env, std::vector<double> &theo_mass, std::vector<double> &theoIntes,
                       double &max_inte,
                       double &theo_mono_mz);

        static void getExpEnvData(const PeakPtrVec &peak_list, MatchEnvPtr &ori_env, std::vector<double> &theo_mass,
                                  std::vector<double> &peak_mass, std::vector<double> &peak_intes);
    };

}

#endif
