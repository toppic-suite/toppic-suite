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
#include <iostream>
#include <math.h> 
#include <exception>
#include <iomanip> 


#include "fdeep/fdeep.hpp"
#include "spec/baseline_util.hpp"
#include "deconv/env/match_env.hpp"
#include "common/util/logger.hpp"
#include "deconv/env/match_env_filter.hpp" 

namespace toppic {

	MatchEnvPtrVec MatchEnvFilter::filter(MatchEnvPtrVec &ori_envs, double prec_mass,
		EnvParaPtr env_para_ptr, int ms_level_) {
		MatchEnvPtrVec low_mass_envs;
		MatchEnvPtrVec high_mass_envs;
		if (ms_level_ > 1) 
		 	std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpScoreDec);
		else
			std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpScoreDec_topfd);
		
		int low_mass_num = (int)(env_para_ptr->low_high_dividor_ / env_para_ptr->aa_avg_mass_ * env_para_ptr->peak_density_);
		int high_mass_num = (int)((prec_mass - env_para_ptr->low_high_dividor_)
			/ env_para_ptr->aa_avg_mass_ * env_para_ptr->peak_density_);
		for (size_t i = 0; i < ori_envs.size(); i++) {
			if (ori_envs[i]->getRealEnvPtr()->getMonoMass() <= env_para_ptr->low_high_dividor_) {
				if ((int)low_mass_envs.size() < low_mass_num) {
					low_mass_envs.push_back(ori_envs[i]);
				}
			}
			else {
				if ((int)high_mass_envs.size() < high_mass_num) {
					high_mass_envs.push_back(ori_envs[i]);
				}
			}
		}
		MatchEnvPtrVec result;
		result.insert(std::end(result), std::begin(low_mass_envs), std::end(low_mass_envs));
		result.insert(std::end(result), std::begin(high_mass_envs), std::end(high_mass_envs));
		if (ms_level_ > 1) 
			std::sort(result.begin(), result.end(), MatchEnv::cmpScoreDec);
		else
			std::sort(result.begin(), result.end(), MatchEnv::cmpScoreDec_topfd);
		
		return result;
	}
	
	MatchEnvPtrVec MatchEnvFilter::filter_using_cnn(MatchEnvPtrVec &ori_envs, PeakPtrVec &peak_list) {
		// Compute baseline intensity
		std::vector<double> intes;
		for (size_t i = 0; i < peak_list.size(); i++) {
			intes.push_back(peak_list[i]->getIntensity());
		}
		double baseline_intensity = baseline_util::getBaseLine(intes);
		std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpScoreDec_topfd);
				
		// Load model and initialize the tensor for evaluation
		const auto model = fdeep::load_model("/data/abbash/TopFD_Variants/Version_1_TopFd/proteomics_cpp_Deep_v2/toppic_resources/my_model.json");
		std::vector<fdeep::tensor5> tensors;

		// Compute prediction Score
		for (size_t i = 0; i < ori_envs.size(); i++) {
			RealEnvPtr real_env = ori_envs[i]->getRealEnvPtr();
			EnvelopePtr theo_env = ori_envs[i]->getTheoEnvPtr();
			int charge = theo_env->getCharge();
			double theo_mono_mz = theo_env->getMonoMz();
			double theo_mono_mass = theo_env->getMonoMass();
			double real_mono_mz = real_env->getMonoMz();
			double real_mono_mass = real_env->getMonoMass();

			// extract theoretical mass and intensities into separate vectors
			std::vector<double> theo_mass;
			std::vector<double> theo_intes;
			for (int i = 0; i < theo_env->getPeakNum(); i++) {
				theo_mass.push_back(theo_env->getMz(i));
				theo_intes.push_back(theo_env->getIntensity(i));
			}

			// Compute max theo inte
			double max_inte = *std::max_element(std::begin(theo_intes), std::end(theo_intes));

			// Compute max and min theo mass for defining intervaled peaks
			double max_theo_mass = *std::max_element(std::begin(theo_mass), std::end(theo_mass));
			double min_theo_mass = *std::min_element(std::begin(theo_mass), std::end(theo_mass));
			
			// extract intervaled experimental mass and intensities into separate vectors
			std::vector<double> peak_mass;
			std::vector<double> peak_intes;
			for (size_t k = 0; k < peak_list.size(); k++) {
				double temp_mz = peak_list[k]->getPosition();
				if (temp_mz >= real_mono_mz - 5 && temp_mz <= real_mono_mz + 5) {
					if ((temp_mz >= min_theo_mass - 0.1) && (temp_mz <= max_theo_mass + 0.1)){
						peak_mass.push_back(temp_mz);
						peak_intes.push_back(peak_list[k]->getIntensity());
					}
				}
			}

			// Normalize Inte
			for (size_t k = 0; k < theo_intes.size(); k++) {
				theo_intes[k] = theo_intes[k] / max_inte;
			}

			// Initialize Matrix
			int row = 300;
			int col = 7;
			double tolerance = 0.02;
			std::vector<std::vector<double>> matrix(row);
			for (int i = 0; i < row; i++) {
				matrix[i] = std::vector<double>(col);
				for (int j = 0; j < col; j++)
					matrix[i][j] = 0;
			}

			// min theoretical mass for computing the bin index in the matrix
			double min_mz = theo_mono_mz - 0.1;
			//min_mz = floor(min_mz * 100)/100;
			
			// Extract feature, compute bin index and populate the matrix
			for (size_t k = 0; k < theo_mass.size(); k++) {
				double t_peak_mass = theo_mass[k];
				double t_peak_inte = theo_intes[k];
				
				// Find experimental peak index 
				int exp_peak_idx = -1;
				for (size_t exp_id = 0; exp_id < peak_mass.size(); exp_id++) {
					double exp_peak = peak_mass[exp_id];
					double mass_diff = t_peak_mass - exp_peak;
					if (mass_diff <= tolerance && mass_diff >= -tolerance) {
						exp_peak_idx = exp_id;
						break;
					}
				}
				
				// Extract Feature
				double exp_inte = 0;
				double mass_diff = -t_peak_mass;
				double inte_diff = t_peak_inte;
				if (exp_peak_idx > -1) {
					double exp_peak = peak_mass[exp_peak_idx];
					exp_inte = peak_intes[exp_peak_idx] / max_inte;
					mass_diff = t_peak_mass - exp_peak;
					inte_diff = t_peak_inte - exp_inte;
				}

				// Populate Matrix
				int bin_index = int((t_peak_mass - min_mz) * 100);
				if (bin_index < 300 && bin_index >= 0 ){
					//std::cout << "exp_peak_idx: " << exp_peak_idx << ", ";
					matrix[bin_index][0] = t_peak_inte;
					matrix[bin_index][1] = exp_inte;
					matrix[bin_index][2] = charge;
					matrix[bin_index][3] = mass_diff;
					matrix[bin_index][4] = inte_diff;
					matrix[bin_index][5] = log10(max_inte);
					matrix[bin_index][6] = log10(baseline_intensity);
				}
			}
			
			// Add noise experimental peaks into the matrix
			/// Note: Case where we have already added an experimental peak with the noise!
			for (size_t exp_id = 0; exp_id < peak_mass.size(); exp_id++) {
				double exp_peak = peak_mass[exp_id];
      			int bin_index = int((exp_peak - min_mz) * 100);
				// Evaluate Peak Condition
      			bool peak_condition = 0;
	  			for (int i = 0; i < 3; i++) {
					// to accomodate +2 and -2 bins - reason tolerance of 0.02
       				if ((bin_index + i < 300) && (bin_index - i > -1) && (matrix[bin_index][1] == 0))
        				if ((matrix[bin_index - i][1] == 0) && (matrix[bin_index + i][1] == 0))
          					peak_condition = 1;
        				else{
          					peak_condition = 0;
							break;
						}
				}
	  			if (peak_condition == 1)
        			matrix[bin_index][1] = (peak_intes[exp_id] / max_inte);
			}

			// Initilize a tensor5 element and copy the values into tensor5
			fdeep::shape5 tensor5_shape(1, 1, 1, 300, 7);
			fdeep::tensor5 t(tensor5_shape, 0.0f);
			for (int y = 0; y < 300; ++y)
				for (int x = 0; x < 7; ++x)
					t.set(0, 0, 0, y, x, matrix[y][x]);
			tensors.push_back(t);

			// Compute prediction score and assign its value to the match_env
			const double pred_score = model.predict_single_output({ t });
			ori_envs[i]->setPredictionScore(pred_score);
			
			// Print Stats
			
			/*std::cout << std::endl << std::endl << "***********************************************" << std::endl;
			std::cout << std::setprecision(17);
			
			std::cout << "Experimental Peak: ";
			for (size_t exp_id = 0; exp_id < peak_mass.size(); exp_id++) 
				std::cout <<  peak_mass[exp_id] << ", ";
			std::cout << std::endl;	

			std::cout << "Theoretical Peaks Peak: ";
			for (size_t id = 0; id < theo_mass.size(); id++) 
				std::cout <<  theo_mass[id] << ", ";
			std::cout << std::endl;	

			std::cout << "Monoisotopic MZ: " << real_mono_mz << std::endl;
			std::cout << "Monoisotopic Mass: " << theo_mono_mass << std::endl;
			std::cout << "Pred Score: " << pred_score << std::endl;
			std::cout << "**** Matrix ****" << std::endl;
			for (int i = 0; i < matrix.size(); i++) {
			   for (int j = 0; j < matrix[i].size(); j++)
				   std::cout << matrix[i][j] << " ";
			   std::cout << std::endl;
		   }*/
		}
		return ori_envs;
	}
}
