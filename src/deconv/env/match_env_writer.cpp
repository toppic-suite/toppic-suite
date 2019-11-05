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

#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

#include "deconv/env/match_env.hpp"
#include "deconv/env/match_env_writer.hpp"
#include "deconv/env/match_env_writer.hpp"
#include "spec/baseline_util.hpp"

namespace toppic {

namespace match_env_writer {
	

void write_env(std::ofstream &file, MsHeaderPtr header, MatchEnvPtr match_env, PeakPtrVec &peak_list, int peakId, double min_inte) {
  EnvelopePtr theo_env = match_env->getTheoEnvPtr();
  RealEnvPtr real_env = match_env->getRealEnvPtr();
  double peakMZ;
  double Mass = real_env->getMonoMz();
  
  file << std::endl;
  file << "BEGIN ENVELOPE" << std::endl;
  file << "SPEC_ID=" << header->getId() << std::endl;
  file << "PEAK_ID=" << peakId << std::endl;
  file << "SPEC_SCAN=" << header->getScansString() << std::endl;
  file << "MS_LEVEL=" << header->getMsLevel() << std::endl;
  file << "REF_IDX=" << theo_env->getReferIdx() << std::endl;
  file << "CHARGE=" << theo_env->getCharge() << std::endl;
  file << "SCORE=" << match_env->getScore() << std::endl;
  file << "THEO_PEAK_NUM=" << theo_env->getPeakNum() << std::endl;
  file << "REAL_PEAK_NUM=" << (real_env->getPeakNum() - real_env->getMissPeakNum()) << std::endl;
  file << "THEO_MONO_MZ=" << theo_env->getMonoMz() << std::endl;
  file << "REAL_MONO_MZ=" << real_env->getMonoMz() << std::endl;
  file << "THEO_MONO_MASS=" << theo_env->getMonoMass() << std::endl;
  file << "REAL_MONO_MASS=" << real_env->getMonoMass() << std::endl;
  file << "THEO_INTE_SUM=" << theo_env->compIntensitySum() << std::endl;
  file << "REAL_INTE_SUM=" << real_env->compIntensitySum() << std::endl;
  file << "BASELINE_INTE=" << min_inte << std::endl;

  file << "Theoretical Peak MZ values and Intensities" << std::endl;

  for (int i = 0; i < theo_env->getPeakNum(); i++) {
    file << theo_env->getMz(i) << " " << theo_env->getIntensity(i) << std::endl;
  }
 
  file << "Experimental Peak MZ values and Intensities" << std::endl;
    for (size_t k = 0; k < peak_list.size(); k++) {
	  peakMZ = peak_list[k]->getPosition();
	  if (peakMZ >= Mass - 5 && peakMZ <= Mass + 5) {
		  file << peakMZ << " " << peak_list[k]->getIntensity() << std::endl;
	  }
	}
	
	file << "END ENVELOPE" << std::endl;
}

void write_env_vec(std::ofstream &file, MsHeaderPtr header, const MatchEnvPtrVec & envs, PeakPtrVec &peak_list) {
	file << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	std::vector<double> intes;
	for (size_t i = 0; i < peak_list.size(); i++) {
	  intes.push_back(peak_list[i]->getIntensity());
	}
	double min_inte = baseline_util::getBaseLine(intes);
	for (size_t i = 0; i < envs.size(); i++) {
	write_env(file, header, envs[i], peak_list, i, min_inte);
	}
	file << std::endl;
}

void write(const std::string & file, MsHeaderPtr header, const MatchEnvPtrVec & envs, PeakPtrVec &peak_list) {
  std::ofstream of(file, std::ofstream::out|std::ofstream::app);
  write_env_vec(of, header, envs, peak_list);
  of.close();
}

}  // namespace msalign_writer

}  // namespace toppic
