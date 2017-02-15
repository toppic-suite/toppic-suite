// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "feature/text_writer.hpp"

namespace prot {

void TextWriter::writeText(std::ofstream &file, MatchEnvPtrVec &envs, MsHeaderPtr header_ptr) {

  file << "Result of scan(s): " << header_ptr->getScansString() << std::endl;
  file << "Ms level: " << header_ptr->getMsLevel() << std::endl;
  file << "Precursor information: monoisotopic m/z: " << header_ptr->getPrecMonoMz()
      << " charge: " << header_ptr->getPrecCharge()
      << " monoisotopic mass: " << header_ptr->getPrecMonoMass() << std::endl;
  file << "Scans\tMs_Level\tEnvelope\tReference_Idx\tCharge\tScore\tTheo_Peak_Num"
      << "\tReal_Peak_Num\tTheo_Mono_Mz\tReal_Mono_Mz\tTheo_Mono_Mass\tReal_Mono_Mass"
      << "\tTheo_Intensity_Sum\tReal_Intensity_Sum" << std::endl;

  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    EnvelopePtr theo_env = env->getTheoEnvPtr();
    RealEnvPtr real_env = env->getRealEnvPtr();
    file << header_ptr->getScansString() << "\t"
        << header_ptr->getMsLevel() << "\t"
        << (i + 1) << "\t"
        << theo_env->getReferIdx() << "\t"
        << theo_env->getCharge() << "\t"
        << env->getScore() << "\t"
        << theo_env->getPeakNum() << "\t"
        << (real_env->getPeakNum() - real_env->getMissPeakNum()) << "\t"
        << theo_env->getMonoMz() << "\t"
        << real_env->getMonoMz() << "\t"
        << theo_env->getMonoMass() << "\t"
        << real_env->getMonoMass() << "\t"
        << theo_env->compIntensitySum() << "\t"
        << real_env->compIntensitySum() << std::endl;
  }
  file << std::endl; 
}

}
