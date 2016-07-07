#include "feature/msalign_writer.hpp"

namespace prot {

void MsalignWriter::writeText(std::ofstream &file, MatchEnvPtrVec &envs, 
                              MsHeaderPtr header_ptr, int id) {
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << id << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (header_ptr->getMsLevel() > 1) {
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
  }

  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    EnvelopePtr theo_env = env->getTheoEnvPtr();
    RealEnvPtr real_env = env->getRealEnvPtr();
    file << real_env->getMonoMass();
    file << "\t" << theo_env->compIntensitySum();
    file << "\t" << theo_env->getCharge();
    file << std::endl;
  }
  file << "END IONS" << std::endl;
  file << std::endl; 
}

}
