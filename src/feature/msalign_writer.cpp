#include "feature/msalign_writer.hpp"

namespace prot {

void MsalignWriter::writeText(std::ofstream &file, MatchEnvPtrVec &envs, 
                              MsHeaderPtr header_ptr, int id) {
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << id << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  file << "RETENTION_TIME=" << header_ptr->getRetentionTime() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (header_ptr->getMsLevel() > 1) {
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
    file << "PRECURSOR_INTENSITY=" << header_ptr->getPrecInte() << std::endl;
    if (header_ptr->getFeatureId() >= 0) {
      file << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
    }
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

void MsalignWriter::writeText(std::ofstream &file, DeconvMsPtr ms_ptr) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  file << "BEGIN IONS" << std::endl;
  file << "ID=" << header_ptr->getId() << std::endl;
  file << "SCANS=" << header_ptr->getScansString() << std::endl;
  file << "RETENTION_TIME=" << header_ptr->getRetentionTime() << std::endl;
  if (header_ptr->getActivationPtr() != nullptr) {
    file << "ACTIVATION=" << header_ptr->getActivationPtr()->getName() << std::endl;
  }
  if (header_ptr->getMsLevel() > 1) {
    file << "PRECURSOR_MZ=" << header_ptr->getPrecMonoMz() << std::endl;
    file << "PRECURSOR_CHARGE=" << header_ptr->getPrecCharge() << std::endl;
    file << "PRECURSOR_MASS=" <<  header_ptr->getPrecMonoMass() << std::endl;
    file << "PRECURSOR_INTENSITY=" << header_ptr->getPrecInte() << std::endl;
    if (header_ptr->getFeatureId() >= 0) {
      file << "FEATURE_ID=" << header_ptr->getFeatureId() << std::endl;
    }
  }

  for (size_t i = 0; i < ms_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms_ptr->getPeakPtr(i);
    file << peak_ptr->getPosition();
    file << "\t" << peak_ptr->getIntensity();
    file << "\t" << peak_ptr->getCharge();
    file << std::endl;
  }
  file << "END IONS" << std::endl;
  file << std::endl; 
}

}
