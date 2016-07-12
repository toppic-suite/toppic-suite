#include "base/logger.hpp"
#include "feature/real_env.hpp"
#include "feature/prec_env.hpp"
#include "feature/feature_ms_reader.hpp" 

namespace prot {

FeatureMsReader::FeatureMsReader(std::string &file_name) {
  reader_ptr_ = RawMsReaderPtr(new RawMsReader(file_name));
}

RawMsPtr FeatureMsReader::getNextMs(double prec_win_size) {
  reader_ptr_->readNext();
  PeakPtrVec peak_list = reader_ptr_->getPeakList();
  MsHeaderPtr header_ptr = reader_ptr_->getHeaderPtr();
  if (header_ptr == nullptr) {
    return nullptr;
  }
  RawMsPtr ms_ptr(new Ms<PeakPtr>(header_ptr, peak_list));

  // update ms_1 
  if (header_ptr->getMsLevel() == 1) {
    ms_one_ = ms_ptr;
  }
  if (do_refine_prec_mass_ && header_ptr->getMsLevel() == 2 && ms_one_ != nullptr) {
    refinePrecChrg(ms_one_, ms_ptr, prec_win_size);
  } else {
    if (header_ptr->getPrecSpMz() != 0.0) {
      header_ptr->setPrecMonoMz(header_ptr->getPrecSpMz());
    }
  }
  return ms_ptr;
}

// refine precursor charge and mz 
void FeatureMsReader::refinePrecChrg(RawMsPtr ms_one, RawMsPtr ms_two, 
                                     double prec_win_size) {
  MsHeaderPtr header_two = ms_two->getMsHeaderPtr();
  double prec_avg_mz = header_two->getPrecSpMz();
  int prec_charge = header_two->getPrecCharge();

  PeakPtrVec peak_list = ms_one->getPeakPtrVec();
  LOG_DEBUG("start refine precursor " << " peak num " << peak_list.size());
  RealEnvPtr env_ptr = PrecEnv::deconv(prec_win_size, peak_list, prec_avg_mz, 
                                       prec_charge);
  if (env_ptr != nullptr) {
    header_two->setPrecMonoMz(env_ptr->getMonoMz());
    header_two->setPrecCharge(env_ptr->getCharge());
    LOG_DEBUG("prec mz " << env_ptr->getMonoMz() << " prec charge " << env_ptr->getCharge());
  }
  else {
    LOG_DEBUG("EMPTY ENVELOPE POINTER");
  }
}



}
