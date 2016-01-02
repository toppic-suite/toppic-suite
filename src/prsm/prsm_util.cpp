#include "base/logger.hpp"
#include "spec/extend_ms_factory.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm_util.hpp"

namespace prot {

std::string PrsmUtil::getValueStr(std::string line) {
  int start = line.find(">");
  int end = line.find("<", start);
  std::string num_str = line.substr(start + 1, end - start - 1);
  //LOG_DEBUG(line << "  Num str: " << num_str);
  return num_str;
}

std::string PrsmUtil::getXmlLine(const std::vector<std::string> &str_vec,
                                 const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return str_vec[i];
    }
  }
  return "";
}

PrsmPtrVec PrsmUtil::selectSpeciesPrsms(const PrsmPtrVec &prsm_ptrs,int species_id){
  PrsmPtrVec select_prsm_ptrs;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    if(species_id == prsm_ptrs[i]->getProteoformPtr()->getSpeciesId()){
      select_prsm_ptrs.push_back(prsm_ptrs[i]);
    }
  }
  return select_prsm_ptrs;
}

std::vector<int> PrsmUtil::getSpeciesIds(const PrsmPtrVec &prsm_ptrs, std::string &seq_name){
  std::vector<int> species_ids;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    int new_id = prsm_ptrs[i]->getProteoformPtr()->getSpeciesId();
    if(prsm_ptrs[i]->getProteoformPtr()->getSeqName() == seq_name){
      bool flag= false;
      for(size_t j=0;j<species_ids.size();j++){
        if(species_ids[j]==new_id){
          flag=true;
          break;
        }
      }
      if(!flag){
        species_ids.push_back(new_id);
      }
    }
  }
  return species_ids;
}

std::vector<int> PrsmUtil::getSpeciesIds(const PrsmPtrVec &prsm_ptrs){
  std::vector<int> species_ids;
  for(size_t i=0;i<prsm_ptrs.size();i++){
    bool find = false;
    for(size_t j=0;j<species_ids.size();j++){
      if(species_ids[j]==prsm_ptrs[i]->getProteoformPtr()->getSpeciesId()){
        find = true;
        break;
      }
    }
    if(!find){
      species_ids.push_back(prsm_ptrs[i]->getProteoformPtr()->getSpeciesId());
    }
  }
  return species_ids;
}

bool isMatchMs(PrsmPtr prsm_ptr, MsHeaderPtr header_ptr) {
  int id = header_ptr->getId();
  std::string scan = header_ptr->getScansString();
  int prec_id = header_ptr->getPrecId();
  if (prsm_ptr->getSpectrumId() == id && prsm_ptr->getPrecursorId() == prec_id) {
    if (prsm_ptr->getSpectrumScan() != scan) {
      LOG_ERROR("Error in Prsm.");
    }
    return true;
  } else {
    LOG_TRACE("prsm spectrum id " << prec_id << " ms spectrum id " << id);
    return false;
  }
}


void PrsmUtil::addSpectrumPtrsToPrsms(PrsmPtrVec &prsm_ptrs, PrsmParaPtr prsm_para_ptr){
  MsAlignReader reader(prsm_para_ptr->getSpectrumFileName(), prsm_para_ptr->getGroupSpecNum());
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  SpectrumSetPtr spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr);
  //use prsm order information (ordered by spectrum id then prec id)
  int start_prsm = 0;
  while (spec_set_ptr != nullptr) {
    if (spec_set_ptr->isValid()) {
      DeconvMsPtrVec deconv_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
      int spectrum_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getId();
      int prec_id = deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getPrecId();
      LOG_TRACE("spectrum id " << spectrum_id);
      for(size_t i = start_prsm; i<prsm_ptrs.size(); i++){
        if(isMatchMs(prsm_ptrs[i], deconv_ms_ptr_vec[0]->getMsHeaderPtr())){
          prsm_ptrs[i]->setDeconvMsPtrVec(deconv_ms_ptr_vec);
          double new_prec_mass = prsm_ptrs[i]->getAdjustedPrecMass();
          prsm_ptrs[i]->setRefineMsVec(ExtendMsFactory::geneMsThreePtrVec(deconv_ms_ptr_vec, prsm_para_ptr->getSpParaPtr(),
                                                           new_prec_mass));
        }
        if ((spectrum_id == prsm_ptrs[i]->getSpectrumId() && prec_id < prsm_ptrs[i]->getPrecursorId()) ||
            spectrum_id < prsm_ptrs[i]->getSpectrumId()) {
          start_prsm = i;
          break;
        }
      }
    }
    spec_set_ptr = reader.getNextSpectrumSet(sp_para_ptr);
  }
  reader.close();
}


}

