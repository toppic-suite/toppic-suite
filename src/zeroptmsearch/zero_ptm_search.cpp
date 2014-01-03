#include <base/logger.hpp>
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

namespace prot {

void zeroPtmSearch(SpectrumSetPtr spec_set_ptr, int type,
                   ProteoformPtrVec &proteoform_ptr_vec, ZeroPtmMngPtr mng_ptr, 
                   PrSMPtrVec &prsms) {
  ExtendMsPtr ms_three = spec_set_ptr->getSpThree();

  ZpFastMatchPtrVec fast_matches = zeroPtmFastFilter(type, ms_three,
                                                     proteoform_ptr_vec, 
                                                     mng_ptr->zero_ptm_filter_result_num_);

  DeconvMsPtr deconv_ms = spec_set_ptr->getDeconvMs();
  ZpSlowMatchPtrVec slow_matches = zeroPtmSlowFilter(type, deconv_ms, fast_matches, mng_ptr ); 

  for (unsigned int i = 0; i < slow_matches.size(); i++) {
      prsms.push_back(slow_matches[i]->geneResult());
  }

  std::sort(prsms.begin(), prsms.end(), prsm_match_fragment_down);
  prsms.erase(prsms.begin() + 1, prsms.end());
}

void zeroPtmSearchProcess(ZeroPtmMngPtr mng_ptr) {
  BaseDataPtr base_data_ptr = mng_ptr->base_data_ptr_;
  
  ProteoformPtrVec raw_forms = readFastaToProteoform(mng_ptr->search_db_file_name_,
                                                     base_data_ptr->getAcidPtrVec(),
                                                     base_data_ptr->getFixModResiduePtrVec(),
                                                     base_data_ptr->getDefaultProtModPtr());
  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, 
                                  base_data_ptr->getResiduePtrVec(),
                                  base_data_ptr->getAllowProtModPtrVec());

  int spectra_num = countSpNum (mng_ptr->spectrum_file_name_.c_str(), 
                                base_data_ptr->getActivationPtrVec());
  LOG_DEBUG("spectra_number  " << spectra_num);

  MsAlignReader reader(mng_ptr->spectrum_file_name_.c_str(), 
                       base_data_ptr->getActivationPtrVec());
  //std::string output_file_name = basename(mng_ptr->spectrum_file_name_ 
  //                                        + "." + mng_ptr->output_file_ext_);
  //PrSMWriter comp_writer;
  //PrSMWriter pref_writer;
  //PrSMWriter suff_writer;
  //PrSMWriter internal_writer;
  //PrSMWriter all_writer;

  ProtModPtrVec prot_mod_ptr_list = base_data_ptr->getProtModPtrVec();
  double shift = base_data_ptr->getAcetylationProtModPtr()->getProtShift();
  IonTypePtrVec ion_type_ptr_list = base_data_ptr->getIonTypePtrVec();
  LOG_DEBUG("start reading");
  int n = 0;
  DeconvMsPtr ms_ptr = reader.getNextMs();
  LOG_DEBUG("init ms_ptr");

  while (ms_ptr.get() != nullptr) {
    n++;
    SpectrumSetPtr spec_set_ptr = getSpectrumSet(ms_ptr, 0, mng_ptr->sp_para_ptr_, 
                                                 shift, ion_type_ptr_list);
    if (spec_set_ptr.get() != nullptr) {
      PrSMPtrVec comp_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_COMPLETE, prot_mod_forms, 
                    mng_ptr, comp_prsms);
      //comp_writer.write(comp_prsms);
      //all_write.write(comp_prsms);
      PrSMPtrVec pref_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_PREFIX, prot_mod_forms, 
                    mng_ptr, pref_prsms);
      //pref_writer.write(pref_prsms);
      //all_write.write(pref_prsms);
      PrSMPtrVec suff_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_SUFFIX, raw_forms, 
                    mng_ptr, suff_prsms);
      //suff_writer.write(suff_prsms);
      //all_write.write(suff_prsms);
      PrSMPtrVec internal_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_SUFFIX, raw_forms, 
                    mng_ptr, internal_prsms);
      //internal_writer.write(internal_prsms);
      //all_write.write(internal_prsms);
      LOG_DEBUG("zero ptm search complete " << n);
    }
    ms_ptr = reader.getNextMs();
    LOG_DEBUG("spectrum " << n);
  }

  reader.close();
  //comp_writer.close();
  //prec_writer.close();
  //suff_writer.close();
  //internal_writer.close();
  std::cout << "Non-ptm search finished." << std::endl;
}

} // end namespace
