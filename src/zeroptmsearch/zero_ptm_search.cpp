#include <base/logger.hpp>
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"
#include "zeroptmsearch/zero_ptm_search.hpp"

namespace prot {

void zeroPtmSearchProcess(ZeroPtmMngPtr mng_ptr) {
  
  ProteoformPtrVec raw_forms = readFastaToProteoform(mng_ptr->search_db_file_name_,
                                                     mng_ptr->base_data_ptr_->getAcidPtrVec(),
                                                     mng_ptr->base_data_ptr_->getFixModResiduePtrVec());
  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, 
                                  mng_ptr->base_data_ptr_->getResiduePtrVec(),
                                  mng_ptr->base_data_ptr_->getAllowProtModPtrVec());

  ProteoformPtrVec prot_mod_none_forms = getProtModNoneProteoform(prot_mod_forms);

  int spectra_num = countSpNum (mng_ptr->spectrum_file_name_.c_str(), 
                                mng_ptr->base_data_ptr_->getActivationPtrVec());
  LOG_DEBUG("spectra_number  " << spectra_num);

  MsAlignReader reader(mng_ptr->spectrum_file_name_.c_str(), 
                       mng_ptr->base_data_ptr_->getActivationPtrVec());

  ProtModPtrVec prot_mod_ptr_list = mng_ptr->base_data_ptr_->getProtModPtrVec();
  double shift = getProtModAcetylationShift(prot_mod_ptr_list);
  IonTypePtrVec ion_type_ptr_list = mng_ptr->base_data_ptr_->getIonTypePtrVec();
  LOG_DEBUG("start reading");
  int n = 0;
  DeconvMsPtr ms_ptr = reader.getNextMs();
  LOG_DEBUG("init ms_ptr");

  while (ms_ptr.get() != nullptr) {
    n++;
    SpectrumSetPtr spec_set_ptr = getSpectrumSet(ms_ptr, 0, mng_ptr->sp_para_ptr_, shift, 
                                                 ion_type_ptr_list);

    if (spec_set_ptr.get() != nullptr) {
      PrSMPtrVec comp_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_COMPLETE, prot_mod_forms, 
                    mng_ptr, comp_prsms);
      PrSMPtrVec pref_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_PREFIX, prot_mod_forms, 
                    mng_ptr, pref_prsms);
      PrSMPtrVec suff_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_SUFFIX, prot_mod_none_forms, 
                    mng_ptr, suff_prsms);
      PrSMPtrVec internal_prsms;
      zeroPtmSearch(spec_set_ptr, SEMI_ALIGN_TYPE_SUFFIX, prot_mod_none_forms, 
                    mng_ptr, internal_prsms);
      LOG_DEBUG("zero ptm search complete " << n);
    }
    ms_ptr = reader.getNextMs();
    LOG_DEBUG("spectrum " << n);
  }

  /*
  String outputFileName = BioIo.getBaseName(mng.spectrumFileName) + "." + mng.outputFileExt;

  MsAlignReader spReader = new MsAlignReader(spFile);
  PrSMXmlWriter writers[] = new PrSMXmlWriter[4];
  for (int i = 0; i < 4; i++) {
    writers[i] = new PrSMXmlWriter(new File(outputFileName +"_" + SemiAlignType.getAlignmentType(i).getName()));
  }
  PrSMXmlWriter allWriter = new PrSMXmlWriter(new File(outputFileName));

  System.out.println("zero-ptm search started.");
  int cnt = 0;
  long startTime = System.currentTimeMillis();
  Ms<DeconvPeak>[] deconvSp;
  PrSM prsms[][] = new PrSM[4][mng.nReport];
  while ((deconvSp = spReader.getNextMses()) != null) {
    cnt++;
    for (int i = 0; i < deconvSp.length; i++) {
      SpectrumSet spectrumSet = SpectrumSet.getSpectrumSet(
          deconvSp[i], 0, mng.spPara);
      if (spectrumSet != null) {
        String scan = deconvSp[i].getHeader().getScansString();
        String msg = ProcessUtil.updateMsg("Zero-ptm search", scan,
                                           nSpectra, cnt, startTime);
        System.out.print(msg);
        searcher.search(spectrumSet, prsms);
        allWriter.write(prsms);
        for (int j = 0; j < 4; j++) {
          writers[j].write(prsms[j]);
        }
      }
    }
  }
  spReader.close();
  allWriter.close();
  for (int i = 0; i < 4; i++) {
    writers[i].close();
  }
  System.out.println("\nNon-ptm search finished.");
  */
}

void zeroPtmSearch(SpectrumSetPtr spec_set_ptr, int type,
                   ProteoformPtrVec &form_ptr_vec, ZeroPtmMngPtr mng_ptr, 
                   PrSMPtrVec &prsms) {
  ExtendMsPtr ms_three = spec_set_ptr->getSpThree();

  ZpFastMatchPtrVec fast_matches = zeroPtmFastFilter(type, ms_three,
                                                     form_ptr_vec, mng_ptr->zero_ptm_filter_result_num_);

  DeconvMsPtr deconv_ms = spec_set_ptr->getDeconvMs();
  ZpSlowMatchPtrVec slow_matches = zeroPtmSlowFilter(type, deconv_ms, fast_matches, mng_ptr ); 

  for (unsigned int i = 0; i < slow_matches.size(); i++) {
      prsms.push_back(slow_matches[i]->geneResult());
  }

  std::sort(prsms.begin(), prsms.end(), prsm_match_fragment_down);
  prsms.erase(prsms.begin() + 1, prsms.end());
}

} // end namespace
