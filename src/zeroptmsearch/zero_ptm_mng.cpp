#include <base/logger.hpp>
#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "spec/spectrum_set.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"

namespace prot {

ZeroPtmMng::ZeroPtmMng(std::string conf_file_name): 
    base_data_ptr_ (new BaseData(conf_file_name)),
    peak_tolerance_ptr_ (new PeakTolerance(ppo_, use_min_tolerance_, min_tolerance_)),
    extend_sp_para_ptr_ (new ExtendSpPara(extend_min_mass_, ext_offsets_)),
    sp_para_ptr_(new SpPara(min_peak_num_, min_mass_, peak_tolerance_ptr_, 
                           extend_sp_para_ptr_, base_data_ptr_->getActivationPtr())) 
  {}

void zeroPtmSearchProcess(ZeroPtmMngPtr mng_ptr) {
  /*
  ProteoformPtrVec ori_forms = readFastaToProteoform(mng_ptr->search_db_file_name_,
                                                    mng_ptr->base_data_ptr_->getAcidPtrVec(),
                                                    mng_ptr->base_data_ptr_->getFixModResiduePtrVec());
  for (unsigned int i = 0; i < 10; i++) {
    std::cout << ori_forms[i]->toString();
  }
  ProteoformPtrVec prot_mod_forms = generateProtModProteoform(ori_forms, 
                                                              mng_ptr->base_data_ptr_->getResiduePtrVec(),
                                                              mng_ptr->base_data_ptr_->getAllowProtModPtrVec());
  for (unsigned int i = 0; i < 10; i++) {
    std::cout << prot_mod_forms[i]->toString();
  }
  */
  LOG_DEBUG("start reading spectra.");

  int spectra_num = countSpNum (mng_ptr->spectrum_file_name_.c_str(), 
                                mng_ptr->base_data_ptr_->getActivationPtrVec());
  LOG_DEBUG("spectra_number  " << spectra_num);

  MsAlignReader reader(mng_ptr->spectrum_file_name_.c_str(), 
                       mng_ptr->base_data_ptr_->getActivationPtrVec());
  LOG_DEBUG("start reading");

  int n = 0;
  DeconvMsPtr ms_ptr = reader.getNextMs();
  ProtModPtrVec prot_mod_ptr_list = mng_ptr->base_data_ptr_->getProtModPtrVec();
  double shift = getProtModAcetylationShift(prot_mod_ptr_list);

  while (ms_ptr.get() != nullptr) {
    n++;
    //SpectrumSetPtr spec_set_ptr = getSpectrumSet(ms_ptr, 0, mng_ptr->sp_para_ptr_, shift);
    
    ms_ptr = reader.getNextMs();
  }

  /*
  searcher = new ZeroPtmSearcher(seqs, mng);

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

} /* namespace prot */

