#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "spec/msalign_reader.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"

#include <log4cxx/logger.h>

namespace prot {

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("ZeroPtmMng"));

ZeroPtmMng::ZeroPtmMng(std::string conf_file_name) {
  base_data_ptr_ = BaseDataPtr(new BaseData(conf_file_name));
}

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

  int spectra_num = countSpNum (mng_ptr->spectrum_file_name_.c_str(), 
                                mng_ptr->base_data_ptr_->getActivationPtrVec());

  LOG4CXX_DEBUG(logger, "spectra_number  " << spectra_num);


  /*
  searcher = new ZeroPtmSearcher(seqs, mng);

  String outputFileName = BioIo.getBaseName(mng.spectrumFileName) + "." + mng.outputFileExt;
  File spFile = new File(mng.spectrumFileName);
  int nSpectra = MsAlignReader.countSpNum(spFile);
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

