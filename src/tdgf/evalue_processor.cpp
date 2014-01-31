#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "tdgf/evalue_processor.hpp"

namespace prot {

EValueProcessor::EValueProcessor(TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
}

void EValueProcessor::init() {
  BaseDataPtr base_data_ptr = mng_ptr_->base_data_ptr_;
  ProteoformPtrVec raw_forms 
      = readFastaToProteoform(mng_ptr_->search_db_file_name_,
                              base_data_ptr->getFixModResiduePtrVec());
  ProteoformPtrVec prot_mod_forms 
      = generateProtModProteoform(raw_forms, 
                                  base_data_ptr->getAllowProtModPtrVec());
  
  /*
  ResFreq resFreqs[] = ResFreqArrayUtil.getResFreqArray(resList
                                                        .getResidues(20));
  ResFreq nResFreqs[] = resFreqs;
  compP = new CompPValueArray(seqs, nModSeqs, nResFreqs, resFreqs, resFreqs, mng);
  String spFileName = mng.spectrumFileName;
  String inputFileName = BioIo.getBaseName(spFileName) + "."
      + mng.inputFileExt;
  prsms = PrSMReader.readPrSM(inputFileName);
  */
}

}

    /* compute E-value. Separate: compute E-value separately or not */
/*
    public void process(boolean isSeparate) throws Exception {
        String spFileName = mng.spectrumFileName;
        File spFile = new File(spFileName);

        int nSpectra = MsAlignReader.countSpNum(spFile);
        MsAlignReader spReader = new MsAlignReader(spFile);
        String outputFileName = BioIo.getBaseName(spFileName) + "."
                + mng.outputFileExt;
        PrSMXmlWriter writer = new PrSMXmlWriter(new File(outputFileName));
        int cnt = 0;
        Ms<DeconvPeak>[] deconvSp;
        long startTime = System.currentTimeMillis();
        while ((deconvSp = spReader.getNextMses()) != null) {
            cnt++;
            for (int i = 0; i < deconvSp.length; i++) {
                String scan = deconvSp[i].getHeader().getScansString();
                String msg = ProcessUtil.updateMsg("E-value computation", scan,
                        nSpectra, cnt, startTime);
                System.out.print(msg);
                processOneSpectrum(deconvSp[i], isSeparate, writer);
            }
        }
        spReader.close();
        writer.close();
    }

    public void processOneSpectrum(Ms<DeconvPeak> deconvSp, boolean isSeparate,
            PrSMXmlWriter writer) throws Exception {
        SpectrumSet spectrumSet = SpectrumSet.getSpectrumSet(deconvSp, 0,
                mng.spPara);
        if (spectrumSet != null) {

            // System.out.println("prsm number " + prsms.length);
            PrSM selectedPrsms[] = PrSMUtil.getPrsms(prsms,
                    deconvSp.getHeader());
            PrSMUtil.processPrSM(selectedPrsms, deconvSp, seqs);
            if (isSeparate) {
                for (int j = 0; j < selectedPrsms.length; j++) {
                    compP.setPValue(deconvSp, selectedPrsms[j]);
                }
            } else {
                compP.setPValueArray(spectrumSet.getSpSix(),
                        selectedPrsms);
            }
            Arrays.sort(selectedPrsms, new EValueComparator());
            writer.write(selectedPrsms);
        }
    }
}
*/
