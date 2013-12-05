package edu.ucsd.proteomics.search.zeroptmsearch;

import java.io.File;

import com.jap.proteomics.base.seqdb.BpSpecReader;
import com.jap.proteomics.base.theosp.NModBpSpec;
import com.jap.proteomics.base.util.BioIo;
import com.jap.proteomics.base.util.ProcessUtil;
import com.jap.proteomics.spec.deconvsp.DeconvPeak;
import com.jap.proteomics.spec.deconvsp.reader.MsAlignReader;
import com.jap.proteomics.spec.sp.Ms;
import com.jap.proteomics.spec.spset.SpectrumSet;

import edu.ucsd.proteomics.prsm.base.PrSM;
import edu.ucsd.proteomics.prsm.base.SemiAlignType;
import edu.ucsd.proteomics.prsm.writer.PrSMXmlWriter;

public class ZeroPtmProcessor {
    private ZeroPtmMng mng;
    
    private NModBpSpec seqs[];
    private ZeroPtmSearcher searcher;

    public ZeroPtmProcessor(ZeroPtmMng mng) throws Exception {
        this.mng = mng; 
        seqs = BpSpecReader.getNModBpSpecs(mng.searchDbFileName,
                mng.resFileName, mng.allowProtNTermMods);
        searcher = new ZeroPtmSearcher(seqs, mng);
    }

    public void process() throws Exception {
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
    }

}
