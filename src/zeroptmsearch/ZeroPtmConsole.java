package edu.ucsd.proteomics.search.zeroptmsearch;

import com.jap.proteomics.base.ion.EnumActivation;
import com.jap.proteomics.base.util.BioIo;

public class ZeroPtmConsole {
    public static void main(String[] args) {
        try {
            ZeroPtmVersion version = new ZeroPtmVersion();
            System.out.println("ZeroPtmConsole " + version.getVersion() + " "
                    + version.getDate());
            ZeroPtmMng mng = new ZeroPtmMng();
            mng.searchDbFileName = args[0];
            mng.resFileName = args[1];
            mng.spectrumFileName = args[2];
            mng.outputFileExt = args[3];
            mng.spPara.setActivationType(EnumActivation.getActivationType(args[4]));
            ZeroPtmProcessor processor = new ZeroPtmProcessor(mng);
            processor.process();
        } catch (Exception e) {
            System.out.println("\n\n" + BioIo.getStackTraceAsString(e));
            System.exit(1);
        } catch (Error e) {
            System.out.println("\n\n" + BioIo.getStackTraceAsString(e));
            System.exit(1);
        }
    }
}
