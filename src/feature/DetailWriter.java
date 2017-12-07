/**
 * Export Match Envelopes
 *
 * @author  Xiaowen Liu
 * @date    2009-10-9
 */

package edu.ucsd.proteomics.msdeconv.writer;

import java.io.PrintWriter;

import com.jap.proteomics.spec.sp.MsHeader;

import edu.ucsd.proteomics.msdeconv.env.Env;
import edu.ucsd.proteomics.msdeconv.env.MatchEnv;
import edu.ucsd.proteomics.msdeconv.env.RealEnv;

public class DetailWriter {

    public static void writeEnv(MatchEnv env, PrintWriter out) {
        Env theoEnv = env.getTheoEnv();
        RealEnv realEnv = env.getRealEnv();
        out.println();
        out.println("BEGIN ENVELOPE");
        out.println("REF_IDX " + theoEnv.getReferIdx());
        out.println("CHARGE " + theoEnv.getCharge());
        out.println("SCORE " + env.getScore());
        out.println("THEO_PEAK_NUM " + theoEnv.getPeakNum() + " REAL_PEAK_NUM "
                + (realEnv.getPeakNum() - realEnv.getMissPeakNum()));
        out.println("THEO_MONO_MZ " + theoEnv.getMonoMz() + " REAL_MONO_MZ "
                + realEnv.getMonoMz());
        out.println("THEO_MONO_MASS " + theoEnv.getMonoMass()
                + " REAL_MONO_MASS " + realEnv.getMonoMass());
        out.println("THEO_INTE_SUM " + theoEnv.compIntensitySum()
                + " REAL_INTE_SUM " + realEnv.compIntensitySum());
        for (int i = 0; i < theoEnv.getPeakNum(); i++) {
            out.println(theoEnv.getMz(i) + " " + theoEnv.getIntensity(i) + " "
                    + realEnv.isExist(i) + " " + realEnv.getPeakIdx(i) + " "
                    + realEnv.getMz(i) + " " + realEnv.getIntensity(i));
        }
        out.println("END ENVELOPE");
    }

    /**
     * Write an MatchEnv List to a text file.
     */
    public static void writeEnv(PrintWriter out, MatchEnv envs[],
            MsHeader header) throws Exception {
        out.println("BEGIN SPECTRUM");
        out.println("ID " + header.getId());
        out.println("SCANS " + header.getScansString());
        out.println("Ms_LEVEL " + header.getMsLevel());
        out.println("ENVELOPE_NUMBER " + envs.length);
        out.println("MONOISOTOPIC_MASS " + header.getPrecMonoMass());
        out.println("CHARGE " + header.getPrecChrg());
        for (int i = 0; i < envs.length; i++) {
            writeEnv(envs[i], out);
        }
        out.println("END SPECTRUM");
    }
}
