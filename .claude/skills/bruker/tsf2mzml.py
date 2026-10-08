#!/usr/bin/env python3
"""tsf2mzml.py <.d folder> <out.mzML> [--mode line|linesum|profsum] [--sum N]
                [--activation CID|HCD|ECD|ETD] [--first F] [--last L] [--ppm P]

Converts a Bruker TSF run (timsTOF / timsOmni without ion mobility separation:
analysis.tsf + analysis.tsf_bin) to a centroid mzML for TopFD, with Bruker's
timsdata library (tsfdata.py):

  line     one spectrum per frame: Bruker's line (centroided) spectrum.
  linesum  up to N consecutive MS/MS frames of one precursor per spectrum:
           their line spectra merged, peaks within --ppm (default 5) of each
           other joined (intensity summed, m/z intensity-weighted).
  profsum  up to N consecutive MS/MS frames of one precursor per spectrum:
           their profile spectra summed on the TOF index grid and centroided
           here (each local maximum with the neighbouring points at or above
           half its height, at most 8 per side; m/z intensity-weighted,
           intensity the points' sum).

--sum 0 sums all consecutive frames of a precursor. MS/MS frames of one
precursor = the same trigger m/z, isolation width and collision energy in
FrameMsMsInfo; a different precursor or an MS1 frame closes a group. MS1
frames are written one per frame (line spectrum, or the frame's profile
centroided with profsum). TSF does not record the fragmentation method:
--activation sets it (default CID; ECD for the Omnitrap's ECD). The precursor
charge is written when FrameMsMsInfo has one. Ids are scan=<k> (k = 1, 2,
...) and the spectrum title lists the frames.

Tested on a timsOmni ECD run (1,976 MS/MS frames of one precursor, no MS1
frame; toppic_test experiment 148). The MS1 path has not been run on data.
"""
import argparse
import base64
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from tsfdata import TsfData  # noqa: E402


def merge_lines(spectra, ppm):
    mz = np.concatenate([s[0] for s in spectra])
    it = np.concatenate([s[1] for s in spectra])
    order = np.argsort(mz)
    mz, it = mz[order], it[order]
    out_mz, out_it = [], []
    i = 0
    n = len(mz)
    while i < n:
        j = i + 1
        while j < n and (mz[j] - mz[j - 1]) <= mz[j - 1] * ppm * 1e-6:
            j += 1
        w = it[i:j]
        out_it.append(w.sum())
        out_mz.append((mz[i:j] * w).sum() / w.sum())
        i = j
    return np.array(out_mz), np.array(out_it)


def centroid_profile(tsf, frame, prof, reach=8):
    """Centroids of a profile spectrum on the TOF index grid: each local
    maximum with the points on both sides that stay at or above half its
    height and do not rise again (at most `reach` points per side; a TOF peak
    is about 2-8 points wide at half height over m/z 200-5500); m/z the
    intensity-weighted mean index converted to m/z, intensity the points'
    sum."""
    p = prof.astype(np.float64)
    n = len(p)
    inner = (p[1:-1] >= p[:-2]) & (p[1:-1] > p[2:]) & (p[1:-1] > 0)
    apex = np.nonzero(inner)[0] + 1
    top = p[apex]
    half = top / 2
    wsum = top.copy()
    isum = apex * top
    for side in (-1, 1):
        prev = top.copy()
        alive = np.ones(len(apex), dtype=bool)
        for k in range(1, reach + 1):
            j = apex + side * k
            inside = (j >= 0) & (j < n)
            v = np.where(inside, p[np.clip(j, 0, n - 1)], 0.0)
            ok = alive & inside & (v >= half) & (v <= prev)
            wsum += np.where(ok, v, 0.0)
            isum += np.where(ok, j * v, 0.0)
            prev = np.where(ok, v, prev)
            alive = ok
    return tsf.index_to_mz(frame, isum / wsum), wsum


ACTIVATION = {
    "CID": ("MS:1000133", "collision-induced dissociation"),
    "HCD": ("MS:1000422", "beam-type collision-induced dissociation"),
    "ECD": ("MS:1000250", "electron capture dissociation"),
    "ETD": ("MS:1000598", "electron transfer dissociation"),
}


def b64(a, dtype):
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode()


def binary_arrays(mz, it):
    return ('    <binaryDataArrayList count="2">\n'
            f'     <binaryDataArray encodedLength="{len(b64(mz, np.float64))}">\n'
            '      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>\n'
            '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>\n'
            '      <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n'
            f'      <binary>{b64(mz, np.float64)}</binary>\n     </binaryDataArray>\n'
            f'     <binaryDataArray encodedLength="{len(b64(it, np.float32))}">\n'
            '      <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>\n'
            '      <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>\n'
            '      <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of detector counts"/>\n'
            f'      <binary>{b64(it, np.float32)}</binary>\n     </binaryDataArray>\n'
            '    </binaryDataArrayList>\n')


def write_mzml(path, spectra, src_name, activation="CID"):
    """spectra: list of dict(level, mz, it, rt (s), frames (first, last), target, width, charge)."""
    with open(path, "w") as f:
        w = f.write
        w('<?xml version="1.0" encoding="utf-8"?>\n')
        w('<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">\n')
        w(' <cvList count="2">\n'
          '  <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.182" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>\n'
          '  <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>\n'
          ' </cvList>\n')
        w(' <fileDescription>\n  <fileContent>\n'
          + ('   <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>\n'
             if any(s["level"] == 1 for s in spectra) else "") +
          '   <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>\n'
          '   <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>\n'
          '  </fileContent>\n  <sourceFileList count="1">\n'
          f'   <sourceFile id="SF1" name="{src_name}" location="file:///">\n'
          '    <cvParam cvRef="MS" accession="MS:1000776" name="scan number only nativeID format" value=""/>\n'
          '   </sourceFile>\n  </sourceFileList>\n </fileDescription>\n')
        w(' <softwareList count="1">\n  <software id="tsf2mzml" version="1.0">\n'
          '   <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="tsf2mzml.py"/>\n'
          '  </software>\n </softwareList>\n')
        w(' <instrumentConfigurationList count="1">\n  <instrumentConfiguration id="IC1">\n'
          '   <cvParam cvRef="MS" accession="MS:1000031" name="instrument model" value=""/>\n'
          '  </instrumentConfiguration>\n </instrumentConfigurationList>\n')
        w(' <dataProcessingList count="1">\n  <dataProcessing id="DP1">\n'
          '   <processingMethod order="1" softwareRef="tsf2mzml">\n'
          '    <cvParam cvRef="MS" accession="MS:1000035" name="peak picking" value=""/>\n'
          '   </processingMethod>\n  </dataProcessing>\n </dataProcessingList>\n')
        w(' <run id="run1" defaultInstrumentConfigurationRef="IC1">\n')
        w(f'  <spectrumList count="{len(spectra)}" defaultDataProcessingRef="DP1">\n')
        for k, s in enumerate(spectra):
            mz, it = s["mz"], s["it"]
            n = len(mz)
            w(f'   <spectrum index="{k}" id="scan={k + 1}" defaultArrayLength="{n}">\n')
            if s["level"] == 1:
                w('    <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>\n'
                  '    <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>\n')
            else:
                w('    <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>\n'
                  '    <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>\n')
            w('    <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>\n'
              '    <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>\n')
            if n:
                w(f'    <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{it.sum():.1f}"/>\n')
            w(f'    <cvParam cvRef="MS" accession="MS:1000796" name="spectrum title" value="frames {s["frames"][0]}-{s["frames"][1]}"/>\n')
            w('    <scanList count="1">\n'
              '     <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>\n'
              '     <scan>\n'
              f'      <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{s["rt"] / 60.0:.6f}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>\n'
              '     </scan>\n    </scanList>\n')
            if s["level"] == 1:
                w(binary_arrays(mz, it))
                w('   </spectrum>\n')
                continue
            half = s["width"] / 2.0
            w('    <precursorList count="1">\n     <precursor>\n      <isolationWindow>\n'
              f'       <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{s["target"]:.4f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n'
              f'       <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{half:.4f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n'
              f'       <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{half:.4f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n'
              '      </isolationWindow>\n      <selectedIonList count="1">\n       <selectedIon>\n'
              f'        <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{s["target"]:.4f}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>\n'
              + (f'        <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{s["charge"]}"/>\n'
                 if s.get("charge") else "") +
              '       </selectedIon>\n      </selectedIonList>\n      <activation>\n'
              f'       <cvParam cvRef="MS" accession="{ACTIVATION[activation][0]}" name="{ACTIVATION[activation][1]}" value=""/>\n'
              '      </activation>\n     </precursor>\n    </precursorList>\n')
            w(binary_arrays(mz, it))
            w('   </spectrum>\n')
        w('  </spectrumList>\n </run>\n</mzML>\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("d")
    ap.add_argument("out")
    ap.add_argument("--mode", choices=["line", "linesum", "profsum"], default="line")
    ap.add_argument("--sum", type=int, default=1)
    ap.add_argument("--activation", choices=sorted(ACTIVATION), default="CID")
    ap.add_argument("--first", type=int, default=None)
    ap.add_argument("--last", type=int, default=None)
    ap.add_argument("--ppm", type=float, default=5.0)
    a = ap.parse_args()
    tsf = TsfData(a.d)
    rows = tsf.conn.execute(
        "select f.Id, f.Time, f.MsMsType, m.TriggerMass, m.IsolationWidth, m.PrecursorCharge, "
        "m.CollisionEnergy from Frames f left join FrameMsMsInfo m on m.Frame = f.Id "
        "order by f.Id").fetchall()
    if a.first is not None:
        rows = [r for r in rows if r[0] >= a.first]
    if a.last is not None:
        rows = [r for r in rows if r[0] <= a.last]
    # groups: one MS1 frame, or consecutive MS/MS frames of one precursor (at most --sum)
    limit = 1 if a.mode == "line" else (0 if a.sum == 0 else a.sum)
    groups = []
    for r in rows:
        ms2 = r[2] != 0 and r[3] is not None
        key = (r[3], r[4], r[6]) if ms2 else None
        g = groups[-1] if groups else None
        if (ms2 and g is not None and g["key"] == key
                and (limit == 0 or len(g["rows"]) < limit)):
            g["rows"].append(r)
        else:
            groups.append({"key": key, "rows": [r]})
    spectra = []
    for g in groups:
        grp = g["rows"]
        frames = [r[0] for r in grp]
        level = 2 if g["key"] is not None else 1
        if a.mode == "line" or (a.mode == "linesum" and level == 1):
            mz, it = tsf.line_spectrum(frames[0])
        elif a.mode == "linesum":
            mz, it = merge_lines([tsf.line_spectrum(fr) for fr in frames], a.ppm)
        else:
            prof = None
            for fr in frames:
                p = tsf.profile_spectrum(fr).astype(np.uint64)
                prof = p if prof is None else prof + p
            mz, it = centroid_profile(tsf, frames[0], prof)
        spectra.append(dict(level=level, mz=mz, it=it, rt=float(np.mean([r[1] for r in grp])),
                            frames=(frames[0], frames[-1]), target=grp[0][3], width=grp[0][4],
                            charge=grp[0][5]))
    write_mzml(a.out, spectra, os.path.basename(os.path.normpath(a.d)), a.activation)
    for level in (1, 2):
        npk = [len(s["mz"]) for s in spectra if s["level"] == level]
        if npk:
            print(f"{a.out}: MS{level} {len(npk)} spectra, peaks per spectrum min {min(npk)} "
                  f"median {int(np.median(npk))} max {max(npk)}")
    tsf.close()


if __name__ == "__main__":
    main()
