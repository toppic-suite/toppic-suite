---
name: bruker
description: Process Bruker TOF top-down data with the TopPIC Suite - timsTOF / timsOmni .d folders (analysis.tsf or analysis.tdf) and maXis / impact Q-TOF .d folders (analysis.baf): identify the format, convert it to mzML on Linux (TopFD reads only mzML / mzXML), choose the TopFD settings (frame summing, -o, -a, -c) and prepare the msalign for TopPIC. Use when asked to deconvolute or search Bruker data, to convert a .d folder or a .baf / .tsf / .tdf file, or when TopFD / TopPIC results on Bruker data look wrong.
---

# Bruker TOF data in the TopPIC Suite

TopFD reads only mzML / mzXML, and the vendored ProteoWizard (`ext/pwiz`) has no
Bruker reader on Linux (ProteoWizard's Bruker support needs Bruker's libraries under
Windows). Bruker data therefore goes through a conversion; this skill covers the
formats, the converter in this folder (`tsf2mzml.py`, `tsfdata.py`) and the TopFD and
TopPIC settings measured on Bruker TOF data (toppic_test experiments 148 and 149).

## 1. What is in the .d folder

| file in `<run>.d/` | format | instruments | reader |
|---|---|---|---|
| `analysis.tsf` + `analysis.tsf_bin` | TSF: spectra without ion mobility separation | timsTOF / timsOmni with TIMS off, MALDI | Bruker timsdata library (`tsf_*`), `tsf2mzml.py` here |
| `analysis.tdf` + `analysis.tdf_bin` | TDF: frames of TIMS scans (PASEF) | timsTOF with TIMS | Bruker timsdata library (`tims_*`); no converter here |
| `analysis.baf` (+ `.baf_idx`, `.baf_xtr`) | BAF | maXis, impact, compact Q-TOF (otofControl) | Bruker baf2sql (needs the index files); toppic_test `bafread.py` without them |

Archives from Box often come as `<run>.d.7z`; without a `7z` binary use py7zr
(`pip install py7zr`, then `py7zr.SevenZipFile(path).extractall(dir)`).

TSF and TDF metadata are SQLite databases; look at them before converting:

```sh
sqlite3 <run>.d/analysis.tsf "select Key, Value from GlobalMetadata"   # InstrumentName, MzAcqRange*, HasLineSpectra, HasProfileSpectra
sqlite3 <run>.d/analysis.tsf "select ScanMode, MsMsType, count(*), min(Time), max(Time) from Frames group by 1, 2"
sqlite3 <run>.d/analysis.tsf "select TriggerMass, IsolationWidth, PrecursorCharge, CollisionEnergy, count(*) from FrameMsMsInfo group by 1, 2, 3, 4"
```

`MsMsType` 0 = MS1 frames; FrameMsMsInfo holds each MS/MS frame's precursor (trigger
m/z, isolation width, often no charge). The fragmentation method is not in TSF: for
a timsOmni read the method folder (`<run>.d/<method>.m/omnitrapApplication.bdl` is a
zip of zips with JSON "bundles"; an "eXd Process" bundle with `EI_SOURCE` and
`electron-source-configuration.json` electron energy 0 eV means ECD).

## 2. Converting TSF runs: `tsf2mzml.py`

Bruker's timsdata library (TDF-SDK) reads TSF; the PyPI package
**opentims-bruker-bridge** ships its Linux and Windows builds:

```sh
pip install opentims-bruker-bridge numpy        # libtimsdata.so / timsdata.dll
# or point to a copy: export TIMSDATA_LIB=/path/to/libtimsdata.so
python3 .claude/skills/bruker/tsf2mzml.py <run>.d <out>.mzML --mode profsum --sum 100 --activation ECD
```

`tsfdata.py` finds the library through `TIMSDATA_LIB`, the opentims_bruker_bridge
package, or `libtimsdata.so` next to it. Modes of `tsf2mzml.py`:

- `--mode line`: one spectrum per frame, Bruker's line (centroided) spectrum.
- `--mode linesum --sum N`: up to N consecutive MS/MS frames of one precursor, line
  spectra merged (peaks within `--ppm`, default 5, joined).
- `--mode profsum --sum N`: up to N consecutive MS/MS frames of one precursor, the
  profile spectra summed on the TOF index grid and centroided by the script.
  **Use this one** for repeated MS/MS of a precursor (infusion, or long PRM-like
  acquisition): from about 100 frames on it beats merging Bruker's line spectra.
- `--sum 0`: all consecutive frames of a precursor in one spectrum. A different
  precursor (trigger m/z, isolation width, collision energy) or an MS1 frame closes a
  group; MS1 frames are written one per frame.
- `--activation CID|HCD|ECD|ETD` (default CID) is written into the mzML, since TSF
  does not record it.

The MS/MS path is tested (experiment 148: byte-identical to the experiment's
conversions); the MS1 path (TSF runs with MS1 frames) has not been run on data yet,
so check its first output (spectra counts per level, `topfd` parameter block
"Number of MS1 scans").

Single frames of a timsOmni are weak (median line intensity about 320 against a
maximum about 44,000): in experiment 148 single frames gave the fewest real fragments,
and among each run's 500 most intense fragment masses only 7-24 were fragments,
against 65-76 for the sum of all frames.

## 3. TopFD settings for Bruker TOF MS/MS

From experiment 148 (timsOmni ECD, one infused precursor, no MS1; scored by fragment
masses matching the known proteoform above the matches of decoy ion sets, and by
TopPIC):

- **`-a ETD` for ECD data.** TopFD has no ECD activation; ETD has the same ion types
  (c, z-dot). With `--activation ECD` in the mzML, pass `-a ETD` (not `-a FILE`).
- **`-o` when the run has no MS1 spectra.** TopFD then writes the `-m` / `-c` limits as
  every spectrum's PRECURSOR_MASS / PRECURSOR_CHARGE (50,000 Da, 30 by default).
- **With `-o`, set `-c` to the precursor's charge + 1.** The `-c` limit bounds the
  fragment charges: with the default 30, 1,011 of the 1,657 masses that passed
  TopPIC's EnvCNN filter were fitted at charges 11-30, impossible for a 9+ precursor,
  and matched only at random. `-c 10` removed them and TopPIC's E-value fell from
  8.3e-8 to 9.0e-29. Do not lower `-m` to the precursor mass (it lost fragments).
  ECD spectra show the precursor's charge directly: the charge-reduced species
  (M+zH)^(z-1)+. , (z-2)+. ... dominate the spectrum; the largest deconvoluted mass
  is the precursor (23,466.61 Da at 8+, 7+, 9+ in experiment 148).
- **Other parameters matter little.** `-e` (0.005-0.05), `-s` (0.5-5) and
  `--mass-group-boundary` (2000-8000) changed the real fragments by about 10% at most
  (`-c 10 -s 5` gave the best coverage); `-n` only reorders the masses;
  `--max-miss-peak-num` and `--disable-filter-by-mz` changed nothing with the greedy
  MS/MS deconvolution (the default), which does not read `--ms2-rank-peak-distance` or
  `--disable-ms2-harmonic-check` (those act with `-q msdeconv`).
- **The noise filters trade fragments for purity:** `-v 0.1` kept 12% of the masses
  and 70% of the fragments, `-g` 3% and 47%, `-q msdeconv` 17% and 80%,
  `-q flashdeconv` 0.3% and 41% (no random match; the best TopPIC E-values when the
  sequence is in the database, too few masses for a search with unknown shifts).
  TopPIC filters masses at EnvCNN 0.4 itself (`-F`), so `topfd -v` below 0.4 does not
  change TopPIC's results.
- Mass accuracy of the fragments: median +0.4 to +0.7 ppm, interquartile range about
  3 ppm, so `-e 0.02` and TopPIC's 10 ppm tolerance are adequate.

Recommended for a timsOmni ECD run of one precursor:

```sh
python3 .claude/skills/bruker/tsf2mzml.py run.d run.mzML --mode profsum --sum 0 --activation ECD
topfd -o -a ETD -c 10 -s 5 run.mzML                 # -c = precursor charge + 1
```

(`--sum 100`, about 12 s of acquisition per spectrum, gave the best single TopPIC
E-values for a search with two unknown shifts, 16 of 20 spectra at the TopFD defaults.)

## 4. TopPIC on spectra deconvoluted with `-o`

- Replace the msalign's precursor before searching, or TopPIC searches a 50,000 Da
  precursor:

  ```sh
  sed -e 's/^PRECURSOR_MASS=.*/PRECURSOR_MASS=23466.611/' -e 's/^PRECURSOR_CHARGE=.*/PRECURSOR_CHARGE=9/' \
      -e 's/^PRECURSOR_MZ=.*/PRECURSOR_MZ=2608.408/' run_ms2.msalign > run_fixed_ms2.msalign
  ```

- `-a ETD` for ECD; `-x` (no feature file without MS1) and `-D` (no post mass matching,
  which reads the TopFD SQLite database by spectrum) when the msalign is edited.
- Shift ranges: a tag or a bound ligand can exceed TopPIC's default +-500 Da (`-m` /
  `-M`); in experiment 148 a 1,910 Da His tag and GDP (443 Da) needed `-s 2 -M 5000`
  against the bare sequence, or the construct's sequence in the database and `-s 1`.
- Native ECD fragments mostly cover the termini; in experiment 148 the c ions covered
  the tag and residues 1-30 and the z ions about 100-185 of KRAS, so shifts in the core
  are placed only roughly.

## 5. TDF and BAF

- **TDF (timsTOF with TIMS)**: the same library reads it (`tims_open`,
  `tims_read_scans_v2`, `tims_index_to_mz`, `tims_scannum_to_oneoverk0`); a converter
  would sum each MS/MS frame's scans (or a PASEF precursor's scans) into one spectrum
  before writing the mzML. Not written or tested yet; ProteoWizard msconvert under
  Windows (`--combineIonMobilitySpectra`) is the alternative.
- **BAF (maXis / impact)**: Bruker's baf2sql library (Linux build from the pyBaf2Sql
  project) needs the `.baf_idx` / `.baf_xtr` files of the `.d` folder ("Error reading
  BAF files" without them; PRIDE often holds only the `.baf`). toppic_test experiment
  149 (`tests/202610/scripts/bruker_qtof_deconv/`: `bafread.py`, `baf2mzml.py`,
  `check_bafread.py`) reads a lone `.baf` directly (Bruker's line spectra, the
  acquisition calibration; profile spectra are VLI2-compressed and not read) and can
  merge the MS/MS spectra of one precursor (`--merge precursor`); its settings study
  was in progress when this skill was written - see its report 0149.
- Tabb et al., J. Proteome Res. 2023, 22, 2199 found on maXis II data that TopFD on
  one mzML spectrum per MS/MS scan impedes identification and that summing replicated
  MS/MS scans of a precursor gains signal-to-noise, as experiment 148 found for the
  timsOmni.

## References

- toppic_test report 0148 (`tests/202610/reports/0148_10bf4788_20261008_tims_tof_deconv.md`):
  timsOmni ECD of KRAS G12C, the conversion, the TopFD grid, TopPIC.
- toppic_test instruction 149 (`tests/202610/instructions/149_53ffa694_bruker_qtof_deconv.txt`):
  maXis II BAF, PXD010825.
- Bruker TDF-SDK via opentims-bruker-bridge (PyPI); `tsfdata.py` uses `tsf_open`,
  `tsf_read_line_spectrum_v2`, `tsf_read_profile_spectrum_v2`, `tsf_index_to_mz`.
