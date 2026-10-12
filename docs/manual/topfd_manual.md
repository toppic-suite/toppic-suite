# TopFD manual

TopFD (Top-down mass spectral Feature Detection) is the spectral
deconvolution tool of the TopPIC Suite. It groups the peaks of top-down mass
spectra into isotopic envelopes and converts the envelopes into monoisotopic
neutral masses, and it extracts proteoform features from LC-MS or CE-MS
data. Its output is the input of TopPIC and TopMG.

TopFD has two main functions:

1. **Deconvoluting mzML/mzXML files** (the normal use): a whole LC-MS/MS run
   is deconvoluted, MS1 spectra are used to detect proteoform features, and
   MS/MS spectra are written as `msalign` files for database search.
2. **Deconvoluting a peak list file**: a single centroided MS/MS spectrum,
   given as a text file of m/z and intensity pairs, is deconvoluted on its
   own. This is meant for testing and for inspecting the deconvolution of one
   spectrum in detail.

Both are run from the command line with the `topfd` executable (`topfd.exe`
on Windows). The general form is

```sh
topfd [options] spectrum-file ...
```

`topfd -h` prints the option list. The GUI tool `topfd_gui` offers both
functions: the second one with the checkbox "Text peak list input" under
"Additional settings" (`-T`).

**Paths must not contain spaces.** TopFD refuses to run if the directory of
the executable contains a space, and input files whose path contains spaces
or quotation marks are not found. Put the data in a directory without spaces
(e.g. `C:\data\` rather than `C:\My Documents\`).

## 1. Deconvoluting mzML/mzXML files

### 1.1 Input

TopFD reads **centroided** mzML or mzXML files (the extension must be
`.mzML`, `.mzml`, `.mzXML` or `.mzxml`; other files are skipped with the
message "is not a valid mass spectral file"). Vendor raw files must be
converted first, e.g. with ProteoWizard's `msconvert`, using peak picking so
the spectra are centroided:

```sh
msconvert --mzML --filter "peakPicking true 1-" sample.raw
```

Several files can be given on one command line; they are processed one after
another with the same parameters, and each gets its own set of output files.

FAIMS data (one file with several compensation voltages) is recognised
automatically. Each voltage level is deconvoluted as a separate fraction and
the integer voltage is appended to the output file names (see 1.3).

### 1.2 Basic usage

```sh
# default parameters: fragmentation method taken from the file,
# max charge 30, max mass 50,000 Da, m/z error 0.02, 1 thread
topfd sample.mzML

# several files, 8 threads
topfd -u 8 sample_1.mzML sample_2.mzML
```

For each file TopFD prints the parameters it uses, then runs three steps:

1. **MS1 deconvolution**: every MS1 spectrum is deconvoluted with a greedy
   search: the charge state of the most intense peak is decided with all the
   charge states of its mass, the mass's envelopes are subtracted, and the
   search repeats on the rest of the spectrum. With `--deconv-ms1 msdeconv`
   the MS-Deconv dynamic programming, the default method of earlier
   versions, is used instead. With either method, envelopes whose EnvCNN
   score is below `--envcnn-cutoff-ms1` (`-f`, default 0.2) are not
   reported.
2. **MS1 feature detection**: the deconvoluted MS1 envelopes are linked
   across scans into proteoform features (LC-MS features), each scored with
   the ECScore neural network. Every MS/MS spectrum is then matched to the
   feature(s) in its precursor isolation window, which gives its precursor
   monoisotopic mass, charge and intensity.
3. **MS/MS deconvolution**: every MS/MS spectrum is deconvoluted, using its
   precursor information from step 2, with the greedy search of step 1
   adapted to fragment ions. With `--deconv-ms2 msdeconv` the MS-Deconv
   dynamic programming, the default method of earlier versions, is used
   instead, and with `--deconv-ms2 flashdeconv` the FLASHDeconv algorithm
   (experimental). With `greedy` and `msdeconv`, envelopes whose EnvCNN
   score is below `--envcnn-cutoff-ms2` (`-g`, default 0: all reported) are
   not reported.

If the file has no MS1 spectra, pass `--missing-level-one` (`-o`): steps 1 and
2 are skipped and no precursor is determined. Every MS/MS spectrum is then
deconvoluted up to the maximum mass and charge given by `-m` and `-c`, and
its `PRECURSOR_MASS` and `PRECURSOR_CHARGE` are written as those limits, with
feature id -1 and intensity 0. Only the isolation window is taken from the
file.

### 1.3 Output files

All output files are written next to the input file, with the input file
name minus its extension as the base name (`sample` for `sample.mzML`; for
FAIMS data `sample_<voltage>` per voltage level, e.g. `sample_-40`).

| File | Contents |
|---|---|
| `sample_ms1.msalign` | Deconvoluted MS1 spectra (monoisotopic masses). |
| `sample_ms2.msalign` | Deconvoluted MS/MS spectra with precursor information (monoisotopic masses; average masses with `-A`). **This is the input of TopPIC and TopMG.** |
| `sample_ptcr.msalign` | Only when the file has PTCR spectra (see below): their deconvoluted species, one line per species, in the `msalign` format with `ACTIVATION=PTCR`. Not read by TopPIC or TopMG. |
| `sample_ptcr.tsv` | The same species with their charge ladders: one line per species with its spectrum, monoisotopic mass, intensity, charges, the charges at which it could have been isolated, the intensity and EnvCNN score per charge, whether it is the precursor, and the MS1 feature it matches (`Feature ID`, `Feature mass`; -1 and empty without one). |
| `sample_ms3.msalign` | Only when the file has MS3 scans (the dissociation of MS2 fragment ions; see below): the deconvoluted MS3 spectra in the `msalign` format with `LEVEL=3`, the parent MS2 spectrum (`MS_TWO_ID`, `MS_TWO_SCAN`, its precursor as `MS_TWO_PRECURSOR_*`) and the MS3 precursor (the MS2 fragment ion). |
| `sample_ms3.tsv` | One line per MS3 scan: its MS2 and MS1 scans, the MS2 precursor, the isolation window, how the precursor was determined (`ms2_species`: a species of the parent MS2 spectrum's deconvolution; `selected_ion`: the file's selected ion at every charge up to the MS2 precursor's), the precursor and the number of masses. |
| `sample_ms1.feature` | Proteoform features detected in the LC-MS map: one line per feature with its mass, intensity, retention-time and scan range, charge range, apex and ECScore. |
| `sample_ms2.feature` | The features assigned to each MS/MS spectrum (precursor mass, m/z, charge and intensity), one line per spectrum and feature. |
| `sample_feature.xml` | The proteoform features in XML, with their per-charge envelope information; read by TopMG when it combines fractions (`-c`). |
| `sample.sqlite` | An SQLite database with the deconvoluted MS1 and MS/MS spectra and their peaks, for spectrum visualisation and for TopPIC's post mass matching (see the `-N` option). Not written with `--no-sql`. With `--sql-3d` it also holds the raw MS1 peak tables for 3D visualisation (`CONFIG` and `PEAKS0`, `PEAKS1`, ...). |

The `msalign` format is a text format. Each spectrum is a `BEGIN IONS` ...
`END IONS` block with header lines (`FILE_NAME`, `MSALIGN_FILE_NAME` (the
msalign file itself, per voltage level for FAIMS data), `SPECTRUM_ID`,
`TITLE`, `SCANS`, `RETENTION_TIME` in minutes, `LEVEL`, and for MS/MS spectra
`MS_ONE_ID`, `MS_ONE_SCAN`, `PRECURSOR_WINDOW_BEGIN/END`, `ACTIVATION`,
`PRECURSOR_MZ`, `PRECURSOR_CHARGE`, `PRECURSOR_MASS`, `PRECURSOR_INTENSITY`,
`PRECURSOR_FEATURE_ID`, `DECONVOLUTED_MASS_NUMBER`) followed by one line per
deconvoluted mass with four tab-separated columns: monoisotopic neutral
mass, intensity, charge and score. With `--average-mass` the MS/MS file
reports average masses instead: each fragment mass, `PRECURSOR_MASS` and
`PRECURSOR_MZ` are the intensity-weighted average mass (m/z) of the
theoretical isotopic distribution of the monoisotopic mass, from the
distribution table in `res/base_data/theo_patt.txt`; the parameter block
says `Reported precursor and fragment masses: Average`, and TopPIC reads
such a file with `--average-mass`. When several features fall in the
isolation window of an MS/MS spectrum, those with at least a tenth of the
strongest feature's intensity are all kept and the `PRECURSOR_*` values are
separated by `:`. The parameters used are recorded as `#` comment lines at
the top of the file, starting with the input file (`File name`) and the
msalign file's own name (`Msalign file name`, e.g. `sample_-40_ms1.msalign`
for a FAIMS voltage level). With `--missing-level-one` no `_ms1.msalign`
and no feature files are produced.

With more than one thread, per-thread partial `msalign` files
(`sample_ms1.msalign_0`, `sample_ms2.msalign_0`, `_1`, ...) exist while
TopFD runs; they are merged into the final files and deleted at the end.

**MS3 spectra.** An MS3 scan (the dissociation of an MS2 fragment ion,
e.g. CID of a peptide and then CID of one of its fragments) is deconvoluted
after the MS/MS spectra, when the file has such scans (not for FAIMS
data): its precursor is the species of the parent MS2 spectrum's
deconvolution whose isotopic envelope falls in the MS3 isolation window
(the most intense one, its monoisotopic mass checked against the MS2
spectrum's peaks one isotope below and above), or the selected ion of the
file at every possible charge when the MS2 spectrum has no
species there; the fragments are deconvoluted with the MS/MS method bounded
by that precursor's mass and charge and written to `sample_ms3.msalign` and
`sample_ms3.tsv`. The hidden option `--disable-ms3` turns the stage off.
TopPIC and TopMG do not read `sample_ms3.msalign` by themselves; TopPIC
can search it as an MS/MS file with `-x` and `--disable-post-match`.

**PTCR spectra.** An MS/MS spectrum whose activation in the input file is
proton transfer charge reduction (PSI-MS MS:1003249, or proton transfer
reaction MS:1003248) without a dissociation term is a PTCR spectrum: the
isolated ions were charge-reduced, not fragmented. Proton transfer only
removes charges, so every ion of such a spectrum is an intact species
isolated in the precursor window, seen at the charge it was isolated at and
at lower charges (a charge ladder). TopFD deconvolutes these spectra with
the charge-state search of its MS1 deconvolution, allowing a mass only at
charges at or below a charge at which its isotopic envelope falls in the
isolation window, reports one mass per species (its charges merged, the
mass refined from all its isotope peaks), and takes the precursor from the
species at the selected ion's charge and m/z: the species whose
monoisotopic m/z is the selected ion's (Thermo converters write the
instrument's monoisotopic m/z when it has one), else the most intense one
whose isotopic envelope, from its monoisotopic peak on, covers it. They are written to
`sample_ptcr.msalign` and `sample_ptcr.tsv`, not to `sample_ms2.msalign` or
the SQLite database. They are also read from files without MS1 spectra
(`-o`). In a run with MS1 spectra, a species whose mass matches an MS1
feature of the spectrum's precursor window (within 10 ppm, or one isotope
apart) is linked to it in `sample_ptcr.tsv`, and a precursor species with a
feature gets the feature's id and intensity. With `-a` set to a
fragmentation method the spectra are fragment PTCR spectra of that method
(next paragraph).

**Fragment PTCR spectra.** A spectrum with a dissociation and a later (or
the same) proton transfer stage, such as the Thermo MS3 scan `ms3
1402.50@cid60.00 1100.00@ptr10.00` (collision-induced dissociation of the
precursor, then proton transfer charge reduction of its fragments), is a
fragment PTCR spectrum. Its activation is the dissociation's (the activation terms of all
its precursors count; msconvert writes one precursor per stage, the
dissociation first). TopFD deconvolutes it as an MS/MS spectrum of that
activation and writes it to `sample_ms2.msalign` (`LEVEL=3` for an MS3
scan), so TopPIC and TopMG search it with that activation's ion types; the
console reports the number of such spectra. A fragment has no more charges
or mass than its precursor, so the fragments are bounded by the charge of
the selected ion and the mass it has at the upper end of the precursor
window, not by `-c` and `-m` (with MS1 features of the window, by the
larger of these and the features' bounds). Without an MS1 feature in the
precursor window, the selected ion (its m/z and charge as the file gives
them) is the precursor. An MSn scan (n > 2) is grouped with the last MS1
scan before it. MSn scans without a proton transfer stage (a dissociation
of fragments) are not read as MS/MS spectra; the MS3 ones are the MS3
spectra above. ThermoRawFileParser writes only the proton transfer stage
of an MS3 scan as its precursor; TopFD then takes the first
stage's isolation target and dissociation from the scan's filter string,
with the default precursor window (`-w`). The same holds for a scan without
any precursor (ThermoRawFileParser writes none when the scan's precursor
scan is not in the file, e.g. a single scan exported from a run): the
filter string gives its isolation target, default window and activation.

**Dissociation of charge-reduced precursors.** A spectrum whose
dissociation stage comes after its proton transfer stage, such as the
Thermo MS3 scan `ms3 942.55@ptr5.00 997.93@hcd80.00` (proton transfer
charge reduction of the ions isolated at m/z 942.55, then the isolation of
one charge-reduced ion at m/z 997.93 and its higher-energy collisional
dissociation), is not a fragment PTCR spectrum: it is an MS/MS spectrum of
the dissociation's activation (`LEVEL=3`), written to `sample_ms2.msalign`,
whose precursor window is the first stage's. The order of the stages comes
from the scan's filter string, else from the order of its precursors in the
file. Its precursor is the charge-reduced ion: TopFD deconvolutes the PTCR
spectrum of the same isolation read before it (the scan the MS3 scan
references, else the last PTCR scan with the same isolation target) and
takes its most intense species with an isotopic envelope at the second
stage's target, with the charge it has there. Without such a PTCR spectrum
the precursors are the MS1 features of the window with an envelope at that
target at a lower charge (all features of the window when none has one),
and without MS1 spectra the selected ion when the file gives its charge.
The console reports the number of these spectra.

### 1.4 Options

Parameters that apply to both MS1 and MS/MS deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-c`, `--max-charge <int>` | 30 | Maximum charge state of precursor and fragment ions. |
| `-m`, `--max-mass <number>` | 50000 | Maximum monoisotopic mass (Da) of precursor and fragment ions. |
| `-e`, `--mz-error <number>` | 0.02 | Error tolerance of peak m/z values (m/z units). |
| `-u`, `--thread-number <int>` | 1 | Number of threads. Must not exceed the number of hardware threads; TopFD warns (but still runs) when the available memory looks too small for the requested number. |
| `-o`, `--missing-level-one` | off | The file has no MS1 spectra: skip MS1 deconvolution and feature detection. |
| `-N`, `--no-sql` | off | Do not write the `.sqlite` database. **Do not use it if TopPIC will search the spectra**: TopPIC's post mass matching, which is on by default, reads the centroided MS/MS peaks from this database and stops with an error when it is missing (see the [TopPIC manual](toppic_manual.md); `toppic --disable-post-match` is the alternative). |
| `-A`, `--average-mass` | off | Report average masses instead of monoisotopic masses for the precursor and fragment ions in `sample_ms2.msalign` (section 1.3). The MS1 msalign file, the feature files and the SQLite database keep monoisotopic masses. Search such a file with `toppic --average-mass`. Also applies to `sample_ptcr.msalign`, `sample_ms3.msalign` and `-T`. |
| `-D`, `--sql-3d` | off | Also store the raw MS1 peaks for 3D visualisation in the `.sqlite` database: `PEAKS0` holds every MS1 peak and `PEAKS1`, `PEAKS2`, ... progressively down-sampled copies, with one `CONFIG` row per table. Cannot be combined with `--no-sql`; has no effect with `-T` or `-o`. In `topfd_gui`, the checkbox "MS1 3D visualization data" under "Additional settings" turns this on (it is greyed out while "Do not generate SQLite database" is checked). |
| `-T`, `--text-peak-list` | off | The input is a text peak list (one MS/MS spectrum), not an mzML file; lines that do not start with a number are skipped, so an MGF file of one spectrum can be used; see section 2. In `topfd_gui`, the checkbox "Text peak list input" under "Additional settings" turns this on; the Add button then lists `.txt` and `.mgf` files (or all files), and Start refuses to run while an mzML or mzXML file is in the list. |

MS1 deconvolution and proteoform feature detection:

| Option | Default | Meaning |
|---|---|---|
| `-r`, `--ms-one-sn-ratio <number>` | 3 | Signal-to-noise ratio for MS1 spectra. In MS1 deconvolution it sets the reference-peak threshold (as `-s` does for MS/MS); in feature detection peaks below `ratio × noise level` are excluded from the LC-MS map. |
| `-p`, `--deconv-ms1 <greedy\|msdeconv>` | greedy | MS1 deconvolution method (section 1.2): `greedy`, the greedy search, or `msdeconv`, the MS-Deconv dynamic programming, the default method of earlier versions. In `topfd_gui`: "MS1 deconvolution". |
| `-f`, `--envcnn-cutoff-ms1 <0..1>` | 0.2 | MS1 envelopes whose EnvCNN score is below the cutoff are not reported, so they do not enter the MS1 msalign file, the SQLite database or the feature detection. Used with both MS1 deconvolution methods (`-p greedy` and `msdeconv`). Not filtered: the envelopes added for the precursor isolation windows of the MS/MS spectra (the greedy search adds one to a window in which no reported envelope has its apex; MS-Deconv deconvolutes the windows before the whole spectrum). The cutoff compares the EnvCNN score whatever `-n` ranks by. `-f 0` reports every envelope. In `topfd_gui`: "MS1 EnvCNN cutoff". |
| `-t`, `--ecscore-cutoff <0..1>` | 0.1 | Features with an ECScore below the cutoff are removed. |
| `-b`, `--min-scan-number <1\|2\|3>` | 1 | Minimum number of MS1 scans a feature must be detected in. |
| `-l`, `--split-intensity-ratio <number>` | 2.5 | Intensity ratio required to split one feature into two. |
| `-i`, `--single-scan-noise` | off | Use the noise level of each MS1 scan to filter low-intensity peaks, instead of the noise level of the whole LC-MS map. |
| `-d`, `--disable-additional-feature-search` | off | By default, an MS/MS spectrum whose isolation window contains no detected feature triggers an extra search of the LC-MS map with the S/N ratio set to 0, the minimum scan number to 1 and the ECScore cutoff to 0. This option disables the extra search. |

MS/MS deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|ETCID\|ETHCD\|HCD\|MPD\|UVPD\|FILE>` | FILE | Fragmentation method. `FILE` takes it from the activation terms of each spectrum in the input file, whatever their order: ETD with beam-type or supplemental beam-type CID (or the EThcD term) is ETHCD, ETD with CID or supplemental CID (or the ETciD term) is ETCID, otherwise ETD, HCD, CID or MPD. ETCID (ETciD) and ETHCD (EThcD) have the ion types of ETD (c, z-dot). Give a method explicitly when the file does not record it or records it wrongly. |
| `-s`, `--ms-two-sn-ratio <number>` | 1 | Signal-to-noise ratio for MS/MS spectra. Only peaks with intensity at least `ratio × noise level` can be the reference peak of an isotopic envelope; the other peaks of an envelope must be above the noise level. Values below 1, down to 0, lower both thresholds to `ratio × noise level`; the noise level itself still bounds how far isotopic envelopes extend. |
| `-w`, `--precursor-window <number>` | 3.0 | Default precursor isolation window width (m/z). Ignored when the file contains isolation window information. |
| `-q`, `--deconv-ms2 <greedy\|msdeconv\|flashdeconv>` | greedy | MS/MS deconvolution method (section 1.2): `greedy`, the greedy search adapted to fragment ions; `msdeconv`, the MS-Deconv dynamic programming, the default method of earlier versions; or `flashdeconv`, the FLASHDeconv algorithm (experimental; the output files are the same; not with `-T`). In `topfd_gui`: "MS2 deconvolution" (it replaces the former option `--ms2-flashdeconv`). |
| `-g`, `--envcnn-cutoff-ms2 <0..1>` | 0 | MS/MS envelopes whose EnvCNN score is below the cutoff are not reported. Used with `-q greedy` and `msdeconv`, not with `flashdeconv`; also for MS3 spectra and `-T`. In `topfd_gui`: "MS2 EnvCNN cutoff". |
| `-n`, `--env-rank <envcnn\|msdeconv>` | envcnn | The score that ranks isotopic envelopes: `envcnn`, the EnvCNN neural-network score, or `msdeconv`, the MS-Deconv score. The ranking decides which envelopes `-j` keeps; the EnvCNN cutoffs `-f` and `-g` compare the EnvCNN score with either value, and the score column of the msalign files is the EnvCNN score. In `topfd_gui`: "Sort using MS-Deconv score" (`msdeconv` when checked). |
| `-j`, `--frag-num-filtering` | off | Limit the number of fragment envelopes in an MS/MS spectrum based on the estimated number of fragment ions: the envelopes ranked highest by the score of `-n` are kept. Used with `-q greedy` and `msdeconv`, not with `flashdeconv`. In `topfd_gui`: "Enable frag num-based filtering". |

### 1.5 Examples

```sh
# 16 threads, no SQLite database
topfd -u 16 -N sample.mzML

# Low-mass proteins: limit charge and mass to speed up deconvolution
topfd -c 20 -m 30000 sample.mzML

# Stricter features: at least 2 MS1 scans and ECScore >= 0.5
topfd -b 2 -t 0.5 sample.mzML

# MS/MS-only file (no MS1 scans)
topfd -o sample.mzML
```

After TopFD finishes, `sample_ms2.msalign` (with `sample_ms1.msalign` and the
feature files in the same directory) is passed to TopPIC or TopMG.

## 2. Deconvoluting a peak list file

### 2.1 Input

The input is a plain-text file with one **centroided** peak per line: the
m/z value and the intensity, separated by spaces or a tab; further columns
are ignored. Blank lines and lines that do not start with a number are
skipped (TopFD prints how many), and the peaks need not be sorted: TopFD
sorts them by increasing m/z before deconvolution. The file describes a
single MS/MS spectrum; TopFD treats it as an MS level 2 spectrum without
precursor information.

```text
500.2513 12034.5
500.5857 25311.0
500.9199 21877.2
501.2542 11504.8
826.3891 8032.1
```

Use a space or a tab, not a comma, between the two columns. Any extension
may be used for the file name.

An MGF file holding one spectrum can be used as it is: its `BEGIN IONS` and
`END IONS` lines, its parameter lines (`TITLE=`, `PEPMASS=`, `CHARGE=`, ...)
and comments do not start with a number and are skipped, and the charge
column some peak lines carry is ignored. The precursor information of the
MGF file is not used. An MGF file with several spectra is read as one
spectrum with the peaks of all of them (TopFD warns when the file has more
than one `BEGIN IONS` line), so split it into one file per spectrum first.

```text
BEGIN IONS
TITLE=scan=1817
PEPMASS=1000.5123
CHARGE=6+
500.2513 12034.5
500.5857 25311.0
826.3891 8032.1
END IONS
```

### 2.2 Usage

Add `-T` (`--text-peak-list`) to tell TopFD that the input is a peak list:

```sh
topfd -T [options] peaks.txt
```

Only the parameters of MS/MS deconvolution apply. The useful ones are
`-c`/`--max-charge`, `-m`/`--max-mass`, `-e`/`--mz-error`,
`-s`/`--ms-two-sn-ratio`, `-n`/`--env-rank`, `-g`/`--envcnn-cutoff-ms2`,
`-j`/`--frag-num-filtering` and `-q`/`--deconv-ms2` (`greedy` or `msdeconv`;
`flashdeconv` is refused with `-T`). The MS1 and feature detection options
(`-r`, `-p`, `-f`, `-t`, `-b`, `-l`, `-i`, `-d`), `-o` and `-w` have
no effect; `-u` only sets the number of threads ONNX Runtime uses for EnvCNN
scoring. `-a`/`--activation` only sets the activation recorded
in the SQLite database; when it is not given (or is `FILE`), HCD is
recorded.

```sh
# deconvolute one spectrum with a maximum charge of 10 and a 0.01 m/z
# tolerance, without the SQLite database
topfd -T -c 10 -e 0.01 -N peaks.txt
```

### 2.3 Output files

As in mzML mode, the output files are written next to the input file, with
the input file name minus its extension as the base name (`peaks` for
`peaks.txt`). Each peak list given on the command line gets its own set:

| File | Contents |
|---|---|
| `peaks_ms2.msalign` | The deconvoluted spectrum in the `msalign` format (spectrum ID 0, one line per monoisotopic mass: mass, intensity, charge, score). |
| `peaks_ms2.env` | The matched isotopic envelopes, peak by peak (see below). |
| `peaks.sqlite` | The spectrum and its envelopes in the SQLite database format. Not written with `--no-sql`. |

`peaks_ms2.env` is a tab-separated table with one line per **matched
peak**, i.e. per input peak that was assigned to an isotopic envelope:

| Column | Meaning |
|---|---|
| `PEAK_IDX` | 0-based index of the peak after sorting the input by m/z (the index in the file when the file is already sorted). |
| `ORIG_MZ`, `ORIG_INTE` | The m/z and intensity of that input peak. |
| `THEO_MONO_MZ`, `THEO_MONO_MASS` | Monoisotopic m/z and neutral mass of the envelope the peak belongs to. |
| `THEO_INTE_SUM` | Total intensity of the theoretical envelope. |
| `THEO_CHARGE` | Charge state of the envelope. |
| `ENV_CNN_SCORE` | EnvCNN score of the envelope. |
| `THEO_MZ`, `THEO_NEUTRAL_MASS`, `THEO_INTE` | The m/z, neutral mass and intensity of the theoretical peak matched to this input peak. |

Peaks that were not assigned to any envelope do not appear in the file. The
lines of one envelope share the same `THEO_MONO_MASS` and `THEO_CHARGE`, so
the table can be grouped by those columns to see which input peaks make up
each reported mass.

The theoretical intensities (`THEO_INTE`, and `THEO_INTE_SUM`, which is also
the intensity reported in the `msalign` file) are scaled to the envelope's
core peaks, the most abundant isotopic peaks around the reference peak. Where
isotopic envelopes overlap, the observed intensities in an envelope's tails
can therefore be higher than the theoretical ones; the surplus belongs to
the neighbouring envelope.

### 2.4 Example

```sh
topfd -T -c 20 -m 30000 data/spectrum_1.txt
```

prints

```text
TopFD 2.1.0
Total thread number: 16
Total memory: 31.07 GiB
Available memory: 21.36 GiB

Processing data/spectrum_1.txt started.
Processing data/spectrum_1.txt finished.
Timestamp: ...
TopFD single finished.
```

and leaves `spectrum_1_ms2.msalign`, `spectrum_1_ms2.env` and
`spectrum_1.sqlite` in `data/`. The monoisotopic masses are in
`spectrum_1_ms2.msalign`; open `spectrum_1_ms2.env` in a spreadsheet to see
which input peaks support each mass.
