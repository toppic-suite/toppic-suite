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

`topfd -h` prints the option list. The GUI tool `topfd_gui` offers the first
function only.

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

# several files, 8 threads, HCD fragmentation for all MS/MS spectra
topfd -u 8 -a HCD sample_1.mzML sample_2.mzML
```

For each file TopFD prints the parameters it uses, then runs three steps:

1. **MS1 deconvolution**: every MS1 spectrum is deconvoluted.
2. **MS1 feature detection**: the deconvoluted MS1 envelopes are linked
   across scans into proteoform features (LC-MS features), each scored with
   the ECScore neural network. Every MS/MS spectrum is then matched to the
   feature(s) in its precursor isolation window, which gives its precursor
   monoisotopic mass, charge and intensity.
3. **MS/MS deconvolution**: every MS/MS spectrum is deconvoluted, using its
   precursor information from step 2.

If the file has no MS1 spectra, pass `--missing-level-one` (`-o`): steps 1 and
2 are skipped, and the precursor of each MS/MS spectrum is taken from the
file (or bounded by the maximum charge and mass when the file has none).

### 1.3 Output files

All output files are written next to the input file, with the input file
name minus its extension as the base name (`sample` for `sample.mzML`; for
FAIMS data `sample_<voltage>` per voltage level, e.g. `sample_-40`).

| File | Contents |
|---|---|
| `sample_ms1.msalign` | Deconvoluted MS1 spectra (monoisotopic masses). |
| `sample_ms2.msalign` | Deconvoluted MS/MS spectra with precursor information. **This is the input of TopPIC and TopMG.** |
| `sample_ms1.feature` | Proteoform features detected in the LC-MS map: one line per feature with its mass, intensity, retention-time and scan range, charge range, apex and ECScore. |
| `sample_ms2.feature` | The feature assigned to each MS/MS spectrum (precursor mass, m/z, charge and intensity). |
| `sample_feature.xml` | The proteoform features in XML, with their per-charge envelope information; used by TopDiff. |
| `sample.sqlite` | An SQLite database with the deconvoluted MS1 and MS/MS spectra and their peaks, for spectrum visualisation. Not written with `--no-sql`. |
| `sample_ms1.csv`, `sample_frac_ms1.mzrt.csv` | Only with `--output-batmass-feature`: the ECScore table and the features in the BatMass CSV format. |

The `msalign` format is a text format. Each spectrum is a `BEGIN IONS` ...
`END IONS` block with header lines (`SPECTRUM_ID`, `SCANS`,
`RETENTION_TIME`, `LEVEL`, and for MS/MS spectra `MS_ONE_ID`,
`MS_ONE_SCAN`, `PRECURSOR_WINDOW_BEGIN/END`, `ACTIVATION`, `PRECURSOR_MZ`,
`PRECURSOR_CHARGE`, `PRECURSOR_MASS`, `PRECURSOR_INTENSITY`,
`PRECURSOR_FEATURE_ID`) followed by one line per deconvoluted mass with four
tab-separated columns: monoisotopic neutral mass, intensity, charge and
score. The parameters used are recorded as `#` comment lines at the top of
the file. With `--missing-level-one` no `_ms1.msalign` and no feature files
are produced.

With more than one thread, per-thread partial `msalign` files
(`sample_ms2.msalign_0`, `_1`, ...) exist while TopFD runs; they are merged
into `sample_ms2.msalign` and deleted at the end.

### 1.4 Options

Parameters that apply to both MS1 and MS/MS deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-c`, `--max-charge <int>` | 30 | Maximum charge state of precursor and fragment ions. |
| `-m`, `--max-mass <number>` | 50000 | Maximum monoisotopic mass (Da) of precursor and fragment ions. |
| `-e`, `--mz-error <number>` | 0.02 | Error tolerance of peak m/z values (m/z units). |
| `-u`, `--thread-number <int>` | 1 | Number of threads. TopFD checks that the machine has enough memory for the requested number. |
| `-o`, `--missing-level-one` | off | The file has no MS1 spectra: skip MS1 deconvolution and feature detection. |
| `-N`, `--no-sql` | off | Do not write the `.sqlite` database. |

MS1 deconvolution and proteoform feature detection:

| Option | Default | Meaning |
|---|---|---|
| `-r`, `--ms-one-sn-ratio <number>` | 3 | Signal-to-noise ratio for MS1 spectra; peaks below it are discarded. |
| `-t`, `--ecscore-cutoff <0..1>` | 0.1 | Features with an ECScore below the cutoff are removed. |
| `-b`, `--min-scan-number <1\|2\|3>` | 1 | Minimum number of MS1 scans a feature must be detected in. |
| `-l`, `--split-intensity-ratio <number>` | 2.5 | Intensity ratio required to split one feature into two. |
| `-i`, `--single-scan-noise` | off | Use the noise level of each MS1 scan to filter low-intensity peaks, instead of the noise level of the whole LC-MS map. |
| `-f`, `--disable-additional-feature-search` | off | By default, an MS/MS spectrum whose isolation window contains no detected feature triggers an extra search of the LC-MS map with the S/N ratio set to 0, the minimum scan number to 1 and the ECScore cutoff to 0. This option disables the extra search. |

MS/MS deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|HCD\|MPD\|UVPD\|FILE>` | FILE | Fragmentation method. `FILE` takes it from each spectrum in the input file; give a method explicitly when the file does not record it or records it wrongly. |
| `-s`, `--ms-two-sn-ratio <number>` | 1 | Signal-to-noise ratio for MS/MS spectra. |
| `-w`, `--precursor-window <number>` | 3.0 | Default precursor isolation window width (m/z). Ignored when the file contains isolation window information. |
| `-n`, `--msdeconv` | off | Rank isotopic envelopes with the MS-Deconv score instead of the EnvCNN neural-network score. |
| `-v`, `--env-cnn-cutoff <0..1>` | 0 | Remove MS/MS envelopes whose EnvCNN score is below the cutoff. |
| `-d`, `--disable-frag-num-filtering` | off | Skip the filter that limits the number of fragment envelopes in an MS/MS spectrum based on the estimated number of fragment ions. |

Advanced options (accepted but not shown by `-h`):

| Option | Meaning |
|---|---|
| `-T`, `--text-peak-list` | The input is a text peak list, not an mzML file (see section 2). |
| `-k`, `--keep` | Also report monoisotopic masses from low-quality envelopes. |
| `-M`, `--multiple-mass` | Output several candidate monoisotopic masses per envelope for MS/MS spectra. |
| `-O`, `--output-batmass-feature` | Also write the feature files in the BatMass CSV format. |
| `--max-miss-peak-num <int>` | Maximum number of missing peaks allowed in a matched envelope (default 1). |
| `--disable-filter-by-mz` | Skip the step that removes an envelope outranked by a higher-scoring neighbour with the same charge. |
| `--output-dp-envs` | Dump the candidate envelopes of every spectrum to `win_envs.txt` / `dp_envs.txt` (debugging). |

### 1.5 Examples

```sh
# Thermo HCD data, 16 threads, no SQLite database
topfd -a HCD -u 16 -N sample.mzML

# Low-mass proteins: limit charge and mass to speed up deconvolution
topfd -c 20 -m 30000 sample.mzML

# Stricter features: at least 2 MS1 scans and ECScore >= 0.5
topfd -b 2 -t 0.5 sample.mzML

# MS/MS-only file (no MS1 scans), ETD fragmentation
topfd -o -a ETD sample.mzML
```

After TopFD finishes, `sample_ms2.msalign` (with `sample_ms1.msalign` and the
feature files in the same directory) is passed to TopPIC or TopMG.

## 2. Deconvoluting a peak list file

### 2.1 Input

The input is a plain-text file with one **centroided** peak per line: the
m/z value and the intensity, separated by a space. Blank lines are ignored.
The file describes a single MS/MS spectrum; TopFD treats it as an MS level 2
spectrum without precursor information.

```text
500.2513 12034.5
500.5857 25311.0
500.9199 21877.2
501.2542 11504.8
826.3891 8032.1
```

Use a space, not a tab or a comma, between the two columns. Any extension
may be used for the file name.

### 2.2 Usage

Add `-T` (`--text-peak-list`) to tell TopFD that the input is a peak list:

```sh
topfd -T [options] peaks.txt
```

Only the parameters of MS/MS deconvolution apply. The useful ones are
`-c`/`--max-charge`, `-m`/`--max-mass`, `-e`/`--mz-error`,
`-s`/`--ms-two-sn-ratio`, `-n`/`--msdeconv`, `-v`/`--env-cnn-cutoff`,
`-d`/`--disable-frag-num-filtering`, `-k`/`--keep`, `-M`/`--multiple-mass`,
`--max-miss-peak-num` and `--disable-filter-by-mz`. The MS1 and feature
detection options (`-r`, `-t`, `-b`, `-l`, `-i`, `-f`), `-o`, `-w` and
`-u` have no effect. `-a`/`--activation` only sets the activation recorded
in the SQLite database; when it is not given (or is `FILE`), HCD is
recorded.

```sh
# deconvolute one spectrum with a maximum charge of 10 and a 0.01 m/z
# tolerance, without the SQLite database
topfd -T -c 10 -e 0.01 -N peaks.txt
```

### 2.3 Output files

The output files are written to the **current working directory** with the
fixed base name `deconv`, whatever the input file is called:

| File | Contents |
|---|---|
| `deconv_ms2.msalign` | The deconvoluted spectrum in the `msalign` format (spectrum ID 0, one line per monoisotopic mass: mass, intensity, charge, score). |
| `deconv_ms2.env` | The matched isotopic envelopes, peak by peak (see below). |
| `deconv.sqlite` | The spectrum and its envelopes in the SQLite database format. Not written with `--no-sql`. |

Because the base name is fixed, running TopFD on a second peak list in the
same directory overwrites the previous output. Run each peak list from its
own directory, or rename the output files between runs. Several peak list
files given on one command line are processed in order and all write the
same `deconv_*` files, so only the last one survives.

`deconv_ms2.env` is a tab-separated table with one line per **matched
peak**, i.e. per input peak that was assigned to an isotopic envelope:

| Column | Meaning |
|---|---|
| `PEAK_IDX` | 0-based index of the peak in the input file. |
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

### 2.4 Example

```sh
mkdir spectrum_1 && cd spectrum_1
topfd -T -c 20 -m 30000 ../spectrum_1.txt
```

prints

```text
TopFD 1.9.0
Processing ../spectrum_1.txt started.
Processing ../spectrum_1.txt finished.
Timestamp: ...
TopFD single finished.
```

and leaves `deconv_ms2.msalign`, `deconv_ms2.env` and `deconv.sqlite` in
`spectrum_1/`. The monoisotopic masses are in `deconv_ms2.msalign`; open
`deconv_ms2.env` in a spreadsheet to see which input peaks support each mass.
