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

# several files, 8 threads
topfd -u 8 sample_1.mzML sample_2.mzML
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
| `sample_ms2.msalign` | Deconvoluted MS/MS spectra with precursor information. **This is the input of TopPIC and TopMG.** |
| `sample_ms1.feature` | Proteoform features detected in the LC-MS map: one line per feature with its mass, intensity, retention-time and scan range, charge range, apex and ECScore. |
| `sample_ms2.feature` | The features assigned to each MS/MS spectrum (precursor mass, m/z, charge and intensity), one line per spectrum and feature. |
| `sample_feature.xml` | The proteoform features in XML, with their per-charge envelope information; read by TopPIC and TopMG when they combine fractions (`-c`). |
| `sample.sqlite` | An SQLite database with the deconvoluted MS1 and MS/MS spectra and their peaks, for spectrum visualisation and for TopPIC's post mass matching (see the `-N` option). Not written with `--no-sql`. With `--sql-3d` it also holds the raw MS1 peak tables for 3D visualisation (`CONFIG` and `PEAKS0`, `PEAKS1`, ...). |

The `msalign` format is a text format. Each spectrum is a `BEGIN IONS` ...
`END IONS` block with header lines (`FILE_NAME`, `SPECTRUM_ID`, `TITLE`,
`SCANS`, `RETENTION_TIME` in minutes, `LEVEL`, and for MS/MS spectra
`MS_ONE_ID`, `MS_ONE_SCAN`, `PRECURSOR_WINDOW_BEGIN/END`, `ACTIVATION`,
`PRECURSOR_MZ`, `PRECURSOR_CHARGE`, `PRECURSOR_MASS`, `PRECURSOR_INTENSITY`,
`PRECURSOR_FEATURE_ID`, `DECONVOLUTED_MASS_NUMBER`) followed by one line per
deconvoluted mass with four tab-separated columns: monoisotopic neutral
mass, intensity, charge and score. When several features fall in the
isolation window of an MS/MS spectrum, those with at least a tenth of the
strongest feature's intensity are all kept and the `PRECURSOR_*` values are
separated by `:`. The parameters used are recorded as `#` comment lines at
the top of the file. With `--missing-level-one` no `_ms1.msalign` and no
feature files are produced.

With more than one thread, per-thread partial `msalign` files
(`sample_ms1.msalign_0`, `sample_ms2.msalign_0`, `_1`, ...) exist while
TopFD runs; they are merged into the final files and deleted at the end.

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
| `-D`, `--sql-3d` | off | Also store the raw MS1 peaks for 3D visualisation in the `.sqlite` database: `PEAKS0` holds every MS1 peak and `PEAKS1`, `PEAKS2`, ... progressively down-sampled copies, with one `CONFIG` row per table. Cannot be combined with `--no-sql`; has no effect with `-T` or `-o`. In `topfd_gui`, the checkbox "Add MS1 peaks for 3D visualization" under "Additional settings" turns this on (it is greyed out while "Do not generate SQLite database" is checked). |
| `-T`, `--text-peak-list` | off | The input is a text peak list (one MS/MS spectrum), not an mzML file; see section 2. |

MS1 deconvolution and proteoform feature detection:

| Option | Default | Meaning |
|---|---|---|
| `-r`, `--ms-one-sn-ratio <number>` | 3 | Signal-to-noise ratio for MS1 spectra. In MS1 deconvolution it sets the reference-peak threshold (as `-s` does for MS/MS); in feature detection peaks below `ratio × noise level` are excluded from the LC-MS map. |
| `-t`, `--ecscore-cutoff <0..1>` | 0.1 | Features with an ECScore below the cutoff are removed. |
| `-b`, `--min-scan-number <1\|2\|3>` | 1 | Minimum number of MS1 scans a feature must be detected in. |
| `-l`, `--split-intensity-ratio <number>` | 2.5 | Intensity ratio required to split one feature into two. |
| `-i`, `--single-scan-noise` | off | Use the noise level of each MS1 scan to filter low-intensity peaks, instead of the noise level of the whole LC-MS map. |
| `-f`, `--disable-additional-feature-search` | off | By default, an MS/MS spectrum whose isolation window contains no detected feature triggers an extra search of the LC-MS map with the S/N ratio set to 0, the minimum scan number to 1 and the ECScore cutoff to 0. This option disables the extra search. |

MS/MS deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|HCD\|MPD\|UVPD\|FILE>` | FILE | Fragmentation method. `FILE` takes it from each spectrum in the input file; give a method explicitly when the file does not record it or records it wrongly. |
| `-s`, `--ms-two-sn-ratio <number>` | 1 | Signal-to-noise ratio for MS/MS spectra. Only peaks with intensity at least `ratio × noise level` can be the reference peak of an isotopic envelope; the other peaks of an envelope must be above the noise level. Values below 1, down to 0, lower both thresholds to `ratio × noise level`; the noise level itself still bounds how far isotopic envelopes extend. |
| `-w`, `--precursor-window <number>` | 3.0 | Default precursor isolation window width (m/z). Ignored when the file contains isolation window information. |
| `-n`, `--msdeconv` | off | Rank isotopic envelopes with the MS-Deconv score instead of the EnvCNN neural-network score. |
| `-v`, `--env-cnn-cutoff <0..1>` | 0 | Remove MS/MS envelopes whose EnvCNN score is below the cutoff. |
| `-g`, `--frag-num-filtering` | off | Limit the number of fragment envelopes in an MS/MS spectrum based on the estimated number of fragment ions. |

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
m/z value and the intensity, separated by a space. Blank lines are ignored,
and the peaks need not be sorted: TopFD sorts them by increasing m/z before
deconvolution. The file describes a single MS/MS spectrum; TopFD treats it as
an MS level 2 spectrum without precursor information.

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
`-s`/`--ms-two-sn-ratio`, `-n`/`--msdeconv`, `-v`/`--env-cnn-cutoff` and
`-g`/`--frag-num-filtering`. The MS1 and feature
detection options (`-r`, `-t`, `-b`, `-l`, `-i`, `-f`), `-o` and `-w` have
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
TopFD 1.9.0
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
