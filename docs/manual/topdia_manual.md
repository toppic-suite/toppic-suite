# TopDIA manual

TopDIA analyzes top-down data-independent-acquisition (TD-DIA) mass
spectrometry data. In DIA, each MS/MS scan fragments every precursor in a
wide isolation window together, so the fragment spectra cannot be assigned
to one precursor directly. TopDIA deconvolutes the MS1 and MS/MS scans,
detects proteoform features in the MS1 data and fragment features in each
isolation window, and matches fragments to precursors by their elution
profiles to generate **pseudo-MS/MS spectra**, one per precursor feature.
The pseudo spectra are then searched with TopPIC or TopMG like ordinary
deconvoluted spectra.

```sh
topdia [options] spectrum-file ...
```

`topdia -h` prints the option list. The GUI tool `topdia_gui` offers the
same function. Paths must not contain spaces.

## 1. Input

TopDIA reads **centroided** mzML or mzXML files (`.mzML`, `.mzml`,
`.mzXML`, `.mzxml`), converted from the vendor format with peak picking,
for example:

```sh
msconvert --mzML --filter "peakPicking true 1-" sample.raw
```

The file must contain both MS1 scans and DIA MS/MS scans whose isolation
windows are recorded in the file: the set of distinct windows drives the
whole analysis. Several files on one command line are processed one after
another. FAIMS data is recognized automatically and each compensation voltage
is processed as a separate fraction, with the integer voltage appended to the
output file names (`sample_-40_ms2.msalign`).

The `--missing-level-one` option is accepted for compatibility with TopFD,
but pseudo spectra need the MS1 features, so TopDIA cannot produce them
without MS1 scans.

## 2. What TopDIA does

For each file TopDIA prints its parameters, the TopFD deconvolution
parameters followed by the TopDIA feature and pseudo-spectrum parameters,
and runs five steps:

1. **MS1 deconvolution.** Every MS1 scan is deconvoluted, as in TopFD.
2. **MS1 feature detection.** The MS1 envelopes are linked across scans into
   proteoform features, scored with ECScore and filtered by
   `--ms1-ecscore-cutoff` and `--ms1-min-scan-number`.
3. **MS/MS deconvolution.** Every DIA MS/MS scan is deconvoluted on its own,
   without precursor information, into `sample_ms2_raw.msalign`.
4. **MS/MS feature detection.** The deconvoluted MS/MS scans are grouped by
   isolation window, and within each window the fragment masses are linked
   across scans into fragment features, scored with ECScore and filtered by
   `--ms2-ecscore-cutoff` and `--ms2-min-scan-number`.
5. **Pseudo spectrum generation.** For each isolation window, the MS1
   features whose envelope lies mostly inside the window are taken in order
   of decreasing intensity. For each such precursor, the fragment features of
   the window that are lighter than the precursor and elute at the same time
   (apex within a few MS1 cycles) are candidates. Each candidate gets a
   **pseudo score** in 0 to 1 from a logistic model of its intensity rank,
   the ratio of its elution length to the precursor's, and the overlap of the
   two elution profiles. Candidates are added to the precursor's pseudo
   spectrum in decreasing score order while the score is at least
   `--pseudo-cutoff` or fewer than `--pseudo-peak-number` fragments have been
   accepted, and a fragment feature assigned to one precursor is not reused
   for another.

## 3. Output files

For an input `sample.mzML`:

| File | Contents |
|---|---|
| `sample_ms2.msalign` | **The pseudo-MS/MS spectra, the file to search with TopPIC or TopMG.** One spectrum per precursor feature and isolation window, with the precursor mass, charge, intensity and feature id in its header, and one line per fragment mass with its intensity, charge and pseudo score (followed by the elution-profile statistics behind the score). |
| `sample_ms2.feature` | The precursor feature of each pseudo spectrum, in the format TopPIC and TopMG read; no special option is needed to search the pseudo spectra with them. |
| `sample_ms1.msalign` | The deconvoluted MS1 spectra. |
| `sample_ms1.feature`, `sample_feature.xml` | The proteoform features detected in the MS1 data, as text and as XML (the XML is used by TopDiff and by the visualization tools). |
| `sample_ms2_raw.msalign` | The deconvoluted DIA MS/MS scans before pseudo-spectrum generation; intermediate. |
| `sample_ms1.csv`, `sample_frac_ms1.mzrt.csv` | The MS1 feature ECScore table and the MS1 features with their elution profiles (BatMass format); intermediate. |
| `sample_<window>_ms2.csv`, `sample_<window>_frac_ms2.mzrt.csv` | The same two tables for the fragment features of each isolation window, named by the window's lower m/z bound; intermediate. |
| `sample.sqlite` | The SQLite database with the deconvoluted MS1 and MS/MS scans and their peaks, as written by TopFD. It is always written. |

The intermediate files are left in place. When searching the pseudo spectra
with TopPIC, add `--disable-post-match`: post mass matching looks up the
centroided peaks of a single MS/MS scan, and a pseudo spectrum does not
correspond to one scan.

## 4. Options

TopDIA reuses TopFD's deconvolution and MS1 feature detection, so most
options are TopFD's (see the [TopFD manual](topfd_manual.md)). The defaults
that differ from TopFD are marked.

Deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|HCD\|MPD\|UVPD\|FILE>` | FILE | Fragmentation method; `FILE` takes it from the file. |
| `-c`, `--max-charge <int>` | 30 | Maximum charge state of precursor and fragment ions. |
| `-m`, `--max-mass <number>` | 50000 | Maximum monoisotopic mass (Da). The help text says 70,000, but the value in effect is TopFD's 50,000. |
| `-e`, `--mz-error <number>` | 0.02 | Error tolerance of peak m/z values. |
| `-r`, `--ms-one-sn-ratio <number>` | 3 | Signal-to-noise ratio for MS1 spectra. |
| `-s`, `--ms-two-sn-ratio <number>` | 1 | Signal-to-noise ratio for MS/MS spectra. |
| `-w`, `--precursor-window <number>` | 4.0 (TopFD: 3.0) | Default isolation window width (m/z), used only when the file records none. |
| `-n`, `--msdeconv` | off | Rank isotopic envelopes with the MS-Deconv score instead of EnvCNN. |
| `-d`, `--final-filtering` | off | Filter the envelopes of MS/MS scans by the estimated number of fragment ions. |
| `-o`, `--missing-level-one` | off | Accepted, but pseudo spectra cannot be generated without MS1 scans (section 1). |
| `-u`, `--thread-number <int>` | 1 | Number of threads. |

Feature detection:

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--ms1-ecscore-cutoff <0..1>` | 0 (TopFD: 0.1) | ECScore cutoff for MS1 proteoform features. |
| `-b`, `--ms1-min-scan-number <1\|2\|3>` | 2 (TopFD: 1) | Minimum number of MS1 scans of a proteoform feature. |
| `-T`, `--ms2-ecscore-cutoff <0..1>` | 0 | ECScore cutoff for fragment features. |
| `-B`, `--ms2-min-scan-number <1\|2\|3>` | 1 | Minimum number of MS/MS scans of a fragment feature. |
| `-i`, `--single-scan-noise` | off | Use each MS1 scan's own noise level instead of the whole map's when filtering low-intensity peaks. |
| `-p`, `--ms1-intensity-correlation-cutoff <0..1>`, `-P`, `--ms2-intensity-correlation-cutoff <0..1>` | 0.5, 0 | Accepted and printed, but not used by the command-line pipeline in this version; feature detection applies a fixed correlation cutoff of 0.5. |

Pseudo spectra:

| Option | Default | Meaning |
|---|---|---|
| `-v`, `--pseudo-cutoff <0..1>` | 0.55 | Minimum pseudo score of a fragment added to a pseudo spectrum, once the minimum number of fragments is reached. |
| `-V`, `--pseudo-peak-number <int>` | 25 | Number of fragments (at least 10) always added to a pseudo spectrum, regardless of their score, when available. |

Advanced option (accepted but not shown by `-h`): `-k`, `--keep`, report
masses from low-quality envelopes as well. TopDIA has no `--no-sql`,
`--sql-3d` or text-peak-list mode; those are TopFD only.

## 5. Example

```sh
topdia -u 8 sample.mzML
toppic --disable-post-match -d -t FDR -v 0.01 proteins.fasta sample_ms2.msalign
```

The first command writes `sample_ms2.msalign` (pseudo spectra) and
`sample_ms2.feature` next to `sample.mzML`; the second searches them.
