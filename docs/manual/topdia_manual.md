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

The `--missing-level-one` option (`-o`; TopFD's `--missing-ms1`) is accepted
for compatibility with TopFD,
but pseudo spectra need the MS1 features: with it TopDIA still deconvolutes
the MS/MS scans and then stops with an error (`msalign file sample_ms1.msalign
does not exist.`) at the pseudo-spectrum step. Do not use it.

## 2. What TopDIA does

For each file TopDIA prints its parameters, the TopFD deconvolution
parameters followed by the TopDIA feature and pseudo-spectrum parameters,
and runs five steps:

1. **MS1 deconvolution.** Every MS1 scan is deconvoluted, by default with
   the MS-Deconv dynamic programming (`--ms1-deconv msdeconv`); `--ms1-deconv
   greedy` selects TopFD's default greedy search. With either method, MS1
   envelopes whose EnvCNN score is below `--ms1-envcnn-cutoff` (`-f`,
   default 0.2) are not reported, except the envelopes of the precursor
   isolation windows.
2. **MS1 feature detection.** The MS1 envelopes are linked across scans into
   proteoform features, scored with ECScore and filtered by
   `--ms1-ecscore-cutoff` and `--ms1-min-scan-number`.
3. **MS/MS deconvolution.** Every DIA MS/MS scan is deconvoluted on its own,
   without precursor information, into `sample_ms2_raw.msalign`, by
   default with the MS-Deconv dynamic programming (`--ms2-deconv msdeconv`;
   `greedy` selects TopFD's default greedy search and `flashdeconv` the
   FLASHDeconv algorithm). With `msdeconv` and `greedy`, envelopes whose
   EnvCNN score is below `--ms2-envcnn-cutoff` (`-g`, default 0: all
   reported) are not reported. For the dynamic programming TopDIA
   uses the parameters of version 1.9.0.0: a 12
   Dalton envelope ranking window, no harmonic check and a 1,500 Dalton mass
   group boundary. The 2.0.0.0 defaults of TopFD (2 Dalton, harmonic check,
   4,250 Dalton) were measured on DDA data; a DIA MS/MS scan holds the
   fragments of every co-isolated precursor of a wide window, where they
   roughly double the reported masses and the fragment features that step 5
   must hold in memory.
4. **MS/MS feature detection.** The deconvoluted MS/MS scans are grouped by
   isolation window, and within each window the fragment masses are linked
   across scans into fragment features, scored with ECScore and filtered by
   `--ms2-ecscore-cutoff` and `--ms2-min-scan-number`.
5. **Pseudo spectrum generation.** For each isolation window, the MS1
   features whose envelope lies mostly inside the window are taken in order
   of decreasing intensity, each charge state of a feature on its own. For
   each such precursor, the fragment features of the window that are lighter
   than the precursor and elute at the same time (apex within a few MS1
   cycles) are candidates. Each candidate gets a
   **pseudo score** in 0 to 1 from a logistic model of its intensity rank,
   the ratio of its elution length to the precursor's, and the overlap of the
   two elution profiles. Candidates are added to the precursor's pseudo
   spectrum in decreasing score order while the score is at least
   `--pseudo-cutoff` or fewer than `--pseudo-peak-number` fragments have been
   accepted, and a fragment feature assigned to one precursor is not reused
   for another. In addition, at most 25 fragments with mass up to 1,500 Da
   and at most 2 × (precursor mass − 1,500) / 120 fragments above it are
   accepted.

## 3. Output files

For an input `sample.mzML`:

| File | Contents |
|---|---|
| `sample_ms2.msalign` | **The pseudo-MS/MS spectra, the file to search with TopPIC or TopMG.** One spectrum per charge state of a precursor feature and isolation window, with the precursor mass, charge, intensity and feature id in its header, and one line per fragment mass with its intensity, charge and pseudo score (followed by the elution-profile statistics behind the score). |
| `sample_ms2.feature` | TopFD's feature file for the DIA MS/MS scans: for each scan, the most intense MS1 feature in its isolation window. TopPIC and TopMG use it to attach feature intensities and elution times to the pseudo spectra through the `PRECURSOR_FEATURE_ID` in the msalign header; no special option is needed. |
| `sample_ms1.msalign` | The deconvoluted MS1 spectra. |
| `sample_ms1.feature`, `sample_feature.xml` | The proteoform features detected in the MS1 data, as text and as XML (the XML is read by TopMG when it combines fractions with `-c`). |
| `sample_ms2_raw.msalign` | The deconvoluted DIA MS/MS scans before pseudo-spectrum generation; intermediate. |
| `sample_ms1.csv`, `sample_frac_ms1.mzrt.csv` | The MS1 feature ECScore table and the MS1 features with their elution profiles (BatMass format); intermediate. |
| `sample_<window>_ms2.csv`, `sample_<window>_frac_ms2.mzrt.csv` | The same two tables for the fragment features of each isolation window, named by the window's lower m/z bound written with six decimals (`sample_400.000000_ms2.csv`); intermediate. |
| `sample.sqlite` | The SQLite database with the deconvoluted MS1 and MS/MS scans and their peaks, as written by TopFD, with the levels of the 3D visualisation of the MS1 peaks (the `LEVEL` column of `ms1_peak` and the table `CONFIG`, see the [TopFD manual](topfd_manual.md)). Not written with `-z` (`--no-sql`). |

The intermediate files are left in place. When searching the pseudo spectra
with TopPIC, add `-D`: post mass matching looks up the
centroided peaks of a single MS/MS scan, and a pseudo spectrum does not
correspond to one scan.

## 4. Options

TopDIA reuses TopFD's deconvolution and MS1 feature detection, so most
options are TopFD's (see the [TopFD manual](topfd_manual.md)). The defaults
that differ from TopFD are marked.

Deconvolution:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|ETCID\|ETHCD\|HCD\|MPD\|UVPD\|FILE>` | FILE | Fragmentation method; `FILE` takes it from the file. ETCID (ETciD) and ETHCD (EThcD) have the ion types of ETD (c, z-dot). |
| `-c`, `--max-charge <int>` | 30 | Maximum charge state of precursor and fragment ions. |
| `-m`, `--max-mass <number>` | 50000 | Maximum monoisotopic mass (Da) of precursor and fragment ions. |
| `-e`, `--mz-error <number>` | 0.02 | Error tolerance of peak m/z values. |
| `-r`, `--ms1-sn-ratio <number>` | 3 | Signal-to-noise ratio for MS1 spectra. |
| `-s`, `--ms2-sn-ratio <number>` | 1 | Signal-to-noise ratio for MS/MS spectra. |
| `-w`, `--precursor-window <number>` | 4.0 (TopFD: 3.0) | Default isolation window width (m/z), used only when the file records none. |
| `-p`, `--ms1-deconv <greedy\|msdeconv>` | msdeconv (TopFD: greedy) | MS1 deconvolution method (section 2, step 1). In `topdia_gui`: "MS1 deconvolution". |
| `-q`, `--ms2-deconv <greedy\|msdeconv\|flashdeconv>` | msdeconv (TopFD: greedy) | MS/MS deconvolution method (section 2, step 3); `flashdeconv` is experimental. In `topdia_gui`: "MS2 deconvolution". |
| `-f`, `--ms1-envcnn-cutoff <0..1>` | 0.2 | MS1 envelopes whose EnvCNN score is below the cutoff are not reported (section 2, step 1); not the envelopes of the precursor isolation windows. Used with both MS1 deconvolution methods. `-f 0` reports every envelope. |
| `-g`, `--ms2-envcnn-cutoff <0..1>` | 0 | MS/MS envelopes whose EnvCNN score is below the cutoff are not reported (section 2, step 3). Used with `-q msdeconv` and `greedy`, not with `flashdeconv`. |
| `-n`, `--env-rank <envcnn\|msdeconv>` | envcnn | The score that ranks isotopic envelopes: `envcnn`, the EnvCNN score, or `msdeconv`, the MS-Deconv score. The ranking decides which envelopes `-j` keeps; the EnvCNN cutoffs `-f` and `-g` compare the EnvCNN score with either value. In `topdia_gui`: "Use MS-Deconv score" (`msdeconv` when checked). |
| `-j`, `--frag-num-filtering` | off | Filter the envelopes of MS/MS scans by the estimated number of fragment ions: the envelopes ranked highest by the score of `-n` are kept. Not used with `-q flashdeconv`. In `topdia_gui`: "Enable frag num filtering". |
| `-o`, `--missing-level-one` | off | Accepted, but the run fails at the pseudo-spectrum step without MS1 scans (section 1). |
| `-u`, `--thread-number <int>` | 1 | Number of threads. |

Feature detection:

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--ms1-ecscore-cutoff <0..1>` | 0 (TopFD: 0.1) | ECScore cutoff for MS1 proteoform features. |
| `-b`, `--ms1-min-scan-number <1\|2\|3>` | 2 (TopFD: 1) | Minimum number of MS1 scans of a proteoform feature. |
| `-T`, `--ms2-ecscore-cutoff <0..1>` | 0 | ECScore cutoff for fragment features. |
| `-B`, `--ms2-min-scan-number <1\|2\|3>` | 1 | Minimum number of MS/MS scans of a fragment feature. |
| `-i`, `--single-scan-noise` | off | Use each scan's own noise level instead of the whole map's when filtering low-intensity peaks: the MS1 scans for proteoform features and the MS/MS scans of an isolation window for fragment features. |
| `-y`, `--ms1-correlation-cutoff <0..1>` | 0.5 | Minimum correlation between the intensities of a seed isotopic envelope and the theoretical ones for starting an MS1 proteoform feature. |
| `-Y`, `--ms2-correlation-cutoff <0..1>` | 0.5 | The same cutoff for starting a fragment feature in an isolation window. |

Pseudo spectra:

| Option | Default | Meaning |
|---|---|---|
| `-v`, `--pseudo-cutoff <0..1>` | 0.55 | Minimum pseudo score of a fragment added to a pseudo spectrum, once the minimum number of fragments is reached. |
| `-V`, `--pseudo-peak-number <int>` | 25 | Number of fragments (at least 10) always added to a pseudo spectrum, regardless of their score, when available. |

Output:

| Option | Default | Meaning |
|---|---|---|
| `-z`, `--no-sql` | off | Do not write the `.sqlite` database (section 3); an earlier `sample.sqlite` in the output folder is deleted with its `-wal` and `-shm` files. `toppic -D`, the setting for pseudo spectra, does not read it, and then writes no `sample_toppic.sqlite` either. Without `-D`, TopPIC's post mass matching stops with an error when the database is missing. |

TopDIA has no text-peak-list mode; it is TopFD only.

## 5. Example

```sh
topdia -u 8 sample.mzML
toppic -D -d -t FDR -v 0.01 proteins.fasta sample_ms2.msalign
```

The first command writes `sample_ms2.msalign` (pseudo spectra) and
`sample_ms2.feature` next to `sample.mzML`; the second searches them.
