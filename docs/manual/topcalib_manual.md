# topcalib manual

topcalib recalibrates the masses of a TopFD MS/MS msalign file in two
passes. A first TopPIC search without unexpected shifts gives confident
PrSMs; their precursor and fragment mass errors are fitted as a function
of the m/z and the retention time (`topcalib fit`); the model then
corrects every precursor and fragment mass of the msalign file
(`topcalib apply`), which is searched again, typically with a narrower
fragment tolerance (`toppic -e 5`). The correction works on TopFD's masses,
so TopFD is not run again. It is useful when the instrument's mass
calibration is off by a few ppm or drifts during the run, as on TOF
instruments; on a well-calibrated Orbitrap the correction is small and the
gain comes from the narrower tolerance it allows.

```sh
topcalib fit [options] <PrSM table> <MS/MS msalign file> <model file>
topcalib apply <model file> <MS/MS msalign file> <output msalign file>
```

`topcalib fit -h` and `topcalib apply -h` print the options.

## 1. Input

- **`fit`**:
  1. **A PrSM table** of a TopPIC search of the msalign file, the file
     `<base>_ms2_toppic_prsm_single.tsv` (`<base>_post_ms2_...` without
     `-D`): the PrSMs that pass the spectrum-level cutoff, one row each.
     The calibrants are selected by E-value (`-e`), so the cutoff only
     needs to keep them (`-t FDR` at 1% FDR does). A search without
     unexpected shifts (`toppic -s 0`) is the intended input: it is quick,
     and its confident PrSMs are the calibrants. The post mass matching is
     not needed (`-D`). The table must have the
     columns `Scan(s)`, `Retention time`, `Charge`, `Precursor mass`,
     `Proteoform mass`, `Proteoform`, `Fragmentation`,
     `# Unexpected modifications` and `E-value`; `# Variable PTMs` is used
     when present.
  2. **The MS/MS msalign file** that the search read
     (`<base>_ms2.msalign`), for the fragment masses of the PrSMs'
     spectra (by scan).
- **`apply`**: a model file written by `topcalib fit` and an MS/MS msalign
  file of TopFD (the one the model was fitted on, or another TopFD output
  of the same run, e.g. with other deconvolution options).

## 2. What topcalib does

### `fit`

1. **Reads the calibrants** (`Reading the calibrants`). The PrSMs with an
   E-value of at most `-e` (1e-5), no unexpected shift and no variable PTM
   are used.
   - **Precursors**: the error of the precursor mass against the
     proteoform mass, (precursor mass - proteoform mass) / proteoform mass
     in ppm, for the PrSMs whose difference rounds to 0 isotopes (1.00235
     Da), with the precursor's m/z and the spectrum's retention time.
   - **Fragments**: the N- and C-terminal fragment masses of the
     proteoform (the ion types of the PrSM's activation: b and y for CID
     and HCD, c and z-dot for ETD; residue masses of
     `res/base_data/amino_acid_base.xml`, the fixed modifications written
     in the proteoform string, the protein N-terminal acetylation) are
     matched to the spectrum's masses; each mass within `-t` (20 ppm) of
     a fragment mass is a calibrant, with the error to the nearest one,
     its m/z (at its charge) and the retention time. The proteoform
     string may hold only fixed modifications on single residues
     (`(C)[Carbamidomethyl]`) and an `[Acetyl]-` prefix; other PrSMs give
     no fragment calibrants.

   The console prints the PrSM numbers and the precursor and fragment
   calibrant numbers.
2. **Fits the models** (`Fitting the models`), each by robust least
   squares (10 rounds of Huber-type reweighting: weights
   sqrt(min(1, 1.5 s / |r|)) of the residuals r, s = 1.4826 x their
   median absolute value):
   - MS1, on the precursor errors: e = c0 + c_rt (rt - rt0) / 60;
   - MS/MS, on the fragment errors: e = c0 + c_mz (m/z - 800) / 1000 +
     c_rt (rt - rt0) / 60;

   e in ppm, rt in minutes, rt0 the median retention time of all
   calibrants, so c_rt is in ppm per hour and c_mz in ppm per 1000 m/z.
   At least 3 precursor and 4 fragment calibrants are needed. The console
   prints the coefficients and the robust standard deviations of the
   residuals.
3. **Writes the model** (below).

### `apply`

Corrects the masses of the msalign file spectrum by spectrum, at the
spectrum's `RETENTION_TIME`, by the relative error e of their m/z: a
corrected m/z is m/z (1 - e 1e-6), so a neutral mass M of charge z becomes
M - e 1e-6 (M + z x proton mass).

- **Precursors** (`PRECURSOR_MZ` and `PRECURSOR_MASS`, every precursor of
  a spectrum with its `PRECURSOR_CHARGE`): the MS1 model at the
  precursor's m/z.
- **Fragment masses**: the MS/MS model at the mass's m/z,
  (M + z x proton mass) / z.

The corrected values are written with six decimals, as TopFD writes them;
all other lines (the intensities, charges and EnvCNN scores of the masses,
the other header lines) are copied unchanged. A spectrum without
`RETENTION_TIME` is an error.

## 3. Output

- **`fit`**: the model file, a JSON object:

  | Key | Meaning |
  |---|---|
  | `mz0` | 800, the m/z of the MS/MS model's m/z term. |
  | `rt0_min` | rt0, minutes. |
  | `ms1`, `ms2` | The MS1 and MS/MS models: `c0` (ppm), `c_mz` (ppm per 1000 m/z; 0 for MS1), `c_rt_per_hour` (ppm per hour). |
  | `n_prec`, `n_frag`, `n_frag_prsms` | The precursor and fragment calibrant numbers and the PrSMs with fragment calibrants. |
  | `robust_sd_prec`, `robust_sd_frag` | 1.4826 x the median absolute residual of each fit, ppm. |
  | `note` | The model's formula. |

  `apply` reads only `mz0`, `rt0_min`, `ms1` and `ms2`.
- **`apply`**: the corrected msalign file. Give it the name TopPIC expects
  next to the run's other TopFD files (`<base>_ms2.msalign` in another
  directory, with `<base>_ms2.feature` linked or copied there and
  `<base>.sqlite` copied there), so that TopPIC finds the feature file and
  its post mass matching finds the database. Copy the database rather
  than link it: TopPIC writes its result tables into it, replacing those
  of an earlier search. The post mass matching uses the raw peaks of the
  database, which are not corrected; its fragment tolerance is TopPIC's
  `-e`.

## 4. Options

`topcalib fit`:

| Option | Default | Meaning |
|---|---|---|
| `-h`, `--help` | | Print the help message. |
| `-e`, `--max-evalue <number>` | 1e-5 | The largest E-value of a calibrant PrSM. |
| `-t`, `--tolerance <number>` | 20 | The tolerance (ppm) of the fragment matching. |

`topcalib apply` has only `-h`.

## 5. Examples

TopFD, a zero-shift search, the fit, the correction, and the search with a
5 ppm fragment tolerance (`-E`, the precursor tolerance, stays 10 ppm):

```sh
topfd sample.mzML
toppic -d -s 0 -D -t FDR -T FDR proteins.fasta sample_ms2.msalign
mkdir recal
topcalib fit sample_ms2_toppic_prsm_single.tsv sample_ms2.msalign recal/model.json
topcalib apply recal/model.json sample_ms2.msalign recal/sample_ms2.msalign
ln -s ../sample_ms2.feature recal/
cp sample.sqlite recal/
toppic -d -e 5 -t FDR -T FDR proteins.fasta recal/sample_ms2.msalign
```

The first search's own result files are written next to
`sample_ms2.msalign`; the second's into `recal/`.
