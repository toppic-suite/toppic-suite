# TopDiff manual

TopDiff (Top-down mass spectrometry-based identification of Differentially
expressed proteoforms) compares the abundances of proteoforms across several
samples. It takes the proteoform identifications that TopPIC or TopMG
produced for each sample, aligns the samples' retention times, matches the
same proteoform across samples, and writes one table with the abundance of
every proteoform in every sample.

```sh
topdiff [options] spectrum-file ...
```

`topdiff -h` prints the option list. The GUI tool `topdiff_gui` offers the
same function. Paths must not contain spaces.

## 1. Input

Each command-line argument is the MS/MS spectrum file of one sample, and at
least two samples are required. For every spectrum file TopDiff reads the
proteoform result file that the database search tool wrote for it. The
name is derived from the spectrum file: the extension and a trailing `_ms2`
or `_post_ms2` are removed, and TopDiff then looks for, in this order,

1. `<base>_post_ms2_<tool>_proteoform.xml`, the results of a TopPIC run
   with its default post mass matching, and
2. `<base>_ms2_<tool>_proteoform.xml`, the results of TopMG or of TopPIC
   run with `--disable-post-match`,

where `<tool>` is `toppic` (the default) or `topmg` (`--tool-name topmg`).
So for a sample deconvoluted by TopFD into `sample_ms2.msalign` and
searched by TopPIC, either `sample_ms2.msalign` or `sample_post_ms2.msalign`
can be given, and `sample_post_ms2_toppic_proteoform.xml` is read. TopDiff
prints the name of each result file it reads. When neither file exists, the
sample contributes no identifications and a warning "does not contain any
PrSM identifications" is printed, so check that the search finished for
every sample.

- **Run TopPIC/TopMG with the TopFD feature file** (the default; not with
  `--no-topfd-feature`). Abundances and retention times come from the
  proteoform features stored in the result file; without them the table has
  no usable abundances.

The samples should have been searched with the same settings and against
the same database, since proteoforms are matched across samples by protein
accession, precursor mass and aligned retention time.

## 2. What TopDiff does

1. **Reads** each sample's proteoform identifications: protein, residue
   range, proteoform sequence, precursor mass, spectrum id, the proteoform's
   cluster abundance, and the feature's retention-time range and apex.
2. **Normalizes** each sample's retention times by dividing them by the
   latest feature end time of the sample.
3. **Aligns** the retention times of samples 2 to n onto sample 1 with a
   dynamic-programming alignment, by apex time, of the 1,000 most abundant
   proteoforms of sample 1 and of the other sample (two positions are scored
   as a match when a proteoform of the other sample within 0.1 normalized
   time of that position has a precursor mass within `--error-tolerance` of
   the sample 1 proteoform), then warps every feature's times
   piecewise-linearly onto sample 1's scale.
4. **Matches proteoforms across samples.** PrSMs are visited in order of
   decreasing abundance; each PrSM whose proteoform (cluster) is not yet in
   the table starts a row, and in every other sample the most abundant
   unused proteoform with the same protein accession, a precursor mass
   within `--error-tolerance`, and an aligned apex time within 0.3 is taken
   as the same proteoform. Each proteoform cluster of a sample appears in at
   most one row.
5. **Writes** the table and prints the number of proteoform rows and how
   many of them were found in all samples.

## 3. Output

One tab-separated file, written into the directory of the **first** spectrum
file, named by `--output` (default `sample_diff.tsv`; give a bare file name,
not a path). Each row is one proteoform. The first columns describe it:

| Column | Meaning |
|---|---|
| `Protein accession` | Protein name of the identification. |
| `Protein description` | Protein description, in double quotes. |
| `First residue`, `Last residue` | First and last residue of the proteoform in the protein (1-based). |
| `Proteoform` | The proteoform sequence with its modifications. |
| `Precursor mass` | Adjusted precursor mass of the identification. |
| `# matched samples` | Number of samples in which the proteoform was found. |

They are followed by five columns per sample, in command-line order,
prefixed by the spectrum file name as typed:

| Column | Meaning |
|---|---|
| `<file> Abundance` | Abundance of the proteoform in that sample, as written by TopPIC or TopMG: the sum of the intensities of the distinct TopFD features in the proteoform's cluster, in scientific notation. |
| `<file> Spectrum id` | MS/MS spectrum id of the PrSM matched in that sample. |
| `<file> Retention time begin`, `<file> Retention time end` | Retention-time range of the feature in seconds, in the sample's own time scale (TopFD's `.feature` files list minutes). |
| `<file> Normalized time apex` | Apex time after normalization and alignment onto sample 1 (0 to 1). |

A sample in which the proteoform was not found has five empty fields.

## 4. Options

| Option | Default | Meaning |
|---|---|---|
| `-e`, `--error-tolerance <number>` | 1.2 | Precursor mass tolerance (Da) for matching proteoforms across samples and for the retention-time alignment. |
| `-t`, `--tool-name <toppic\|topmg>` | toppic | Which tool's result files to read. |
| `-o`, `--output <name>` | sample_diff.tsv | Name of the output table, written next to the first spectrum file. |

TopDiff stops with an error when fewer than two spectrum files are given
or when one of them does not exist.

## 5. Example

Three samples searched with TopPIC using its default settings:

```sh
topdiff -o abundance_diff.tsv sample_1_ms2.msalign sample_2_ms2.msalign sample_3_ms2.msalign
```

reads `sample_1_post_ms2_toppic_proteoform.xml` and its two counterparts
and writes `abundance_diff.tsv` next to `sample_1_ms2.msalign`. The same
samples searched with TopMG:

```sh
topdiff -t topmg sample_1_ms2.msalign sample_2_ms2.msalign sample_3_ms2.msalign
```
