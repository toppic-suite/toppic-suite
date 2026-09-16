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
least two samples are required. For every spectrum file TopDiff reads one
result file written by the database search tool, named after the spectrum
file:

| Spectrum file given | Result file read (`--tool-name toppic`, the default) | Result file read (`--tool-name topmg`) |
|---|---|---|
| `sample_ms2.msalign` | `sample_ms2_toppic_proteoform.xml` | `sample_ms2_topmg_proteoform.xml` |
| `sample_post_ms2.msalign` | `sample_post_ms2_toppic_proteoform.xml` | `sample_post_ms2_topmg_proteoform.xml` |

The rule is: strip the extension and a trailing `_ms2`, then append
`_ms2_<tool>_proteoform.xml`. Two consequences:

- **After a default TopPIC run, give TopDiff the `_post_ms2.msalign`
  files.** TopPIC's post mass matching (on by default) names all its result
  files after `sample_post_ms2.msalign`. If you pass the original
  `sample_ms2.msalign` instead, TopDiff looks for
  `sample_ms2_toppic_proteoform.xml`, does not find it, prints "does not
  contain any PrSM identifications" and produces an empty table. TopPIC runs
  with `--disable-post-match` and all TopMG runs keep the plain `_ms2` names.
- **Run TopPIC/TopMG with the TopFD feature file** (the default; not with
  `--no-topfd-feature`). Abundances and retention times come from the
  proteoform features stored in the result file; without them the table has
  no usable abundances.

The samples should have been searched with the same settings and against
the same database, since proteoforms are matched across samples by protein
accession, precursor mass and aligned retention time.

## 2. What TopDiff does

1. **Reads** each sample's proteoform identifications: protein, residue
   range, proteoform sequence, precursor mass, spectrum id, and the feature's
   abundance and retention-time range and apex.
2. **Normalizes** each sample's retention times to the 0 to 1 range of its
   own run.
3. **Aligns** the retention times of samples 2 to n onto sample 1 with a
   dynamic-programming alignment of the 1,000 most abundant features of each
   pair (two features can be paired when their precursor masses agree within
   `--error-tolerance` and their normalized times within 0.1), then warps
   every feature's times piecewise-linearly onto sample 1's scale.
4. **Matches proteoforms across samples.** Features are visited in order of
   decreasing abundance; each identified feature not yet used starts a table
   row, and in every other sample the feature with the same protein
   accession, a precursor mass within `--error-tolerance`, and an aligned
   apex time within 0.3 is taken as the same proteoform. Every feature is
   used at most once.
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
| `<file> Abundance` | Abundance (intensity) of the proteoform feature, in scientific notation. |
| `<file> Spectrum id` | MS/MS spectrum id of the identification; empty when the feature was matched but not identified in that sample. |
| `<file> Retention time begin`, `<file> Retention time end` | Retention-time range of the feature, in the sample's own time scale. |
| `<file> Normalized time apex` | Apex time after normalization and alignment onto sample 1 (0 to 1). |

A sample in which the proteoform was not found has five empty fields.

## 4. Options

| Option | Default | Meaning |
|---|---|---|
| `-e`, `--error-tolerance <number>` | 1.2 | Precursor mass tolerance (Da) for matching proteoforms across samples and for the retention-time alignment. |
| `-t`, `--tool-name <toppic\|topmg>` | toppic | Which tool's result files to read. |
| `-o`, `--output <name>` | sample_diff.tsv | Name of the output table, written next to the first spectrum file. |

TopDiff has no hidden options. It stops with an error when fewer than two
spectrum files are given or when one of them does not exist.

## 5. Example

Three samples searched with TopPIC using its default settings:

```sh
topdiff -o abundance_diff.tsv \
    sample_1_post_ms2.msalign sample_2_post_ms2.msalign sample_3_post_ms2.msalign
```

writes `abundance_diff.tsv` next to `sample_1_post_ms2.msalign`. The same
samples searched with TopMG:

```sh
topdiff -t topmg sample_1_ms2.msalign sample_2_ms2.msalign sample_3_ms2.msalign
```
