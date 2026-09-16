# TopPIC manual

TopPIC (Top-down mass spectrometry-based Proteoform Identification and
Characterization) identifies proteoforms by searching deconvoluted top-down
MS/MS spectra against a protein sequence database. It finds proteoforms with
known modifications, with up to two unexpected mass shifts (mutations,
unknown modifications), and with N-terminal truncations or forms, estimates
the statistical significance of every identification, and characterizes
unexpected shifts with common modifications.

```sh
toppic [options] database-file spectrum-file ...
```

`toppic -h` prints the option list. The GUI tool `toppic_gui` offers the same
function. Paths must not contain spaces or quotation marks, and the database
and spectrum paths must be shorter than 200 characters.

## 1. Input

**Protein database.** A FASTA file (`.fasta` or `.fa`). Every protein needs
a unique accession; duplicated accessions stop the run.

**Spectrum files.** One or more `_ms2.msalign` files written by TopFD (or
by TopDIA, whose pseudo-MS/MS spectra use the same format). Several files
given on one command line are searched one after another with the same
parameters.

**Companion files from TopFD.** For a spectrum file `sample_ms2.msalign`,
TopPIC also reads, from the same directory:

| File | Needed when | Purpose |
|---|---|---|
| `sample_ms2.feature` | Always, unless `--no-topfd-feature` is given | The proteoform feature of each MS/MS spectrum, used to cluster the identifications of one proteoform across spectra and to report abundances. |
| `sample.sqlite` | Always, unless `--disable-post-match` is given | The centroided MS/MS peaks, used by post mass matching (section 2, step 10). TopFD writes it unless run with `-N`. A missing database is an error, reported before the search starts. |

**Database indexes.** TopPIC does not require `topindex`. It prepares the
database itself in a `<database>_idx/` directory next to the FASTA file, and
builds the filtering indexes in memory for each run. Running `topindex` first
with the same `-f`, `-n`, `-e` and `-d` settings stores those indexes on disk
and saves time when the same database is searched repeatedly (see the
[TopIndex manual](topindex_manual.md)).

## 2. What TopPIC does

For each spectrum file TopPIC prints its parameters and runs these steps.
Intermediate files carry the spectrum file name plus a `.toppic_*`
extension and are deleted at the end unless `--keep-temp-files` is given.

1. **Zero unexpected shift filtering and search.** Protein candidates are
   selected with index-based filters and aligned to each spectrum without any
   unexpected mass shift, for complete, N-terminal, C-terminal and internal
   proteoforms (`--proteoform-type`).
2. **Variable PTM filtering and search** (only with `--variable-ptm-file-name`
   and `--variable-ptm-num` of at least 1): the same, allowing the listed
   variable modifications.
3. **One unexpected shift filtering and search** (with `--num-shift` of at
   least 1): alignment allowing one unexpected mass shift between
   `--min-shift` and `--max-shift`.
4. **Multiple unexpected shifts filtering and search** (with `--num-shift`
   2): alignment allowing two unexpected shifts.
5. **Merging PrSMs.** The proteoform-spectrum matches (PrSMs) of all searches
   are merged.
6. **E-value computation.** A generating-function method estimates the
   E-value of each PrSM.
7. **Top PrSM selection.** The best PrSM of each spectrum is kept.
8. **Recounting matched masses and fragments.** The search uses only the
   masses whose EnvCNN score is at least `--filter-by-env-cnn`; matched
   masses and fragments are recounted against the complete spectra.
9. **Post mass matching** (default; off with `--disable-post-match`). For
   each PrSM, the theoretical fragment masses that no deconvoluted mass
   matched are searched in the centroided MS/MS peaks of the spectrum, stored
   in the TopFD SQLite database, using the isotopic-envelope test of
   MSPathFinderT: an envelope is accepted when its most abundant peak and at
   least `--post-min-peak-num` of its peaks are present and its intensities
   fit the theoretical ones. Each accepted fragment adds one mass to the
   spectrum, scored by EnvCNN. The spectra with the added masses are written
   to `sample_post_ms2.msalign`, and **all following steps and result files
   use that name**. The added envelopes are also stored in the database.
10. **Proteoform and protein clustering.** PrSMs of the same proteoform are
    grouped by their TopFD feature and by precursor mass within
    `--proteoform-error-tolerance`, then grouped by protein.
11. **PTM characterization** (only with `--local-ptm-file-name`). Each
    unexpected mass shift is compared with the listed common modifications,
    and the modification and its site are reported when the MIScore reaches
    `--miscore-threshold`.
12. **FDR computation** (only with `--decoy`). Spectrum-, proteoform- and
    protein-level false discovery rates are estimated from the decoy hits.
13. **Filtering and output at three levels.** PrSMs are filtered by the
    spectrum-level cutoff and written to the PrSM tables; the best PrSM per
    proteoform is filtered by the proteoform-level cutoff and written to the
    proteoform tables; the best PrSM per protein is filtered by the
    protein-level cutoff and written to the protein tables.
14. **Writing identifications to the SQLite database**, for visualization.

## 3. Output files

For an input `sample_ms2.msalign` searched with the default settings, the
result files are named after `sample_post_ms2`; with `--disable-post-match`
they are named after `sample_ms2`:

| File | Contents |
|---|---|
| `sample_post_ms2.msalign` | The spectra with the masses added by post mass matching. The input file is left unchanged. |
| `sample_post_ms2_toppic_prsm.tsv`, `..._prsm_single.tsv` | All PrSMs that pass the spectrum-level cutoff. |
| `sample_post_ms2_toppic_proteoform.tsv`, `..._proteoform_single.tsv` | One PrSM per identified proteoform (the best one), passing the proteoform-level cutoff. |
| `sample_post_ms2_toppic_protein.tsv`, `..._protein_single.tsv` | One PrSM per identified protein (the best one), passing the protein-level cutoff. |
| `sample_post_ms2_toppic_prsm.xml`, `..._proteoform.xml`, `..._protein.xml` | The same three sets of PrSMs in XML, read by TopDiff and the visualization tools. |
| `sample.sqlite` | Updated with the tables `prsm`, `proteoform`, `protein`, `prsm_mass_shift`, `prsm_protein_match` and `fasta_seq`, for visualization. |
| `sample_ms2.toppic_raw_prsm` | The PrSMs before recounting and cutoffs; kept for combining fractions. |

The `_single` tables have one line per PrSM. The tables without `_single`
add, after each PrSM, one line for every other protein that contains the
same proteoform sequence, with only the protein columns filled. Each table
starts with the parameters of the run, the numbers of identified PrSMs,
proteoforms and proteins, and then the header line. The columns are:

| Column | Meaning |
|---|---|
| `Data file name` | The spectrum file. |
| `Prsm ID`, `Spectrum ID`, `Fragmentation`, `Scan(s)`, `Retention time` | The PrSM and its spectrum (several values when spectra are combined, see `--num-combined-spectra`); retention time in minutes. |
| `# Masses`, `Charge`, `Precursor mass`, `Adjusted precursor mass` | Number of deconvoluted masses, precursor charge, precursor mass, and the precursor mass adjusted to the proteoform. |
| `Proteoform ID`, `Proteoform intensity` | The proteoform cluster and its abundance. |
| `Feature ID`, `Feature intensity`, `Feature score`, `Feature apex time` | The TopFD proteoform feature of the spectrum: id, abundance, ECScore and apex retention time (minutes). |
| `# Protein hits`, `Protein accession`, `Protein description` | Number of proteins containing the proteoform, and the reported protein. |
| `First residue`, `Last residue`, `Special amino acids`, `Database protein sequence`, `Previous amino acid`, `Proteoform`, `Next amino acid`, `Proteoform mass` | Position of the proteoform in the protein, the protein sequence, the proteoform with its modifications, and its mass. |
| `Protein N-terminal form`, `Fixed PTMs`, `# Unexpected modifications`, `Unexpected modifications`, `# Variable PTMs`, `Variable PTMs`, `MIScore` | The N-terminal form, the modifications of the proteoform, and the modification identification score of the PTM characterization. |
| `# Matched masses`, `# Matched fragments`, `E-value` | The match statistics and the E-value. |
| `Spectrum-level Q-value`, `Proteoform-level Q-value`, `Protein-level Q-value` | The FDR-based q-values (only with `--decoy`; `-` otherwise). |

A `-` marks a value that does not apply.

## 4. Options

Database and modifications:

| Option | Default | Meaning |
|---|---|---|
| `-f`, `--fixed-mod <C57\|C58\|file>` | none | Fixed modifications: `C57` carbamidomethylation of cysteine, `C58` carboxymethylation of cysteine, or a modification file (section 5). |
| `-n`, `--n-terminal-form <list>` | NONE,NME,NME_ACETYLATION,M_ACETYLATION | Allowed protein N-terminal forms, comma-separated: no modification, methionine excision, acetylation after methionine excision, acetylation of the initial methionine. |
| `-R`, `--proteoform-type <list>` | COMPLETE,PREFIX,SUFFIX,INTERNAL | Allowed proteoform types, comma-separated: the complete protein, its N-terminal part, its C-terminal part, or an internal part. |
| `-d`, `--decoy` | off | Search a shuffled decoy database as well, to estimate false discovery rates. Required for the `FDR` cutoff types. |

Mass shifts and variable modifications:

| Option | Default | Meaning |
|---|---|---|
| `-s`, `--num-shift <0\|1\|2>` | 1 | Maximum number of unexpected mass shifts in a proteoform. |
| `-m`, `--min-shift <number>` | -500 | Minimum mass of an unexpected shift (Da). |
| `-M`, `--max-shift <number>` | 500 | Maximum mass of an unexpected shift (Da). |
| `-b`, `--variable-ptm-file-name <file>` | none | A modification file (section 5) with the variable modifications to search. |
| `-S`, `--variable-ptm-num <int>` | 3 | Maximum number of variable modifications in a proteoform. |
| `-B`, `--local-ptm-file-name <file>` | none | A modification file with the common modifications used to characterize unexpected shifts. |
| `-H`, `--miscore-threshold <0..1>` | 0.15 | Minimum MIScore for reporting a characterized modification. |

Spectra and tolerances:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|HCD\|UVPD\|FILE>` | FILE | Fragmentation method; `FILE` takes it from each spectrum. |
| `-e`, `--mass-error-tolerance <int>` | 10 | Error tolerance of precursor and fragment masses (ppm). |
| `-F`, `--filter-by-env-cnn <0..1>` | 0.2 | Masses with an EnvCNN score below this are ignored during filtering, search and E-value computation. |
| `-p`, `--proteoform-error-tolerance <number>` | 1.2 | Precursor mass tolerance (Da) for grouping PrSMs into proteoform clusters. |
| `-r`, `--num-combined-spectra <int>` | 1 | Number of consecutive MS/MS spectra of one precursor searched together, for alternating fragmentation (2 for pairs, 3 for triplets). |
| `-A`, `--approximate-spectra` | off | Use approximate spectra in protein filtering, for higher sensitivity. |
| `-x`, `--no-topfd-feature` | off | Run without the TopFD feature file; proteoforms are then clustered by precursor mass only and no abundances are reported. |

Post mass matching:

| Option | Default | Meaning |
|---|---|---|
| `-E`, `--disable-post-match` | off | Do not run post mass matching; result files keep the `_ms2` name and the SQLite database is not needed. |
| `-I`, `--post-min-peak-num <int>` | 1 | Minimum number of observed isotopic peaks of a fragment accepted by post mass matching. |

Cutoffs and output:

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--spectrum-cutoff-type <EVALUE\|FDR>`, `-v`, `--spectrum-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for PrSMs. |
| `-T`, `--proteoform-cutoff-type <EVALUE\|FDR>`, `-V`, `--proteoform-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteoforms. |
| `-y`, `--protein-cutoff-type <EVALUE\|FDR>`, `-Y`, `--protein-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteins; with `FDR` a protein is represented by its best proteoform. |
| `-c`, `--combined-file-name <name>` | none | Combine several fractions: after the files are searched, their spectra, features and PrSMs are merged under `<name>_ms2.msalign` and the post-search steps are repeated on the merged data. See the note below. |
| `-u`, `--thread-number <int>` | 1 | Number of threads. The filtering steps may use fewer, depending on the available memory. |
| `-k`, `--keep-temp-files` | off | Keep the intermediate files. |
| `-K`, `--keep-decoy-ids` | off | Keep decoy identifications in the result tables. |

Advanced options (accepted but not shown by `-h`):

| Option | Meaning |
|---|---|
| `-P`, `--proteoform-ppm-error` | Interpret `--proteoform-error-tolerance` in ppm instead of Da. |
| `-U`, `--top-prsm-number <int>` | Number of best PrSMs kept per spectrum (default 1). |
| `-C`, `--combine-result-only` | With `-c`: skip the searches and only merge the existing results of the fractions. |
| `-o`, `--output-raw-prsm-table` | Also write the PrSM tables before the cutoffs (`_toppic_raw_prsm.tsv`, `_toppic_raw_prsm_single.tsv`). |
| `-O`, `--output-prsm-coverage` | Write a `.toppic_prsm_coverage` file with the fragment coverage of each PrSM. |
| `--filtering-result-number <int>` | Number of candidates kept by the multiple-shift filter (default 20). |

**Combining fractions.** With `-c`, the merged spectrum file
`<name>_ms2.msalign` has no TopFD SQLite database, so post mass matching
cannot run on it. Use `-c` together with `--disable-post-match`.

## 5. Modification files

The files given to `-f`, `-b` and `-B` are plain text with one modification
per line and five comma-separated fields:

```text
# Name,Mass,Residues,Position,UnimodID
Phospho,79.966331,STY,any,21
Acetyl,42.010565,K,any,1
Oxidation,15.994915,M,any,35
```

`Residues` lists the one-letter codes of the residues that can carry the
modification, or `*` for any residue; `Position` is `any`, `N-term` or
`C-term` (`*` and `any` cannot be combined); `UnimodID` is the Unimod
accession, or `-1` when unknown. Lines starting with `#` are comments.

## 6. Examples

```sh
# Default search: one unexpected shift, E-value cutoffs, post mass matching
toppic proteins.fasta sample_ms2.msalign

# Target-decoy search with 1% FDR at all three levels, 8 threads
toppic -d -t FDR -v 0.01 -T FDR -V 0.01 -y FDR -Y 0.01 -u 8 proteins.fasta sample_ms2.msalign

# Carbamidomethylated cysteines, phosphorylation as a variable PTM, no unexpected shift
toppic -f C57 -b phospho.txt -s 0 proteins.fasta sample_ms2.msalign

# Two unexpected shifts characterized with a list of common modifications
toppic -s 2 -B common_mods.txt proteins.fasta sample_ms2.msalign

# Three fractions searched and combined (post mass matching off, see section 4)
toppic -E -c combined proteins.fasta frac1_ms2.msalign frac2_ms2.msalign frac3_ms2.msalign
```

The identifications of the first example are in
`sample_post_ms2_toppic_prsm_single.tsv` (all PrSMs),
`sample_post_ms2_toppic_proteoform_single.tsv` (one per proteoform) and
`sample_post_ms2_toppic_protein_single.tsv` (one per protein).
