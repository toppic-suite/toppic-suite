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
function. Avoid spaces and quotation marks in paths (the directory of the
`toppic` executable must not contain spaces), and the database and spectrum
paths must be at most 200 characters long.

## 1. Input

**Protein database.** A FASTA file (`.fasta`, `.fa`, `.FASTA` or `.FA`).
If two records share an accession (the header text up to the first
space), the tool writes `<database>_dedup.fasta` next to the input,
keeping the first record of each accession, prints a notice, and searches
that file instead; its indexes go to `<database>_dedup.fasta_idx/`.

**Spectrum files.** One or more `_ms2.msalign` files written by TopFD (or
by TopDIA, whose pseudo-MS/MS spectra use the same format). Several files
given on one command line are searched one after another with the same
parameters. A file written by `topfd --average-mass`, whose precursor and
fragment masses are average masses, must be searched with `--average-mass`
(section 2).

**Companion files from TopFD.** For a spectrum file `sample_ms2.msalign`,
TopPIC also reads, from the same directory:

| File | Needed when | Purpose |
|---|---|---|
| `sample_ms2.feature` | Always, unless `--no-topfd-feature` is given | The proteoform feature of each MS/MS spectrum, used to cluster the identifications of one proteoform across spectra and to report abundances. |
| `sample.sqlite` | Always, unless `--disable-post-match` is given | The centroided MS/MS peaks, used by post mass matching (section 2, step 11). TopFD writes it unless run with `-N`. A missing database is an error, reported before the search starts. TopPIC only reads it: it searches its own copy, `sample_toppic.sqlite` (section 2). |

**Database indexes.** TopPIC does not require `topindex`. It prepares the
database itself in a `<database>_idx/` directory next to the FASTA file, and
builds the filtering indexes in memory for each run. Running `topindex` first
with the same `-f`, `-n`, `-e` and `-d` settings stores those indexes on disk
and saves time when the same database is searched repeatedly (see the
[TopIndex manual](topindex_manual.md)).

## 2. What TopPIC does

For each spectrum file TopPIC prints its parameters and runs these steps.
Intermediate files carry the spectrum file name plus a `.toppic_*`
extension and are deleted at the end unless `--keep-temp-files` is given,
except the `toppic_record_prsm` file (section 3), which is always kept.

**TopPIC's database.** Before a spectrum file `sample_ms2.msalign` is
searched, TopPIC writes its own SQLite database, `sample_toppic.sqlite`: a
copy of TopFD's `sample.sqlite` with all its tables (spectra, peaks,
deconvoluted masses and their envelopes, and the other tables TopFD
wrote). An earlier `sample_toppic.sqlite` is deleted first, so every run
starts from TopFD's data. Post mass matching reads the peaks of the copy
and adds its envelopes to it (step 11), and the identifications are
written into it (step 17); TopFD's database is only read. The name is the
same with every option (`--average-mass`, `--disable-post-match`, `-c`),
so a TopFD and a TopPIC run leave two databases per spectrum file. Without
TopFD's database (TopFD `-N`, possible with `--disable-post-match`) no
copy is written and the identifications are not written to a database.

With `--average-mass`, every spectrum file is first converted: the average
precursor and fragment masses of `sample_ms2.msalign` are converted back to
monoisotopic masses with the theoretical isotopic distributions TopFD used
(the inverse of its conversion), the spectra are written to
`sample_mono_ms2.msalign` (the input file is left unchanged, a comment line
records the conversion), `sample_ms2.feature` is copied to
`sample_mono_ms2.feature`, and the steps below, the result file names and
the merge of a `-c` run all use the converted file. TopFD keeps the
monoisotopic masses in its database, so `sample_toppic.sqlite` is the
same copy as without the option. The result files are named after
`sample_mono_post_ms2`, or `sample_mono_ms2` with `--disable-post-match`.

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

   With `--variable-ptm-file-name` the filters of steps 1-3 run as one step
   over the index files, printed as `Zero unexpected shift, variable PTM and
   one unexpected shift filtering` (or `Zero unexpected shift and variable
   PTM filtering` with `--num-shift 0`), and the searches follow in the order
   above; the candidates are the same as with separate filters.
4. **Merging PrSMs.** The proteoform-spectrum matches (PrSMs) of the
   searches above are merged.
5. **E-value computation.** A generating-function method estimates the
   E-value of each PrSM.
6. **Top PrSM selection.** The best PrSM of each spectrum is kept.
7. **Multiple unexpected shifts filtering and search** (with `--num-shift`
   2): alignment allowing two unexpected shifts, only for the spectra whose
   best PrSM so far has an E-value at or above the smallest E-value of a
   decoy PrSM (every spectrum a decoy could compete with is searched
   again; 0.01 without `--decoy`); the spectra with a better PrSM keep it. The new PrSMs are merged, get E-values, and,
   together with the earlier PrSMs of the searched spectra, the best PrSM per
   spectrum is kept.
8. **Proteoform and protein clustering.** PrSMs whose spectrum has no
   TopFD feature are discarded (unless `--no-topfd-feature`). PrSMs sharing
   a TopFD feature form one proteoform; proteoforms of the same protein
   whose adjusted precursor masses agree within
   `--proteoform-error-tolerance` are merged; proteoforms are then grouped
   by protein.
9. **Cluster protein correction.** In every proteoform cluster the PrSMs
   whose protein differs from the protein of the cluster's PrSM with the
   lowest E-value are given that PrSM's proteoform, so that a proteoform
   mass matched to two proteins is reported once, for the better-supported
   protein; their E-values are recomputed for the new proteoform on the
   masses the search used.
10. **Recounting matched masses and fragments.** The search and the
    cluster protein correction use only the masses whose EnvCNN score is at
    least `--filter-by-env-cnn`; the matched masses and fragments of every
    PrSM are recounted against the complete spectra.
11. **Post mass matching** (default; off with `--disable-post-match`). For
    the best PrSM of every spectrum, target or decoy, after the clustering,
    the cluster protein correction and the recount, the theoretical
    fragment masses that no deconvoluted mass matched are searched in the
    centroided MS/MS peaks of the spectrum, stored in TopPIC's copy of the
    TopFD SQLite database (`sample_toppic.sqlite`), using the isotopic-envelope test of MSPathFinderT: an envelope
    is accepted when its most abundant peak and at least
    `--post-min-peak-num` of its peaks are present and its intensities fit
    the theoretical ones. Each accepted fragment adds one mass to the
    spectrum, scored by EnvCNN. The spectra with the added masses are
    written to `sample_post_ms2.msalign`, **the steps below and the result
    files use that name**, and the added envelopes are also stored in
    `sample_toppic.sqlite`. The E-values and the clusters are not changed, so without
    `--decoy` the cutoffs and the identifications are those of a run with
    `--disable-post-match` and only the matched masses and fragments grow;
    with `--decoy` the rescoring (step 13) uses the matched masses and
    fragments, so the FDRs can differ slightly.
12. **PTM characterization** (only with `--local-ptm-file-name`). Each
    unexpected mass shift is compared with the listed common modifications,
    and the modification and its site are reported when the MIScore reaches
    `--miscore-threshold`.
13. **Rescoring** (only with `--decoy`). The best PrSM of every spectrum,
    target or decoy, gets a score from a linear support vector machine
    trained to separate the targets at 1% FDR from the decoys (the method
    of Percolator, Kall et al., Nat. Methods 2007), with 3-fold
    cross-validation: the spectra are split into three parts, and the PrSMs
    of each part are scored by a model trained on the other two. The
    features describe the PrSM and its spectrum: the E-value, the matched
    masses and fragments, the unexpected shifts and variable PTMs, the
    N-terminal form, the proteoform's length and mass, the precursor, the
    TopFD feature, the fragment mass errors, the cleavage sites covered and
    the retention time against a value predicted from the sequence. The
    E-values are not changed. When a part has fewer than 10 targets at 1%
    FDR to train on, every PrSM keeps its E-value as the score.
14. **FDR computation** (only with `--decoy`). Spectrum-, proteoform- and
    protein-level false discovery rates are estimated from the decoy hits,
    with the PrSMs ranked by the score of step 13, so the `FDR` cutoffs
    apply to it; the `EVALUE` cutoffs still compare E-values.
15. **Record file and decoy removal.** The PrSMs with their FDRs are
    kept as `sample_post_ms2.toppic_record_prsm`, decoy identifications
    included. Unless `--keep-decoy-ids` is given, the decoy identifications
    are then removed from the PrSMs.
16. **Filtering and output at three levels.** The PrSMs passing the
    spectrum-level cutoff are written to the PrSM tables; the best PrSM per
    proteoform passing the proteoform-level cutoff to the proteoform tables;
    the best PrSM per protein passing the protein-level cutoff to the
    protein tables.
17. **Writing identifications to the SQLite database**
    `sample_toppic.sqlite`, for visualization.

## 3. Output files

For an input `sample_ms2.msalign` searched with the default settings, the
result files are named after `sample_post_ms2`; with `--disable-post-match`
they are named after `sample_ms2`, and the combined results of a `-c` run
after `<name>_ms2`. With `--average-mass` the names start with
`sample_mono` instead of `sample` (section 2):

| File | Contents |
|---|---|
| `sample_post_ms2.msalign` | The spectra with the masses added by post mass matching. The input file is left unchanged. |
| `sample_post_ms2_toppic_prsm.tsv`, `..._prsm_single.tsv` | All PrSMs that pass the spectrum-level cutoff. |
| `sample_post_ms2_toppic_proteoform.tsv`, `..._proteoform_single.tsv` | One PrSM per identified proteoform (the best one), passing the proteoform-level cutoff. |
| `sample_post_ms2_toppic_protein.tsv`, `..._protein_single.tsv` | One PrSM per identified protein (the best one), passing the protein-level cutoff. |
| `sample_post_ms2_toppic_prsm.xml`, `..._proteoform.xml`, `..._protein.xml` | The same three sets of PrSMs in XML, read by TopDiff and the visualization tools. Every PrSM records the TopPIC Suite version that wrote it, and TopDiff reads only the files of its own version. Every proteoform element has its classification, `<classification>` (the `Proteoform classification` column). |
| `sample_toppic.sqlite` (the same name with every option) | TopPIC's database (section 2): a copy of TopFD's `sample.sqlite` with the envelopes added by post mass matching and the tables `prsm`, `proteoform`, `protein`, `prsm_mass_shift`, `prsm_protein_match` and `fasta_seq`, for visualization. The first three hold the PrSMs of the tables below, with the columns `precursor_id` and `proforma` for `Precursor ID` and `ProForma`. It is written again by every run on the spectrum file; TopFD's `sample.sqlite` is not changed. |
| `sample_post_ms2.toppic_record_prsm` | The best PrSM of every spectrum after the FDR computation and before the cutoffs, decoy identifications included. Always kept. |

The `_single` tables have one line per PrSM. The tables without `_single`
add, after each PrSM, one line for every other protein that contains the
same proteoform sequence, with only the data file name, PrSM ID, protein
accession and description, first and last residue and special amino acids
filled. Each table starts with the parameters of the run, the numbers of
identified PrSMs, proteoforms and proteins, and then the header line. The
columns are:

| Column | Meaning |
|---|---|
| `Data file name` | The spectrum file. |
| `Prsm ID`, `Spectrum ID`, `Precursor ID`, `Fragmentation`, `Scan(s)`, `Retention time` | The PrSM and its spectrum (several values when spectra are combined, see `--num-combined-spectra`); `Precursor ID` is the position of the PrSM's precursor in the spectrum's precursor list of the msalign file (0, the first); retention time in minutes. |
| `# Masses`, `Charge`, `Precursor mass`, `Adjusted precursor mass` | Number of deconvoluted masses, precursor charge, precursor mass, and the precursor mass adjusted to the proteoform. |
| `N-terminal label mass` | The mass of the N-terminal label (`-N`, for example TMT) that the search assumed on every proteoform, 0 without a label. The adjusted precursor mass, the proteoform and its mass leave the label out; the precursor mass includes it. |
| `Proteoform ID`, `Proteoform intensity` | The proteoform cluster and its abundance. |
| `Feature ID`, `Feature intensity`, `Feature score`, `Feature apex time` | The TopFD proteoform feature of the spectrum: id, abundance, ECScore and apex retention time (minutes). |
| `# Protein hits`, `Protein accession`, `Protein description` | Number of proteins containing the proteoform, and the reported protein. |
| `First residue`, `Last residue`, `Special amino acids`, `Database protein sequence`, `Previous amino acid`, `Proteoform`, `Next amino acid` | Position of the proteoform in the protein, the protein sequence, and the proteoform with its modifications (the annotation is explained in [Proteoform annotation](proteoform.md)). |
| `ProForma` | The proteoform in the ProForma 2.1 notation of the HUPO-PSI (https://github.com/HUPO-PSI/ProForma), with the PTMs as Unimod names and accessions, for example `[Acetyl|UNIMOD:1]-SGREGM[Oxidation|UNIMOD:35]DEE...` (see [Proteoform annotation](proteoform.md), section 6). |
| `Proteoform mass` | The monoisotopic mass of the proteoform. |
| `Protein N-terminal form`, `Fixed PTMs`, `# Unexpected modifications`, `Unexpected modifications`, `# Variable PTMs`, `Variable PTMs`, `MIScore` | The N-terminal form, the modifications of the proteoform as `name:[position]` lists (see [Proteoform annotation](proteoform.md)), and the modification identification score of the PTM characterization. |
| `# Matched exp fragments`, `# Matched theo fragments`, `E-value` | The match statistics and the E-value. |
| `Spectrum-level Q-value`, `Proteoform-level Q-value`, `Protein-level Q-value` | The FDR-based q-values (only with `--decoy`; `-` otherwise). |
| `Proteoform classification` | The level of the proteoform in the five-level classification of proteoform identifications (Smith et al., Nat. Methods 16:939-940, 2019), from the unexpected modifications, the protein hits and the variable PTM sites. Without unexpected modifications: `1` in one protein (`# Protein hits` 1), `2D` in several proteins; `2A` and `3` instead when a variable PTM site is undecided (a variable PTM on a range of residues with more residues that can carry it than PTMs, see `Variable PTMs`). With unexpected modifications: `4` in one protein, `5` in several. Protein N-terminal forms and fixed PTMs do not change the level. |

A `-` marks a value that does not apply.

## 4. Options

Database and modifications:

| Option | Default | Meaning |
|---|---|---|
| `-f`, `--fixed-mod <C57\|C58\|file>` | none | Fixed modifications: `C57` carbamidomethylation of cysteine, `C58` carboxymethylation of cysteine, or a modification file (section 5). |
| `-n`, `--n-terminal-form <list>` | NONE,NME,NME_ACETYLATION,M_ACETYLATION | Allowed protein N-terminal forms, comma-separated: no modification, methionine excision, acetylation after methionine excision, acetylation of the initial methionine. |
| `-N`, `--n-terminal-label <number>` | 0 | Mass (Da) of a label on the N-terminus of every proteoform, for example 229.162932 for TMT. The label is added to the N-terminal fragment ions and removed from the precursor mass; the reported proteoforms and adjusted precursor masses are label-free, and the `N-terminal label mass` column gives the label (section 3). Every N-terminus is labeled, acetylated ones included. In `toppic_gui`: "N-terminal label mass (Da)" on the Advanced Parameters tab. |
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
| `-a`, `--activation <CID\|ETD\|ETCID\|ETHCD\|HCD\|UVPD\|FILE>` | FILE | Fragmentation method; `FILE` takes it from each spectrum. ETCID (ETciD) and ETHCD (EThcD) are searched with the ion types of ETD (c, z-dot). |
| `-e`, `--frag-error-tolerance <int>` | 10 | Error tolerance of fragment masses (ppm). |
| `-E`, `--prec-error-tolerance <int>` | 10 | Error tolerance of precursor masses (ppm). |
| `-F`, `--filter-by-env-cnn <0..1>` | 0.4 | Masses with an EnvCNN score below this are ignored during filtering, search and E-value computation. |
| `-p`, `--proteoform-error-tolerance <number>` | 1.2 | Precursor mass tolerance (Da) for grouping PrSMs into proteoform clusters. |
| `-r`, `--num-combined-spectra <int>` | 1 | Number of consecutive MS/MS spectra of one precursor searched together, for alternating fragmentation (2 for pairs, 3 for triplets). |
| `-A`, `--approximate-spectra` | off | Use approximate spectra in protein filtering, for higher sensitivity. Only affects the variable PTM filtering step (`-b`). |
| `-x`, `--no-topfd-feature` | off | Run without the TopFD feature file; proteoforms are then clustered by protein, precursor mass and sequence overlap (two PrSMs share a proteoform only if their proteoform sequences overlap by at least 70% of the shorter one, so equal masses at different positions of a protein stay apart), the `Feature ...` columns are `-` and the proteoform intensity is `-1`. |
| `-g`, `--average-mass` | off | The spectrum files report average masses (`topfd --average-mass`): each is converted to `<name>_mono_ms2.msalign` with monoisotopic masses before the search (section 2). In `toppic_gui`, the checkbox "Spectrum files with average masses" on the Advanced Parameters tab turns this on. |

Post mass matching:

| Option | Default | Meaning |
|---|---|---|
| `-D`, `--disable-post-match` | off | Do not run post mass matching; result files keep the `_ms2` name and TopFD's SQLite database is not needed (without it no `<name>_toppic.sqlite` is written). |
| `-I`, `--post-min-peak-num <int>` | 1 | Minimum number of observed isotopic peaks of a fragment accepted by post mass matching (at least 1). |

Cutoffs and output:

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--spectrum-cutoff-type <EVALUE\|FDR>`, `-v`, `--spectrum-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for PrSMs. |
| `-T`, `--proteoform-cutoff-type <EVALUE\|FDR>`, `-V`, `--proteoform-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteoforms. |
| `-y`, `--protein-cutoff-type <EVALUE\|FDR>`, `-Y`, `--protein-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteins; with `FDR` a protein is represented by its best proteoform. |
| `-c`, `--combined-file-name <name>` | none | Combine several fractions: after the files are searched, their spectra and best PrSMs per spectrum are merged under `<name>_ms2`, and the clustering (by protein, precursor mass and sequence overlap), rescoring, FDR computation, decoy removal, cutoffs and tables are repeated on the merged data; the files' results are then written with the proteoform ids of the combined results. A relative `<name>` is placed in the directory of the first spectrum file and must differ from the input file names. See the note below. |
| `-u`, `--thread-number <int>` | 1 | Number of threads; must not exceed the number of hardware threads. Each protein filtering step may use fewer: it fits its thread number to the memory available when it starts (about 1.25 GiB per thread for zero-shift filtering, 1 GiB for variable-PTM and one-shift filtering, 0.5 GiB for multiple-shift filtering). The run starts with a warning listing the filtering steps that will use fewer threads and how many, and a step prints its number again if the memory available when it starts gives a different one. All other steps use the full number. |
| `-k`, `--keep-temp-files` | off | Keep the intermediate files. |
| `-K`, `--keep-decoy-ids` | off | Keep decoy identifications in the result tables. |

**Combining fractions.** With `-c`, each fraction is searched,
post-matched, clustered and given its FDRs on its own first (post mass
matching needs the fraction's own TopFD SQLite database, and each
fraction gets its own `<fraction>_toppic.sqlite`). The merge then
takes, for every fraction, the post-matched spectra
(`<fraction>_post_ms2.msalign`, or `<fraction>_ms2.msalign` with
`--disable-post-match`) and the best PrSM of every spectrum after the
fraction's FDR computation
(`<fraction>_post_ms2.toppic_record_prsm`, or
`<fraction>_ms2.toppic_record_prsm`; decoy identifications included). The
merged PrSMs are clustered again into proteoforms and proteins across the
fractions, by protein and precursor mass as with `-x` (the feature files
are not merged), rescored and given one set of FDRs (steps 8, 13 and 14), and the decoy
removal, the cutoffs and the output steps are repeated on them; the
recount, post mass matching and PTM characterization of the fractions are
kept. The combined results are always named `<name>_ms2_toppic_*`, the
merged spectra `<name>_ms2.msalign` and the merged best PrSMs
`<name>_ms2.toppic_record_prsm`. No SQLite database is written for the
combined results. The fractions' own results (tables, XML files and SQLite
databases) are written after the combined results, with the proteoform ids
of the combined results: every PrSM of a fraction carries the
`Proteoform ID` of its PrSM in the combined tables, so a proteoform has
one id in the combined tables and in the tables of every fraction it is
identified in, and the fractions' results can be compared by it. The
fractions' FDRs and cutoffs are their own; a proteoform of the combined
results may hold PrSMs that were two proteoforms in the fraction, or the
other way round, and the fraction's proteoform table lists it accordingly.

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
accession, or `-1` when unknown. Lines starting with `#` are comments. A
name that is not one of the built-in PTMs of `res/base_data/ptm_base.xml`
(such as `Acetyl`, `Phospho`, `Oxidation`) may be used on one line only.
A line whose name is a built-in PTM uses the built-in entry and its mass
(a WARNING is printed when the line's mass differs by more than 0.01 Da);
the built-in PTMs are every Unimod modification under its Unimod (PSI-MS)
name; the names of earlier versions, such as `Carbamidomethylation` (now
`Carbamidomethyl`), are no longer built in (see
[Proteoform annotation](proteoform.md), section 3).

In the variable PTM file (`-b`), `N-term` is the N-terminus of the
proteoform, the alpha-amine of its first residue: a truncation N-terminal
modification. It is searched only on a proteoform without a protein
N-terminal form (`-n`) that starts at the third residue of the protein or
later (the first two are those of the protein N-terminal forms), and only
when the proteoform's first residue is one of the line's residues. The
first residue's side chain can still carry an `any` modification. At most
one such modification is on a proteoform, and it counts toward `-S`.
Examples:

```text
Acetyl,42.010565,*,N-term,1
Gln->pyro-Glu,-17.026549,Q,N-term,28
```

With `Acetyl,42.010565,K,any,1` and `Acetyl,42.010565,*,N-term,1` in one
file, a proteoform starting with K can be acetylated on its N-terminus, on
the side chain of K, or both (`[Acetyl]-(K)[Acetyl]...`). One acetylation
of the first K gives the same fragment masses on the N-terminus and on the
side chain; TopPIC then reports the N-terminal one. In the tables a
truncation N-terminal modification is the prefix `[Name]-` of the
`Proteoform` string and `Name:[N-term]`, first, in `Variable PTMs`; `#
Variable PTMs` counts it (see [Proteoform annotation](proteoform.md)).

The protein N-terminal forms themselves are set with `-n`. The variable PTM
stage does not use `C-term` modifications.

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

# Three fractions searched, post-matched and combined into combined_ms2_toppic_*
toppic -c combined proteins.fasta frac1_ms2.msalign frac2_ms2.msalign frac3_ms2.msalign
```

The identifications of the first example are in
`sample_post_ms2_toppic_prsm_single.tsv` (all PrSMs),
`sample_post_ms2_toppic_proteoform_single.tsv` (one per proteoform) and
`sample_post_ms2_toppic_protein_single.tsv` (one per protein).
