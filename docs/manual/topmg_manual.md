# TopMG manual

TopMG (Top-down mass spectrometry-based proteoform identification using Mass
Graphs) identifies highly modified proteoforms, such as histone proteoforms,
by searching deconvoluted top-down MS/MS spectra against a protein sequence
database. Candidate proteoforms carrying several variable modifications are
represented compactly as proteoform graphs, spectra as mass graphs, and the
two are aligned. The statistical significance of each identification is
estimated with a Markov chain Monte Carlo method (TopMCMC).

```sh
topmg [options] -i modification-file database-file spectrum-file ...
```

`topmg -h` prints the option list. The GUI tool `topmg_gui` offers the same
function. Avoid spaces and quotation marks in paths (the directory of the
`topmg` executable must not contain spaces), and the database and spectrum
paths must be at most 200 characters long.

## 1. Input

**Protein database.** A FASTA file (`.fasta`, `.fa`, `.FASTA` or `.FA`)
with unique protein accessions.

**Modification file (required).** TopMG has no built-in list of
modifications: the variable modifications used to build the proteoform
graphs must be given with `--mod-file-name` (`-i`). The file has the same
format as TopPIC's modification files, one modification per line with five
comma-separated fields:

```text
# Name,Mass,Residues,Position,UnimodID
Acetyl,42.010565,K,any,1
Methyl,14.015650,KR,any,34
Phospho,79.966331,STY,any,21
```

`Residues` lists the residues that can carry the modification or `*` for
any; `Position` is `any`, `N-term` or `C-term` (`*` and `any` cannot be
combined); `UnimodID` is `-1` when unknown. Lines starting with `#` are
comments. The file is checked before the search starts.

**Spectrum files.** One or more `_ms2.msalign` files written by TopFD or
TopDIA, searched one after another with the same parameters.

**Companion file from TopFD.** For `sample_ms2.msalign`, TopMG reads
`sample_ms2.feature` from the same directory to cluster the identifications
of one proteoform and to report abundances; a missing feature file is an
error, and PrSMs whose spectrum has no feature are discarded at the
clustering step. `--no-topfd-feature` (`-x`) runs without it, clustering
proteoforms by protein and precursor mass instead. If `sample.sqlite`
exists, the identifications are also written into it; it is not required.

**Database indexes.** `topindex` is optional. TopMG prepares the database in
`<database>_idx/` next to the FASTA file and builds its filtering indexes in
memory unless `topindex` stored them with the same `-f`, `-n`, `-e` and
`-d` settings (see the [TopIndex manual](topindex_manual.md)).

## 2. What TopMG does

For each spectrum file TopMG prints its parameters and runs these steps.
Intermediate files (`sample_ms2.topmg_*`, `sample_ms2.msalign_*` and two
`_cutoff_xml` directories) are deleted at the end unless `--keep-temp-files`
is given; `sample_ms2.topmg_raw_prsm` is always kept.

1. **ASF-One PTM filtering.** Approximate spectrum-based filtering selects
   candidate proteins allowing one modification from the modification file,
   separately for complete, N-terminal, C-terminal and internal proteoforms.
2. **ASF-Diagonal PTM filtering** (only with `--use-asf-diagonal`): a second,
   diagonal-based filter that adds more candidates.
3. **Combining filtering results.**
4. **Graph alignment.** Each candidate protein is expanded into a proteoform
   graph whose edges span up to `--proteo-graph-gap` residues and carry up
   to `--var-ptm-in-gap` variable modifications, at most `--var-ptm` per
   proteoform. The spectrum's mass graph is aligned to it by dynamic
   programming, allowing `--num-shift` unexpected mass shifts.
5. **Graph alignment post-processing.** The aligned mass differences are
   converted into concrete modification assignments and mass shifts.
6. **E-value computation using MCMC.** TopMCMC estimates a p-value for each
   PrSM by sampling random proteoforms and converts it into an E-value.
   PrSMs with fewer than four matched fragments get the maximum E-value.
7. **Top PrSM selection.** The best PrSM of each spectrum is kept.
8. **Recounting matched masses and fragments** against the complete spectra
   (the search uses only masses with an EnvCNN score of at least
   `--filter-by-env-cnn`).
9. **Proteoform and protein clustering.** PrSMs sharing a TopFD feature form
   one proteoform; proteoforms of the same protein whose precursor masses
   agree within `--proteoform-error-tolerance` are merged; proteoforms are
   then grouped by protein.
10. **FDR computation** (only with `--decoy`).
11. **Filtering and output at three levels.** PrSMs by the spectrum-level
    cutoff, the best PrSM per proteoform by the proteoform-level cutoff, and
    the best PrSM per protein by the protein-level cutoff.
12. **Writing identifications to the SQLite database**, if `sample.sqlite`
    exists.

TopMG has no post mass matching step, so its result files keep the `_ms2`
name.

## 3. Output files

For an input `sample_ms2.msalign`:

| File | Contents |
|---|---|
| `sample_ms2_topmg_prsm.tsv`, `..._prsm_single.tsv` | All PrSMs that pass the spectrum-level cutoff. |
| `sample_ms2_topmg_proteoform.tsv`, `..._proteoform_single.tsv` | One PrSM per identified proteoform, passing the proteoform-level cutoff. |
| `sample_ms2_topmg_protein.tsv`, `..._protein_single.tsv` | One PrSM per identified protein, passing the protein-level cutoff. |
| `sample_ms2_topmg_prsm.xml`, `..._proteoform.xml`, `..._protein.xml` | The same three sets in XML, read by TopDiff and the visualization tools. |
| `sample.sqlite` | Updated with the identification tables when it exists. |
| `sample_ms2.topmg_raw_prsm` | The PrSMs before recounting and cutoffs. Always kept, even without `--keep-temp-files`; a later combined run (`-c`) reads it. |

The tables have the same layout and columns as TopPIC's (see the
[TopPIC manual](toppic_manual.md), section 3): a parameter block, the
numbers of identified PrSMs, proteoforms and proteins, then one line per PrSM
in the `_single` files, with extra lines for other proteins containing the
proteoform in the others. The `MIScore` column is always `-`, since TopMG
has no PTM characterization step; the `Variable PTMs` columns hold the
modifications assigned from the modification file. The annotation of the
`Proteoform` column (parentheses, brackets, named modifications and mass
shifts) is explained in [Proteoform annotation](proteoform.md).

## 4. Options

Database and modifications:

| Option | Default | Meaning |
|---|---|---|
| `-i`, `--mod-file-name <file>` | required | The modification file with the common PTMs used to build proteoform graphs (section 1). |
| `-f`, `--fixed-mod <C57\|C58\|file>` | none | Fixed modifications: `C57` carbamidomethylation of cysteine, `C58` carboxymethylation, or a modification file. |
| `-n`, `--n-terminal-form <list>` | NONE,NME,NME_ACETYLATION,M_ACETYLATION | Allowed protein N-terminal forms. |
| `-d`, `--decoy` | off | Also search a shuffled decoy database, to estimate false discovery rates. Required for the `FDR` cutoff types. |

Proteoform graphs and mass shifts:

| Option | Default | Meaning |
|---|---|---|
| `-P`, `--var-ptm <int>` | 5 | Maximum number of variable modifications in a proteoform. |
| `-j`, `--proteo-graph-gap <int>` | 40 | Maximum number of residues spanned by one edge of a proteoform graph. Larger values consider more combinations at the cost of memory and time. |
| `-G`, `--var-ptm-in-gap <int>` | 5 | Maximum number of variable modifications on one edge (never more than `--var-ptm`). |
| `-s`, `--num-shift <0\|1\|2>` | 0 | Maximum number of unexpected mass shifts in a proteoform, in addition to the variable modifications. |
| `-M`, `--max-shift <number>` | 500 | Maximum absolute mass of an unexpected shift (Da). |
| `-w`, `--whole-protein-only` | off | Report only proteoforms covering whole proteins, not truncated ones. |
| `-D`, `--use-asf-diagonal` | off | Add the ASF-Diagonal filtering step (more candidates, more sensitive, slower). |

Spectra and tolerances:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|HCD\|UVPD\|FILE>` | FILE | Fragmentation method; `FILE` takes it from each spectrum. |
| `-e`, `--mass-error-tolerance <int>` | 10 | Error tolerance of precursor and fragment masses (ppm). |
| `-F`, `--filter-by-env-cnn <0..1>` | 0.2 | Masses with an EnvCNN score below this are ignored during filtering, search and E-value computation. |
| `-p`, `--proteoform-error-tolerance <number>` | 1.2 | Precursor mass tolerance (Da) for grouping PrSMs into proteoforms. |
| `-x`, `--no-topfd-feature` | off | Run without the TopFD feature file. |

Cutoffs and output:

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--spectrum-cutoff-type <EVALUE\|FDR>`, `-v`, `--spectrum-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for PrSMs. |
| `-T`, `--proteoform-cutoff-type <EVALUE\|FDR>`, `-V`, `--proteoform-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteoforms. |
| `-y`, `--protein-cutoff-type <EVALUE\|FDR>`, `-Y`, `--protein-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteins; with `FDR` a protein is represented by its best proteoform. |
| `-c`, `--combined-file-name <name>` | none | Combine fractions: after the files are searched, their spectra, features (unless `-x`) and raw PrSMs are merged into `<name>_ms2.msalign`, `<name>_ms2.feature` and `<name>_ms2.topmg_raw_prsm`, and the post-search steps are repeated on them, giving `<name>_ms2_topmg_*` results. A relative `<name>` is placed in the directory of the first spectrum file and must differ from the input file names. No SQLite database is written for the combined results. |
| `-u`, `--thread-number <int>` | 1 | Number of threads. The filtering steps may use fewer, depending on the available memory. |
| `-k`, `--keep-temp-files` | off | Keep the intermediate files. |
| `-K`, `--keep-decoy-ids` | off | Keep decoy identifications in the result tables. |

## 5. Examples

```sh
# Histone sample: common histone PTMs, up to 5 per proteoform, 1% spectrum-level FDR
topmg -i histone_mods.txt -d -t FDR -v 0.01 -u 8 proteins.fasta sample_ms2.msalign

# Also allow one unexpected mass shift and the diagonal filter
topmg -i mods.txt -s 1 -D proteins.fasta sample_ms2.msalign

# Two fractions searched and combined
topmg -i mods.txt -c combined proteins.fasta frac1_ms2.msalign frac2_ms2.msalign
```

The identifications of the first example are in
`sample_ms2_topmg_prsm_single.tsv`, `sample_ms2_topmg_proteoform_single.tsv`
and `sample_ms2_topmg_protein_single.tsv`.
