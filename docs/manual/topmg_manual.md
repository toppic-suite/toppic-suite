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

**Protein database.** A FASTA file (`.fasta`, `.fa`, `.FASTA` or `.FA`).
If two records share an accession (the header text up to the first
whitespace), the tool writes a copy that keeps the first record of each
accession next to the input (`proteins_dedup.fasta` for `proteins.fasta`),
prints a notice, and searches that file instead; its indexes go to
`proteins_dedup.fasta_idx/`.
A sequence may carry PTMs written in brackets after their residues, by
their names in `res/base_data/ptm_base.xml` (the Unimod names, for example
`M[Oxidation]`); the graph alignment treats them as fixed modifications of
those residues, and `-f` does not modify such a residue again.

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

An `N-term` line is a modification of the proteoform N-terminus (a
truncation N-terminal modification, as in TopPIC's variable PTM file): it
is searched only on a proteoform that starts at the third residue of the
protein or later and whose first residue is one of the line's residues,
and it counts toward `--var-ptm`. The protein N-terminal forms themselves
are set with `-n`. TopMG does not use `C-term` lines.

**Spectrum files.** One or more `_ms2.msalign` files written by TopFD or
TopDIA, searched one after another (sorted by their paths as given) with
the same parameters.

**Companion file from TopFD.** For `sample_ms2.msalign`, TopMG reads
`sample_ms2.feature` from the same directory to cluster the identifications
of one proteoform and to report abundances; a missing feature file is an
error, and PrSMs whose spectrum has no feature are discarded at the
clustering step. `--no-topfd-feature` (`-x`) runs without it, clustering
proteoforms by protein, precursor mass and sequence overlap instead (two
PrSMs share a proteoform only if their proteoform sequences overlap by at
least 70% of the shorter one). If `sample.sqlite`
exists, the identifications are also written into it (not with
`--no-sql`); it is not required.

**Database indexes.** `topindex` is optional. TopMG prepares the database in
`<database>_idx/` next to the FASTA file and builds its filtering indexes in
memory unless `topindex` stored them with the same `-f`, `-n`, `-e` and
`-d` settings (see the [TopIndex manual](topindex_manual.md)).

## 2. What TopMG does

For each spectrum file TopMG prints its parameters and runs these steps.
Intermediate files (`sample_ms2.topmg_*` and `sample_ms2.msalign_*`) are
deleted at the end unless `--keep-temp-files` is given;
`sample_ms2.topmg_raw_prsm` is always kept.

1. **ASF-One PTM filtering.** Approximate spectrum-based filtering selects
   candidate proteins allowing one modification from the modification file,
   separately for complete, N-terminal, C-terminal and internal proteoforms.
2. **ASF-Diagonal PTM filtering** (only with `--use-asf-diagonal`): a second,
   diagonal-based filter that adds more candidates.
3. **Combining filtering results.**
4. **Graph alignment.** Each candidate protein is expanded into a proteoform
   graph, whose vertices are the cleavage sites and whose edges are the
   residues, with or without a variable modification of the modification
   file (and the acetylated protein N-terminal forms of `-n`). The
   spectrum's mass graph is aligned to it by dynamic programming, allowing
   `--num-shift` unexpected mass shifts; a mass difference of the spectrum
   graph matches the mass between two vertices of the proteoform graph at
   most `--proteo-graph-gap` residues apart, with at most
   `--var-ptm-in-gap` variable modifications between them and at most
   `--var-ptm` in the proteoform, within 0.1 Da. The hidden option
   `--align-error-tolerance <number>` changes this tolerance (Da), which the
   post-processing also uses to merge modification masses; the MCMC E-value
   computation keeps its own 0.1 Da. `topmg_gui` does not offer it.
5. **Graph alignment post-processing.** The aligned mass differences are
   converted into the protein N-terminal form (or a truncation N-terminal
   modification), concrete modification assignments and mass shifts.
6. **E-value computation using MCMC.** TopMCMC estimates a p-value for each
   PrSM by sampling random proteoforms and converts it into an E-value.
   PrSMs with fewer than four matched fragments get the maximum E-value.
7. **Top PrSM selection.** The best PrSM of each spectrum is kept.
8. **Recounting matched masses and fragments** against the complete spectra
   (the search uses only masses with a confidence score, such as the EnvCNN
   score, of at least `--filter-by-score`).
9. **Proteoform and protein clustering.** PrSMs sharing a TopFD feature form
   one proteoform; proteoforms of the same protein whose precursor masses
   agree within `--proteoform-error-tolerance` are merged; proteoforms are
   then grouped by protein.
10. **FDR computation** (only with `--decoy`). Spectrum-, proteoform- and
    protein-level false discovery rates are estimated from the decoy hits;
    the decoy identifications are then removed unless `--keep-decoy-ids`
    is given.
11. **Filtering and output at three levels.** PrSMs by the spectrum-level
    cutoff, the best PrSM per proteoform by the proteoform-level cutoff, and
    the best PrSM per protein by the protein-level cutoff.
12. **Writing identifications to the SQLite database**, if `sample.sqlite`
    exists (skipped with `--no-sql`).

TopMG has no post mass matching step, so its result files keep the `_ms2`
name.

## 3. Output files

For an input `sample_ms2.msalign`:

| File | Contents |
|---|---|
| `sample_ms2_topmg_prsm.tsv`, `..._prsm_single.tsv` | All PrSMs that pass the spectrum-level cutoff. |
| `sample_ms2_topmg_proteoform.tsv`, `..._proteoform_single.tsv` | One PrSM per identified proteoform, passing the proteoform-level cutoff. |
| `sample_ms2_topmg_protein.tsv`, `..._protein_single.tsv` | One PrSM per identified protein, passing the protein-level cutoff. |
| `sample_ms2_topmg_prsm.xml`, `..._proteoform.xml`, `..._protein.xml` | The same three sets in XML, read by TopDiff and the visualization tools. Every PrSM records the TopPIC Suite version that wrote it, and TopDiff reads only the files of its own version. |
| `sample.sqlite` | TopFD's database, updated with the identification tables (`prsm`, `proteoform`, `protein`, `prsm_mass_shift`, `prsm_protein_match` and `fasta_seq`) when it exists. TopMG writes into TopFD's file itself, not into a copy as TopPIC does. With `--no-sql` it is not changed. |
| `sample_ms2.topmg_raw_prsm` | The PrSMs before recounting and cutoffs. Always kept, even without `--keep-temp-files`; with `-c` the files' raw PrSMs are merged from it. |

The tables have the same layout and columns as TopPIC's (see the
[TopPIC manual](toppic_manual.md), section 3): a parameter block, the
numbers of identified PrSMs, proteoforms and proteins, then one line per PrSM
in the `_single` files, with extra lines for other proteins containing the
proteoform in the others. The `MIScore` column is always `-` and the
`Alternative PTM explanation` column empty, since TopMG has no PTM
characterization step; the `Variable PTMs` columns hold the modifications
assigned from the modification file. The annotation of the `Proteoform`
column (parentheses, brackets, named modifications and mass shifts) is
explained in [Proteoform annotation](proteoform.md).

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
| `-j`, `--proteo-graph-gap <int>` | 40 | Maximum number of residues between two vertices of a proteoform graph whose mass difference is matched with the spectrum graph (section 2, step 4). Larger values consider more combinations at the cost of memory and time. |
| `-G`, `--var-ptm-in-gap <int>` | 5 | Maximum number of variable modifications between two such vertices (never more than `--var-ptm`). |
| `-s`, `--num-shift <0\|1\|2>` | 1 | Maximum number of unexpected mass shifts in a proteoform, in addition to the variable modifications. TopMG keeps its own number: `toppic --num-shift` does not change it. |
| `-M`, `--max-shift <number>` | 500 | Maximum absolute mass of an unexpected shift (Da). |
| `-w`, `--whole-protein-only` | off | Report only proteoforms covering whole proteins, not truncated ones, that is, allow the proteoform type COMPLETE alone. |
| `-D`, `--use-asf-diagonal` | off | Add the ASF-Diagonal filtering step (more candidates, more sensitive, slower). |

Spectra and tolerances:

| Option | Default | Meaning |
|---|---|---|
| `-a`, `--activation <CID\|ETD\|ETCID\|ETHCD\|HCD\|MPD\|UVPD\|FILE>` | FILE | Fragmentation method; `FILE` takes it from each spectrum. ETCID (ETciD) and ETHCD (EThcD) are searched with the ion types of ETD (c, z-dot). |
| `-e`, `--frag-error-tolerance <int>` | 10 | Error tolerance of fragment masses (ppm). |
| `-E`, `--prec-error-tolerance <int>` | 10 | Error tolerance of precursor masses (ppm). |
| `-F`, `--filter-by-score <0..1>` | 0.4 | Masses whose confidence score (the fourth column of the msalign file: the EnvCNN score in TopFD's files, the pseudo score in TopDIA's) is below this are ignored during filtering, search and E-value computation. In `topmg_gui`: "Confidence score cutoff". |
| `-p`, `--proteoform-error-tolerance <number>` | 1.2 | Precursor mass tolerance (Da) for grouping PrSMs into proteoforms. |
| `-x`, `--no-topfd-feature` | off | Run without the TopFD feature file. In `topmg_gui`: "No TopFD feature file". |

Cutoffs and output:

| Option | Default | Meaning |
|---|---|---|
| `-t`, `--spectrum-cutoff-type <EVALUE\|FDR>`, `-v`, `--spectrum-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for PrSMs. |
| `-T`, `--proteoform-cutoff-type <EVALUE\|FDR>`, `-V`, `--proteoform-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteoforms. |
| `-y`, `--protein-cutoff-type <EVALUE\|FDR>`, `-Y`, `--protein-cutoff-value <number>` | EVALUE, 0.01 | Cutoff for proteins; with `FDR` a protein is represented by its best proteoform. |
| `-c`, `--combined-file-name <name>` | none | Combine fractions: after the files are searched, their spectra (TopFD's `_ms1.msalign` and `_ms2.msalign` files), features (`_feature.xml` and `_ms2.feature`, unless `-x`) and raw PrSMs are merged into `<name>_ms1.msalign`, `<name>_ms2.msalign`, `<name>_feature.xml`, `<name>_ms1.feature`, `<name>_ms2.feature` and `<name>_ms2.topmg_raw_prsm`, and the post-search steps are repeated on them, giving `<name>_ms2_topmg_*` results. Every fraction's `_ms1.msalign` file (and `_feature.xml` file, unless `-x`) must exist. A relative `<name>` is placed in the directory of the first spectrum file (in the sorted order of section 1) and must differ from the input file names. No SQLite database is written for the combined results. |
| `-u`, `--thread-number <int>` | 1 | Number of threads; must not exceed the number of hardware threads. Each protein filtering step may use fewer: it fits its thread number to the memory available when it starts (about 1 GiB per thread for ASF-One PTM filtering, 0.5 GiB for ASF-Diagonal filtering). The run starts with a warning listing the filtering steps that will use fewer threads and how many, and a step prints its number again if the memory available when it starts gives a different one. The graph alignment also fits its thread number to the memory available when it starts (0.5 GiB per thread, times (`-s` + 1) x (`-P` + 1) / 8 when that is more than 1: 0.75 GiB with the defaults) and prints a warning when it uses fewer threads. All other steps use the full number. |
| `-k`, `--keep-temp-files` | off | Keep the intermediate files. |
| `-K`, `--keep-decoy-ids` | off | Keep decoy identifications in the result tables. |
| `-z`, `--no-sql` | off | Do not write the identification tables into TopFD's `sample.sqlite`: TopMG then writes or changes no SQLite file. The other result files are the same as without the option. |

## 5. Examples

```sh
# Histone sample: common histone PTMs, up to 5 per proteoform, 1% spectrum-level FDR
topmg -i histone_mods.txt -d -t FDR -v 0.01 -u 8 proteins.fasta sample_ms2.msalign

# Allow up to two unexpected mass shifts and add the diagonal filter
topmg -i mods.txt -s 2 -D proteins.fasta sample_ms2.msalign

# Two fractions searched and combined
topmg -i mods.txt -c combined proteins.fasta frac1_ms2.msalign frac2_ms2.msalign
```

The identifications of the first example are in
`sample_ms2_topmg_prsm_single.tsv`, `sample_ms2_topmg_proteoform_single.tsv`
and `sample_ms2_topmg_protein_single.tsv`.
