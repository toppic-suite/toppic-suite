# TopIndex manual

TopIndex builds the index files that TopPIC and TopMG use to filter protein
candidates before aligning spectra. Neither search tool requires it: when
the index files are absent they build the indexes in memory for each run.
TopIndex stores them on disk once, so repeated searches against the same
database skip that work.

```sh
topindex [options] database-file
```

`topindex -h` prints the option list. The GUI tool `topindex_gui` offers the
same function. Paths must not contain spaces, and the database path must be
shorter than 200 characters.

## 1. Input

A protein database in FASTA format (`.fasta` or `.fa`). Every protein needs
a unique accession: if two share one, TopIndex reports the duplicate,
removes the index directory it started, and stops.

## 2. What TopIndex does

1. Creates the directory `<database>_idx/` next to the FASTA file (for
   `proteins.fasta`, the directory `proteins.fasta_idx/`).
2. Prepares the database in that directory: a cleaned copy
   (`proteins.fasta_standard`), the search database `proteins.fasta_target`
   or, with `--decoy`, `proteins.fasta_target_decoy` (with shuffled decoy
   proteins whose accessions start with `DECOY_`), its block files and block
   index, and its `.fai` index for random access. These are the same files
   TopPIC and TopMG create when they run without an index.
3. Writes three families of index files, one file per database block, each
   processed in a thread of its own:
   - `Generating non shift index files`: the zero-shift indexes
     `zero_ptm_term_index`, `zero_ptm_diag_index`, `zero_ptm_rev_term_index`
     and `zero_ptm_rev_diag_index`, used by TopPIC's zero-shift and
     variable-PTM filters.
   - `Generating one shift index files`: the one-shift indexes
     `one_ptm_term_index`, `one_ptm_diag_index`, `one_ptm_rev_term_index` and
     `one_ptm_rev_diag_index`, used by TopPIC's one-shift filter and TopMG's
     ASF-One PTM filter.
   - `Generating multiple shift index files`: `multi_ptm_index`, used by
     TopPIC's multiple-shift filter and TopMG's ASF-Diagonal filter.

Each file name carries the parameters the index depends on, in the form
`<family>_<fixed mod>_<N-terminal forms>_<error tolerance>_<no_decoy|decoy><block>`,
for example `zero_ptm_term_index_C57_N_NME_NMEA_MA_10_decoy0`. The N-terminal
forms are abbreviated `N`, `NME`, `NMEA` and `MA`; a fixed-modification file
appears by its file name.

## 3. Using the indexes

TopPIC and TopMG look for index files whose names match their own `-f`,
`-n`, `-e` and `-d` settings. **Run TopIndex with the same four settings
as the searches** that should use it; the thread number does not matter.
Indexes built with other settings are ignored silently, and the search
builds its own in memory. TopIndex skips database files that already exist,
so several parameter sets can be indexed into the same `_idx` directory.

The index files are large: roughly ten thousand times the size of the FASTA
file, spread over the three families, and about twice as much with
`--decoy`. A few-megabyte database can produce tens of gigabytes of index
files. Delete the `_idx` directory to reclaim the space; the searches then
fall back to in-memory indexes.

## 4. Options

| Option | Default | Meaning |
|---|---|---|
| `-f`, `--fixed-mod <C57\|C58\|file>` | none | Fixed modifications: `C57` carbamidomethylation of cysteine, `C58` carboxymethylation, or a modification file (see the [TopPIC manual](toppic_manual.md), section 5). |
| `-n`, `--n-terminal-form <list>` | NONE,NME,NME_ACETYLATION,M_ACETYLATION | Allowed protein N-terminal forms, comma-separated. |
| `-d`, `--decoy` | off | Also index a shuffled decoy database. |
| `-e`, `--mass-error-tolerance <int>` | 10 | Error tolerance of precursor and fragment masses (ppm). |
| `-u`, `--thread-number <int>` | 1 | Number of threads, one database block per thread. Each thread needs about 0.5 GiB of memory; TopIndex warns when the machine has too little for the requested number. |


## 5. Example

```sh
# Index for target-decoy searches with carbamidomethylated cysteines at 15 ppm
topindex -d -f C57 -e 15 -u 8 proteins.fasta

# Searches that use the index: same -d, -f and -e
toppic -d -f C57 -e 15 proteins.fasta sample_ms2.msalign
topmg  -d -f C57 -e 15 -i mods.txt proteins.fasta sample_ms2.msalign
```

The console shows the three phases with one progress line per block and ends
with `TopIndex - finished.`
