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
same function. Avoid spaces in paths (the directory of the `topindex`
executable must not contain spaces), and the database path must be at most
200 characters long.

## 1. Input

A protein database in FASTA format (`.fasta`, `.fa`, `.FASTA` or `.FA`).
If two records share an accession (the header text up to the first
space), TopIndex writes `<database>_dedup.fasta` next to the input,
keeping the first record of each accession, prints a notice, and builds
the indexes for that file in `<database>_dedup.fasta_idx/`. TopPIC and
TopMG do the same when given the original file, so they find the indexes.

## 2. What TopIndex does

1. Creates the directory `<database>_idx/` next to the FASTA file (for
   `proteins.fasta`, the directory `proteins.fasta_idx/`).
2. Prepares the database in that directory: a cleaned copy
   (`proteins.fasta_standard`), the search database `proteins.fasta_target`
   or, with `--decoy`, `proteins.fasta_target_decoy` (with shuffled decoy
   proteins whose accessions start with `DECOY_`), its block files and block
   index, and its `.fai` index for random access. These are the same files
   TopPIC and TopMG create when they run without an index.
3. Writes three families of index files, one file per index name and
   database block, in one pass over the blocks (`Generating index files`);
   blocks are processed in parallel on the `-u` threads:
   - the zero-shift indexes `zero_ptm_term_index`, `zero_ptm_diag_index`,
     `zero_ptm_rev_term_index` and `zero_ptm_rev_diag_index`, used by
     TopPIC's zero-shift and variable-PTM filters, cover proteoform masses
     up to 100 kDa;
   - the one-shift indexes `one_ptm_term_index`, `one_ptm_diag_index`,
     `one_ptm_rev_term_index` and `one_ptm_rev_diag_index`, used by TopPIC's
     one-shift filter and TopMG's ASF-One PTM filter, are the zero-shift
     indexes cut at 50 kDa, so each is written from the index just built
     rather than built a second time;
   - `multi_ptm_index`, used by TopPIC's multiple-shift filter and TopMG's
     ASF-Diagonal filter, is the same index as `one_ptm_diag_index`, so it is
     written as a hard link to that file (a copy on file systems without
     hard links); `du` counts the pair once.

Each file name carries the parameters the index depends on, in the form
`<family>_<fixed mod>_<N-terminal forms>_<error tolerance>_<no_decoy|decoy><block>`,
for example `zero_ptm_term_index_C57_N_NME_NMEA_MA_10_decoy0`. The N-terminal
forms are abbreviated `N`, `NME`, `NMEA` and `MA`; a fixed-modification file
appears by its file name without directory or extension (`-f mods/C57.txt`
gives `C57`).

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
| `-e`, `--frag-error-tolerance <int>` | 10 | Error tolerance of fragment masses (ppm). |
| `-E`, `--prec-error-tolerance <int>` | 10 | Error tolerance of precursor masses (ppm). |
| `-u`, `--thread-number <int>` | 1 | Number of threads; database blocks are processed in parallel. Must not exceed the number of hardware threads. Each thread needs about 0.5 GiB of memory; TopIndex warns (but still runs) when the machine has too little for the requested number. The build is limited by memory and disk bandwidth rather than by cores, so more threads than physical cores gain nothing (on an 8-core machine, `-u 8` was faster than `-u 16`). |


## 5. Example

```sh
# Index for target-decoy searches with carbamidomethylated cysteines at 15 ppm
topindex -d -f C57 -e 15 -u 8 proteins.fasta

# Searches that use the index: same -d, -f and -e
toppic -d -f C57 -e 15 proteins.fasta sample_ms2.msalign
topmg  -d -f C57 -e 15 -i mods.txt proteins.fasta sample_ms2.msalign
```

The console shows one progress line per block and ends with
`TopIndex - finished.`
