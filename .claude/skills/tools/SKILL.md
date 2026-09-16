---
name: tools
description: Run the TopPIC Suite command-line tools (topfd, topdia, topindex, toppic, topmg, topdiff) on top-down mass spectrometry data and read their outputs. Use when asked to deconvolute mzML/mzXML or a peak list, search spectra against a FASTA database, build database indexes, compare proteoform abundances across samples, choose options, chain the tools into a pipeline, or interpret msalign/feature/result tables and the proteoform annotation.
---

# Running the TopPIC Suite tools

Six command-line tools, built into `<repo>/bin/` (see the `compile` skill),
each with a `*_gui` counterpart that shells out to it. Every tool prints its
option list with `-h`. The full manuals are in `docs/manual/`; this skill is
the operational summary. Consult the manual before quoting a default or an
option the summary does not mention.

| Tool | Does | Manual |
|---|---|---|
| `topfd` | Deconvolutes mzML/mzXML (or one peak list) into `msalign` masses and proteoform features | `docs/manual/topfd_manual.md` |
| `topdia` | Deconvolutes top-down DIA data and builds pseudo-MS/MS spectra | `docs/manual/topdia_manual.md` |
| `topindex` | Stores the protein-filtering indexes of a FASTA database on disk | `docs/manual/topindex_manual.md` |
| `toppic` | Identifies proteoforms with up to two unexpected shifts, characterizes them | `docs/manual/toppic_manual.md` |
| `topmg` | Identifies highly modified proteoforms with mass graphs | `docs/manual/topmg_manual.md` |
| `topdiff` | Compares proteoform abundances across samples | `docs/manual/topdiff_manual.md` |

The proteoform string in the result tables is explained in
`docs/manual/proteoform.md`.

## Pipeline

```text
raw -> msconvert (centroided mzML) -> topfd -> sample_ms2.msalign (+ _ms2.feature, .sqlite)
                                              -> toppic | topmg -> *_prsm/proteoform/protein tables
                                              -> topdiff across samples
topindex (optional, once per database + settings) speeds up toppic/topmg
```

Rules that hold across the tools:

- **All output goes next to the input file**, named from the input path
  minus its extension. Nothing has an output-directory option (only
  `topdiff -o` names its table, still in the first sample's directory).
- **No spaces or quotation marks in paths.** The executable's directory
  must not contain a space (hard error); input paths with spaces are simply
  not found. Database and spectrum paths: at most 200 characters.
- **Companion files are found by name** in the spectrum file's directory:
  `sample_ms2.msalign` implies `sample_ms2.feature` and `sample.sqlite`.
  Do not rename or move one without the others.
- **FASTA** extensions `.fasta`, `.fa`, `.FASTA`, `.FA`; accessions must
  be unique. TopPIC/TopMG/TopIndex create `<database>_idx/` next to it.
- Several input files on one command line are processed one after another
  with the same parameters (TopMG sorts them alphabetically).
- Threads: `-u N` on every tool except topdiff; N must not exceed the
  hardware threads, and the search tools may use fewer for filtering.

## topfd

```sh
topfd [options] sample.mzML ...        # centroided .mzML/.mzml/.mzXML/.mzxml
topfd -T [options] peaks.txt           # one MS/MS spectrum as "m/z intensity" lines
```

Convert vendor files with `msconvert --mzML --filter "peakPicking true 1-"`.
Outputs for `sample.mzML`: `sample_ms1.msalign`, `sample_ms2.msalign`
(the search input), `sample_ms1.feature`, `sample_ms2.feature`,
`sample_feature.xml`, `sample.sqlite`. FAIMS files give one set per
voltage, `sample_<voltage>_*`.

Key options (defaults): `-c` max charge 30, `-m` max mass 50000,
`-e` m/z error 0.02, `-a` activation FILE, `-u` threads 1,
`-r`/`-s` MS1/MS2 S/N ratio 3/1, `-t` ECScore cutoff 0.1, `-b` min scans 1,
`-o` no MS1 scans (no features, precursor set to the `-m`/`-c` limits),
`-N` no SQLite, `-D` also store raw MS1 peaks for 3D views.

- **Do not use `-N` if TopPIC will search the spectra**: TopPIC's default
  post mass matching needs `sample.sqlite` and fails without it (or run
  `toppic -E`).
- `-T` mode: only the MS/MS options apply; the peak list need not be
  sorted; it writes `peaks_ms2.msalign`, `peaks_ms2.env` (one line per
  matched peak, `PEAK_IDX` is the index after sorting by m/z) and
  `peaks.sqlite`. Useful for inspecting one spectrum's deconvolution; the
  reference test is the horse myoglobin list
  `~/code/mms/data/myoglobin_peaks.txt` with `topfd -T -g`.

## topdia

```sh
topdia [options] sample.mzML ...
toppic -E ... proteins.fasta sample_ms2.msalign    # search the pseudo spectra
```

Needs MS1 and DIA MS/MS scans with recorded isolation windows. Writes the
pseudo spectra as `sample_ms2.msalign` plus `sample_ms2.feature`, the MS1
files, `sample.sqlite` (always) and intermediates (`_ms2_raw.msalign`,
`*.csv`) that are left in place. Options are TopFD's with different
defaults (`-w` 4.0, `-t` 0, `-b` 2) plus fragment-feature and pseudo-spectrum
options: `-T`/`-B` MS2 ECScore cutoff and min scans (0, 1), `-p`/`-P`
intensity-correlation cutoffs (0.5), `-v` pseudo score cutoff 0.55,
`-V` fragments always added 25. No `-N`, `-D` or `-T`-peak-list mode.
**Search pseudo spectra with `toppic -E`** (post mass matching assumes one
scan per spectrum). Never pass `-o`: the run aborts at the pseudo-spectrum
step.

## topindex

```sh
topindex [-d] [-f C57|C58|file] [-e ppm] [-n forms] [-u N] proteins.fasta
```

Optional; without it the search tools build the same indexes in memory
each run. Writes `<database>_idx/` with the prepared database (`_standard`,
`_target` or `_target_decoy`, block files, `.fai`) and three index families.
**The searches use the index only when their `-f`, `-n`, `-e` and `-d`
match** the ones given to topindex; otherwise they fall back silently.
Index files are large (roughly ten thousand times the FASTA size, double
with `-d`); delete `_idx/` to reclaim space.

## toppic

```sh
toppic [options] proteins.fasta sample_ms2.msalign ...
```

Needs `sample_ms2.feature` (unless `-x`) and `sample.sqlite` (unless `-E`)
next to the spectrum file; a missing database is reported before the
search starts. Steps: zero-shift search, variable-PTM search (`-b`),
one-shift (`-s` >= 1) and two-shift (`-s 2`) searches, E-values, top PrSM,
recount, **post mass matching** (default), clustering, PTM characterization
(`-B`), FDR (`-d`), three-level cutoffs, SQLite tables.

Key options (defaults): `-f` fixed mod none (`C57`, `C58` or file), `-n`
N-terminal forms all four, `-R` proteoform types all four, `-d` decoy,
`-s` unexpected shifts 1 (0..2), `-m`/`-M` shift range -500/500 Da,
`-b` variable-PTM file with `-S` max 3, `-B` common-modification file for
characterization with `-H` MIScore threshold 0.15, `-e` tolerance 10 ppm,
`-F` EnvCNN filter 0.2, `-p` cluster tolerance 1.2 Da, `-r` combined
spectra 1, `-x` no feature file, `-E` disable post matching, `-I` post
min peaks 1, cutoffs `-t/-v`, `-T/-V`, `-y/-Y` (EVALUE 0.01 each; `FDR`
needs `-d`), `-c <name>` combine fractions, `-u` threads, `-k` keep temp
files, `-K` keep decoy ids.

Result naming: with post mass matching (default) everything is named
`sample_post_ms2_toppic_{prsm,proteoform,protein}[_single].tsv` and `.xml`,
plus `sample_post_ms2.msalign`; with `-E` the names are `sample_ms2_toppic_*`.
`_single` tables have one line per PrSM; the others add a line per extra
protein containing the proteoform. `sample_ms2.toppic_raw_prsm` is always
kept. Combined runs (`-c name`) always produce `name_ms2_toppic_*` and no
SQLite output; `name` must differ from the input base names.

Modification file (for `-f`, `-b`, `-B`): one `Name,Mass,Residues,Position,UnimodID`
per line, `#` comments, `Residues` may be `*` (not together with `Position`
`any`), `Position` is `any`/`N-term`/`C-term`. A name that matches a
built-in abbreviation in `res/base_data/ptm_base.xml` uses the built-in
mass, not the file's.

## topmg

```sh
topmg [options] -i mods.txt proteins.fasta sample_ms2.msalign ...
```

`-i` is required (the variable modifications for the proteoform graphs;
same file format as TopPIC). Needs `sample_ms2.feature` unless `-x`;
`sample.sqlite` is optional. No post mass matching, so results are always
`sample_ms2_topmg_*`; `MIScore` is always `-`. Key options (defaults):
`-P` max variable PTMs 5, `-j` graph gap 40, `-G` PTMs per gap 5, `-s`
unexpected shifts 0, `-M` max shift 500, `-D` add the ASF-Diagonal filter,
`-w` whole proteins only, plus the same `-f -n -d -a -e -F -p -x`, cutoff,
`-c`, `-u`, `-k`, `-K` options as TopPIC. PrSMs with fewer than four
matched fragments get the maximum E-value.

## topdiff

```sh
topdiff [-e 1.2] [-t toppic|topmg] [-o sample_diff.tsv] s1_ms2.msalign s2_ms2.msalign ...
```

At least two samples. For each spectrum file it reads
`<base>_post_ms2_<tool>_proteoform.xml`, else `<base>_ms2_<tool>_proteoform.xml`
(so either `_ms2.msalign` or `_post_ms2.msalign` may be given), prints the
name, and warns "does not contain any PrSM identifications" when neither
exists. Samples must have been searched with the feature file (not `-x`)
for abundances. Output: one TSV next to the first sample with the
proteoform, precursor mass, `# matched samples`, then five columns per
sample (abundance, spectrum id, RT begin/end in seconds, normalized apex).

## Reading the outputs

- `msalign`: `BEGIN IONS`/`END IONS` blocks, header `KEY=value` lines
  (`PRECURSOR_MASS`, `PRECURSOR_CHARGE`, `PRECURSOR_FEATURE_ID`, ...), then
  one mass per line: monoisotopic mass, intensity, charge, score. `#` lines
  at the top record the parameters. Several precursors are `:`-separated.
- Result TSV tables start with a parameter block and identification
  counts, then the header. Column meanings: TopPIC manual section 3.
  Q-value columns are `-` without `-d`; feature columns are `-` with `-x`.
- `Proteoform` column: `(X)[Name]` named modification on a residue,
  `(range)[+123.4567]` unexpected shift somewhere in the range, nested
  parentheses when a shift spans fixed-modified residues, `[Acetyl]-`
  N-terminal acetylation. Flanking residues are separate columns (`-` at a
  protein terminus). Full description: `docs/manual/proteoform.md`.
- Shifts near 0 or ±1 Da are precursor or isotope errors, not modifications.

## When the code changes

Options, defaults, output names and columns come from
`src/console/<tool>_argument.cpp`, the `*_process.cpp` drivers and the
writers in `src/prsm`, `src/ms/spec` and `src/merge`. Nothing checks the
manuals against them: after changing any of those, update the tool's
manual (and this skill if the summary is affected) in the same commit.
