# toplib_build and toplib_search notes

Facts to keep straight when changing the C++ ports of TopLib
(`src/toplib`, `console/toplib_build*`, `console/toplib_search*`).
Moved out of the repository's `CLAUDE.md` on 2026-10-12; the code-wide
conventions, the dependencies and the source layout are there.

## toplib_build (`src/toplib`, `console/toplib_build`)

toplib_build is a C++ port of TopLib's `src/ms_library_building.py`
(`~/claude/toplib`, commit 7577b09): it clusters identified MS/MS spectra
and writes target and decoy representatives into the SQLite table
`target_decoy_spectra_representatives`. With the hidden
`--python-compatible` the port reproduces the Python script's library value
for value, and the tests that show it (toppic_test experiments 81 and 82,
reports 0081 and 0082: data/480_lib and 480_F3_02's TopPIC 1.9 and 2.0
output) compare every cell. By default it departs from the script in two
ways, and `--disable-different-proteoform-merge` adds a third (experiment
82):

- **A proteoform is the pair (proteoform id, protein accession)** in the
  step that keeps one cluster's identification per proteoform and
  precursor charge: TopPIC 1.9 gave PrSMs of different proteins sharing a
  precursor feature one proteoform id (43 of 1,818 ids on 480_F3_02; 2.0's
  cluster protein correction removed that), and the script, which groups
  by the id alone, dropped one protein's identification when two proteins
  shared an id.
- **The PrSM table's optional columns are also found under TopPIC 1.9's and
  2.0's names** (`idColumnMap`: TopLib's name, TopPIC 1.7's, first; then
  e.g. `Unexpected modifications`, `# Masses`, `# Matched masses` /
  `# Matched exp fragments`). The script looks for 1.7's names only, so
  from 1.9 and 2.0 tables it copied none of the eight renamed columns.
- **`--disable-different-proteoform-merge`** (visible, argument key
  `disableDifferentProteoformMerge`, default false; `build`'s
  `disable_different_proteoform_merge`): spectra of different proteoforms
  (the same pair) are never clustered together. In the spectral
  clustering a pair of identified spectra of different proteoforms gets
  the distance n * n (n the precursor cluster's size) instead of the
  cosine distance (`clusterGroup`'s `cannot_link`): every average linkage
  merge containing such a pair then averages at least 4, above the 0.7
  cutoff, so no flat cluster holds it, while unidentified spectra can still
  join any cluster. This was the default up to main 47b9fc57; the user made
  the script's behaviour the default (2026-10-04). By default every pair
  gets its cosine distance, so the clusters are the script's: they join
  different proteoforms (46 of 47 multi-member clusters on data/480_lib,
  34 of them different proteins), and as each representative keeps one
  identification, 50 of 4,630 proteoforms get none (8 of 1,908 on
  480_F3_02's TopPIC 2.0 output, 21 of 1,868 on its 1.9 output). With the
  option the library and cluster files equal main 47b9fc57's default byte
  for byte on experiment 82's inputs; without it the cluster files equal
  `--python-compatible`'s. The parameter block prints
  `Different proteoform merge: enabled|disabled`. Do not replace the
  option by a stricter cosine cutoff: that changes the clustering of every
  spectrum, not only of the mixed ones, and still does not guarantee the
  separation.

`--python-compatible` turns the first two off (it is what the comparisons
with the Python script use); it is independent of
`--disable-different-proteoform-merge`, so without that option it
reproduces the script. Facts to keep straight:

- **The Python choices are kept on purpose**, including odd ones: the
  top-50 peaks, precursor clusters by charge then complete linkage within
  2.2 Da, the 10 ppm cosine distance with an isotope shift (TopPIC's
  `mass_constant::getIsotopeMass()`, 1.00228 Da; the script's 1.00235 Da
  is not kept) and its
  early return when the sum passes 1, average linkage at 0.7, the
  `weighted_mass_sum_inte_cal` merge (a peak of the first spectrum takes
  every matching peak of the second; a second-spectrum peak may be merged
  into several), the precursor representative from the truncated integer
  mass differences, the lowest-E-value identification and the removal of
  all but the best cluster per (precursor charge, proteoform id), the decoy
  loops (redraw while negative, then while above the precursor mass, the
  second loop not rechecking the first), and pandas' column types. Each
  function's comment names the Python function it ports; keep that.
- **Exact arithmetic of the Python libraries** (`lib_numeric`,
  `lib_cluster`): SciPy's `nn_chain` linkage (with its tie rule, the stable
  sort by distance and `label`), `fcluster(criterion='distance')`'s
  traversal order (the cluster ids), numpy's pairwise sum, CPython's
  `sum()` (Neumaier, Python >= 3.12) for the PPIR, CPython's Mersenne
  Twister with `random.seed(int)` and `random.sample` (the decoys), pandas'
  `precise_xstrtod` for the PrSM table's numbers (not always correctly
  rounded: 395 of the data set's 5,360 E-values differ from `strtod` by an
  ulp), and `np.linalg.norm` as OpenBLAS's AVX-512 (SkylakeX) `ddot`
  computes it, with `std::fma`. The `src/toplib` sources are compiled with
  `-ffp-contract=off` (CMakeLists.txt) so that nothing else is fused.
- **Not reproduced:** numpy's `log2` (Intel SVML on AVX-512) differs from
  glibc's by an ulp for about 0.006% of values, so a few normalized
  intensities differ by an ulp; numpy's default sort is not stable and the
  port uses stable sorts, so equal values (two peaks of one mass with
  different charges, equal E-values) may come out in another order. The
  norm matches numpy only where numpy's OpenBLAS picks the SkylakeX kernel.
  Everything else was bit-identical in experiment 81.
- **Inputs:** a TopFD msalign file plus a PrSM table (the script's input;
  the header is the first line with a `Spectrum ID` column, so TopPIC's
  parameter block is skipped, which TopLib does with a separate script), or
  a library msalign file (v1.2.1 keys, identifications in the headers;
  `lib_msalign_reader::readLibraryMsalign` maps them onto the table's
  columns). The optional table columns are read under TopLib's names or
  TopPIC 1.9's and 2.0's (see above).
- The hidden option `--cluster-file <file>` writes every spectrum's cluster
  (for comparing the clustering with the Python pipeline's, which it does
  not write).

## toplib_search (`src/toplib/lib_search*`, `console/toplib_search`)

toplib_search is a C++ port of TopLib's `src/ms_library_query.py` (the
same commit): it searches the spectra of a TopFD msalign file against the
library table and writes `toplib_output/query_res.tsv`. It reproduces the
script's table byte for byte (toppic_test experiment 83, report 0083:
480_F3_02 against five libraries of data/480_lib and 480_F3_02, the
parameter variants and a stress query); there is no compatibility switch
because nothing in the search was changed. Facts to keep straight:

- **The script's choices are kept on purpose**: the library rows with a
  proteoform, targets and decoys alike; only the query's first precursor;
  the ppm test in ppm of the *library* mass with the isotope shift
  (`mass_constant::getIsotopeMass()`, 1.00228 Da, not the script's
  1.00235; Da mode has no shift); the cosine loop of
  `peak_50_cosine_mass_ppm` with the library spectrum first and its early
  return; the *first* candidate of the smallest distance (`np.argmin`,
  a NaN distance first) kept when Python's `min()` is at most 1; the
  cutoffs as `np.arange(0.1, 1.0, 0.001)` computes them (start + i *
  (the second value - start), not i * 0.001); the cutoff whose FDR is
  *closest* to 1% (it may be above 1%); the similarity filter of 0.3
  afterwards. The decoys share their target's precursor mass and rarely
  win (40 of 3,036 best matches on 480_F3_02), so the chosen distance
  cutoff is about 0.98 and the 0.3 similarity filter, not the FDR,
  decides which matches are written (none of the 40 decoys reaches it).
- **pandas' output**: the library columns keep the dtype `read_sql_query`
  gives them over the whole table (integers only -> int64; any REAL value,
  or an integer column with a NULL -> float64, written `924.0`; text or
  all NULL -> object), floats are written as Python's repr
  (`lib_numeric::pythonRepr`, `std::to_chars` shortest digits laid out as
  CPython's `format_float_short`), NULL as an empty cell, text quoted as
  `csv.QUOTE_MINIMAL` quotes.
- **Deliberate differences**, none of which changes a result on a single
  TopFD file: the library path is used as given (the script prepends
  `./TopLib/`); `-T` and `-c` are case-insensitive and checked (the
  script's `-T da` takes the 2.2 default but tests in ppm); tolerances
  must be nonnegative numbers; spectra of a repeated `SPECTRUM_ID` (a
  combined msalign file) are searched and reported one by one, where the
  script merges their masses and then fails ("All arrays must be of the
  same length").
- **Not reproduced**, as in toplib_build: numpy's `log2` (SVML) differs
  from glibc's by an ulp for a few values, which moves a distance by an
  ulp (1 query of 4,257 on 480_F3_02); pandas' descending `sort_values`
  and `np.argsort` are not stable, so a tie of intensities across the
  50-peak cut or two equal masses may be ordered differently (stable sorts
  here; none in 480_F3_02). The norm is `lib_numeric::numpyNorm` (OpenBLAS
  SkylakeX `ddot`), and `-ffp-contract=off` keeps the loops unfused as
  numba's are.
- The library search needs no index: candidates come from the library
  sorted by precursor mass with a window that covers the exact test, then
  in library order; a search of 4,257 spectra takes about 0.5 s.
- The hidden option `--match-file <file>` writes every query's match
  before the FDR step (`spectrum_id`, `file_name`, `lib_index`,
  `representative_id`, `flag`, `calculated_cosine_dist`), the table the
  script computes in `ms_query` and does not store.
