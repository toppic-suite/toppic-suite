# TopDiff notes

Facts to keep straight when changing topdiff (`src/merge`,
`console/topdiff_*`) and topdiff_gui.
Moved out of the repository's `CLAUDE.md` on 2026-10-12; the code-wide
conventions, the dependencies and the source layout are there.

## TopDiff (`merge/feature_sample_merge`, `console/topdiff_process`)

TopDiff reads one `<base>_post_ms2_<tool>_proteoform.xml` (or `_ms2_`)
per sample into `FeaturePrsm`s (protein, residue range, adjusted precursor
mass, proteoform cluster id and abundance, feature min/max/apex retention
times), normalizes each sample's times by its latest feature end, aligns
samples 2..n onto sample 1 (`sampleAlignTime`: a dynamic programme over
the 1,000 most abundant proteoforms of each, sorted by apex time; a cell
matches when a feature of the other sample within 0.1 normalized time has
a precursor mass within `--error-tolerance`, scored by log10 abundance;
the backtracked matches are the piecewise-linear warp `setAlignTime`
applies), then builds the table (`getFeatureTable`): all proteoforms of
all samples by decreasing abundance, each not yet in the table starts a
row and takes, in every other sample, the most abundant unused proteoform
of the same protein within `--error-tolerance` and 0.3 aligned apex time
(`findMatchFeature`). Facts to keep straight:

- **The matching is by protein accession, precursor mass and aligned
  time only**; the residue range and the proteoform string are not
  compared. On 480_F3_02 result sets of the same TopFD features
  (toppic_test experiment 65) that matches every proteoform of a sample
  with itself and the same spectrum (1,765 of 1,765), and across search
  settings 1,780 of 1,799 matched rows join the same MS/MS scan; 4 rows
  join different features of one protein within 1.2 Da.
- **The hot spots were exact rewrites** (experiment 65, byte-identical
  tables): the alignment's per-cell `findMatch` scan became per-row sorted
  match-position lists with a binary search (`matchPositions` /
  `findMatch`; three quarters of the run before), the matrix is a flat
  `std::vector<Cell>` with the cell matches kept for the backtracking,
  `featureUsed` is a set of (sample, cluster id), and `findMatchFeature`
  looks only at the protein's features through a per-sample index in the
  intensity order of the scan it replaced (the dead `prot == ""` branch is
  kept, with the strict tolerances, so the result is the same). 32 copies
  of one result set: 14.8 s -> 2.4 s; two samples: 0.56 s -> 0.17 s, most
  of it now the XML reading. `initCells` used to leave cell (0, b_size)
  uninitialized (fixed; no output change on that data).
- **Result files of other versions are refused**, not parsed: every
  `<prsm>` starts with a `<version>` element (`Prsm::toXmlElement`,
  `Version::getVersion()`), and `PrsmReader::readOnePrsmStr` exits with
  "was written by version X of TopPIC Suite and cannot be read by this
  version" when it differs from the running version, or "by an earlier
  version" when it is absent (1.9 and earlier wrote `match_peak_num` /
  `match_fragment_num` and no `rt`; `PrsmStr` threw `std::invalid_argument`
  on them). Every tool that reads PrSM files gets this check; a result
  file needs the search rerun after a version change.
- The console prints `Aligning retention times of sample <k>` and
  `Matching proteoforms across samples` as started/finished steps (the
  manual lists them); the table is written next to the **first** input
  file, so give inputs from one directory or expect the output there.
- Not done: matching an identified proteoform to an unidentified TopFD
  feature of another sample (the commented-out "Match without MS/MS ID"
  path) - the `_ms1.feature` files are not read.
