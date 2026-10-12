# TopPIC notes

Facts to keep straight when changing toppic (`console/toppic_*`,
`src/filter`, `src/search`, `src/stat`, `src/prsm`, `src/para`,
`src/ms/calib`) and toppic_gui: the filters and searches, the E-values,
post mass matching, the FDRs and the rescoring, the PTM localization,
the second precursor search, the m/z recalibration (also topcalib), the
index files (topindex writes them) and the PrSM tables and ProForma
strings that topmg writes too.
Moved out of the repository's `CLAUDE.md` on 2026-10-12; the code-wide
conventions, the dependencies and the source layout are there.

## Post mass matching (`prsm/prsm_post_mass_match`, toppic only)

After the search, the clustering, the cluster protein correction and the
match recount, `console/toppic_process` runs
`prsm_post_mass_match::process`, which matches the theoretical fragment
masses of each PrSM that no deconvoluted mass matched against the
**centroided MS/MS peaks of TopFD's SQLite database, read from toppic's copy
of it** (`<base>_toppic.sqlite`; MSPathFinderT's envelope test), scores the
added masses with EnvCNN, and writes
`<base>_post_ms2.msalign` plus PrSM/feature files under that name, which every
downstream step (PTM characterization, FDR, tables, SQL output) then reads
instead of the input. Facts to keep straight:

- **It is on by default.** `initArguments()` sets `postMassMatch = "true"`;
  `-D` / `--disable-post-match` sets it to false (there is no enabling flag
  any more), and `-I` / `--post-min-peak-num` (default 1) is the minimum
  number of observed isotopic peaks per accepted fragment.
- **toppic's database, `<base>_toppic.sqlite`** (branch `toppic_sqlite`,
  toppic_test experiment 156, the user's design; merged into main on
  2026-10-11 as 8b1f61ac): one database per
  spectrum file given to toppic, whatever the options, so that a TopFD and
  a TopPIC run leave two databases per file and toppic never writes into
  TopFD's. `TopPICProgress_multi_file` names it for each input file
  (`prsm_post_mass_match::toppicSqlFileName` of the given file, the base of
  TopFD's `<base>.sqlite`) and, before the file is searched, deletes an
  earlier one (`sql_db_copy::removeDb`, with its -journal / -wal / -shm
  files; "Deleted the earlier database ...") and writes a new one: a copy
  of TopFD's database with every table (`writeToppicDb`,
  `sql_db_copy::copy`, SQLite's backup API, toppic's tables of an older
  version dropped, closed in WAL mode as TopFD writes its own; "Copied the
  database ... (n s)"), or with
  `--mass-calib` the corrected copy that `TopPIC_mass_calib` writes after
  the calibration search (a plain copy when there are too few
  calibrants). The name reaches every step through the argument map
  (`toppicDbFileName`, and `topfdDbFileName` for the calibration), not
  through the name of the file a step reads: `prsm_post_mass_match::
  process` and `prsm_sql_output::write` take it as an argument, the second
  pass of `--second-precursor` gets it with its copy of the arguments, `-c`
  sets each fraction's for its report and "" for the combined result
  (which has no database: "No SQLite database for ..."). `-A` (TopFD's
  database holds monoisotopic masses) and `-D` use the same copy. `-C` (no
  search) neither deletes nor writes it: it rewrites toppic's tables in the
  earlier run's database, which holds that run's post-matched envelopes; a
  missing one gives the "not found" note. Without TopFD's database (TopFD
  `-N`, possible with `-D`) nothing is written and the PrSM tables are not
  written to a database. The envelopes that toppic versions before this
  change appended to TopFD's database cannot be told apart from TopFD's and
  are copied with them; TopFD again gives a clean database. TopFD writes
  its database in WAL mode, so the read-only connection of the copy leaves
  an empty `<base>.sqlite-wal` and a `<base>.sqlite-shm` next to it (a
  read-only connection cannot remove them; main up to b34fec7e opened it
  read-write, which removed them on close); TopFD's file keeps its md5.
  Checked on the merge against main b34fec7e on 300 spectra of 480_F3_02
  (toppic_test experiment 141's dev subset, the full TopFD database), -d
  -t FDR -T FDR -y FDR: main's result files, `<base>_toppic.sqlite` equal
  in all 15 tables to the database main wrote into, the same again in a
  second run. Main up to
  e5c15398 wrote into TopFD's database without `-A` and `--mass-calib`
  (every run appended its envelopes again) and into `<base>_mono.sqlite`,
  `<base>_calib.sqlite` or `<base>_mono_calib.sqlite` with them. topmg is
  unchanged: it writes its PrSM tables into TopFD's database
  (`sqlFileName` of its spectrum file, `_post_ms2` read as `_ms2`).
  Results (toppic_test experiment 156, report 0156; 8f6cae30 against main
  e5c15398 on experiment 154's subsets of 480_F3_02, setting d, and the
  full 480_F3_02 and Bruker runs): every branch folder holds TopFD's
  database, whose md5 no run changes, and `<base>_toppic.sqlite` only;
  the default, `-D`, `--mass-calib`, `-A`, `-A --mass-calib`, `-A -D`,
  `-c AB -A` (also with `--mass-calib`), `-C` and `--second-precursor -A
  --mass-calib` write main's result files, and `<base>_toppic.sqlite`
  equals, table for table, the database main wrote into (TopFD's after a
  default or `-D` run, `<base>_mono.sqlite`, `<base>_calib.sqlite`,
  `<base>_mono_calib.sqlite`; the same `mass_calib` model). Main's TopFD
  database gained the 3,209 post-matched envelopes of the subset again
  with every run; the branch's copy after a rerun equals the first run's.
  480_F3_02 3,581 / 2,040 / 445 and Bruker 806 / 83 / 72 PrSMs /
  proteoforms / proteins as main; the copy takes 0.4-2 s (875 MB-2.5 GB),
  the runs as long as main's. A read transaction held on the copy during a
  run made the first commit (4916bc98, the copy in rollback mode) stop with
  "database is locked"; with the WAL copy the run completes, as main does
  with the same reader on TopFD's database.
- **A missing database is an error, not a skip.** TopFD's database of a
  spectrum file is `prsm_post_mass_match::sqlFileName(sp_file)`
  (`<base>.sqlite` for `<base>_ms2.msalign`; use it, do not re-derive the
  name). It is checked in three places with the same message ("run TopFD
  without disabling its SQLite database output (without -N) ..., or run
  TopPIC with --disable-post-match"): `ToppicArgument::validateArguments`
  (TopFD's database, so a run fails before the search), `process()` itself
  (toppic's copy it is given; abort), and `toppic_gui`'s `checkError()`
  (TopFD's database; a warning dialog before toppic is launched). Keep all
  three in step.
- **It runs on the best PrSM of every spectrum, after the clustering, the
  cluster protein correction and the recount.** `TopPIC_postprocess`
  clusters and corrects the PrSMs (`toppic_raw_prsm` → `toppic_cluster_prot`),
  recounts them against the full spectra (`toppic_recount`), then
  post-matches `toppic_recount` (targets and decoys alike) into
  `<base>_post_ms2.toppic_recount`,
  and the PTM characterization (`-l`), FDR, cutoffs, tables and SQL output
  read the spectra with the added masses under the new name. Post matching
  changes neither E-values nor cluster ids, so the FDRs and the
  identifications equal those of a `-D` run; only the matched masses and
  fragments differ. (It used to run before the clustering, which let the
  cluster protein correction recompute E-values on the post-matched
  spectra: 42 E-values and 6 PrSMs at the 1% cutoff differed from `-D` on
  480_F3_02, toppic_test report 0062.) The PrSM XML carries the cluster
  ids (`proteo_cluster_id`, `prot_id`, `prot_cluster_id` in the proteoform
  element), so the read-modify-write in `process()` keeps them.
- **It reads each PrSM's peaks with one query on `ms2_peak` by `spec_id`**
  (`RawPeakReader`), and TopFD writes that table without an index, so the
  reader first runs `CREATE INDEX IF NOT EXISTS ms2_peak_spec_id_index`.
  Without the index every query scans the whole table: on 480_F3_02
  (4.5 M peaks, ~4,000 PrSMs) the step took 11 of the 18 minutes of a
  default run. Keep any new per-spectrum query on the SQLite file indexed.
- **Cluster protein correction runs right after the clustering, before
  the match recount** (`stat/tdgf/cluster_prot_corrector`; it lives in
  `stat` because it needs `EValueProcessor`): in every proteoform cluster
  the PrSMs whose protein differs from the protein of the cluster's
  lowest-E-value PrSM (the same proteoform mass matched to another
  protein) get that PrSM's proteoform (shared pointer; its mass becomes
  their adjusted precursor mass), and their matched fragment numbers and
  E-values (`EValueProcessor::compEvalues`) are recomputed on the spectra
  read as the search read them: `MsAlignReader` with the EnvCNN cutoff and
  `SpPara::getFragDedupPpo()` (the combined masses of `--frag-mass-dedup`);
  the recount that follows replaces the matched numbers with full-spectrum
  counts (a `TdgfMng` built from the same values as the search's, with
  `varPtmTypeNum` from the argument map, 0 if absent, and, with
  `--topmg pipeline` or `topmg-only`, max(`-s`, `--topmg-num-shift`)
  unexpected shifts, since a graph-alignment proteoform is corrected like
  any other). A corrected PrSM
  also takes the best PrSM's `<tool>` (toppic or topmg), so the per-search
  FDR counts it with that search. Only with `--topmg-use-mcmc` are the other
  proteins' PrSMs of a cluster whose best PrSM comes from the graph alignment
  dropped instead (their generating-function E-value would be counted
  against MCMC ones; `Cluster protein correction - removed N PrSMs`).
  Otherwise nothing is removed; a corrected PrSM that matches at most
  `comp_evalue_min_match_theo_frag_num_` fragments gets the maximum
  E-value and fails the cutoffs. The output is `toppic_cluster_prot`
  (deleted without `-k`); the recount reads it and writes `toppic_recount`
  with the matched numbers of every PrSM, corrected ones included, counted
  on the full spectra, which the post mass matching, the PTM
  characterization (`-l`) and the FDR use.
  Combined runs (`TopPIC_multi_postprocess`) do not run it. The earlier
  `prsm_cluster_prot_filter`, which dropped those PrSMs, is gone.
- **`<base>.toppic_record_prsm` is where `TopPIC_postprocess` ends and
  `TopPIC_report` starts.** It is a copy, made inline at the end of
  `TopPIC_postprocess`, of the PrSM file after the FDR computation
  (`toppic_cluster_fdr` with `-d`; `toppic_recount`, or
  `toppic_cluster_local` with `-l`, without): the best PrSM of every
  spectrum with the added masses, the cluster ids and the FDRs, before the
  cutoffs, **decoys included** (`prsm_fdr::process` and
  `prsm_fdr_groups::process` always write the decoys; their former
  keep-decoy parameter is gone). `cleanToppicDir` never deletes it. `TopPIC_report`
  (the decoy removal, the cutoffs, the tables and the SQL output) reads
  only the argument map and re-derives the file names from it, so the two
  stages can be run separately, and a combined run (`-c`) merges the
  fractions' record files and runs the report alone. Without `-K`,
  `prsm_decoy_filter::process` writes the PrSMs whose sequence name does
  not start with `DECOY_` to `toppic_prsm_target`, which the spectrum-level
  cutoff reads (with `-K` it reads `toppic_record_prsm` itself); the
  proteoform and protein cutoffs and their best-PrSM selections chain from
  `toppic_prsm_cutoff`, and `cleanToppicDir` deletes `toppic_prsm_target`
  and the cutoff files with the other intermediates.
- **The PrSM-only output steps take the spectrum base name, not the
  spectrum file.** `prsm_decoy_filter`, `prsm_cutoff_selector`,
  `prsm_form_filter`, `prsm_prot_filter` and `prsm_sql_output::write` only
  ever name `<base>.<ext>` files (the SQL writer takes its database as an
  argument: toppic's `toppicDbFileName`, topmg's `sqlFileName` of
  `<base>_ms2.msalign`, `_post_ms2` read as `_ms2`), so their name argument
  is
  `base_name`, the spectrum file name without its extension
  (`file_util::basename(sp_file)`); `TopPIC_report` and `TopMG_post`
  compute it once. Do not pass the spectrum file name to them. The steps
  that need the spectra themselves — `PrsmMatchTableWriter` through
  `prsm_reader_util::readPrsmsWithSpectra`, and `PrsmCoverage` — open the
  msalign named in `PrsmPara`, so the report still needs `<base>.msalign`
  on disk for the tables (the PrSM XML stores no peaks).
- **Result file names carry `_post_ms2`** whenever post matching ran
  (`postBaseName`: the spectrum file's own base, without `_ms2` and
  `_resid`: `<base>_mono_post_ms2` for the `_mono` copy of `-A`,
  `<base>_calib_post_ms2` for the corrected spectra of `--mass-calib`,
  `<base>_mono_calib_post_ms2` for both; with the `_resid` suffix of
  `--second-precursor`'s residual spectra kept: `<base>_resid_post_ms2`);
  the
  cleanup in `toppic_process` (`cleanToppicDir`) takes that name to find the
  intermediate files. With `-D` the names are the plain `<base>_ms2` ones.
- **Combined runs (`-c`) merge the spectra and record files, re-cluster,
  compute one FDR and report.** `TopPICProgress_multi_file` merges each
  fraction's post-matched spectra (`<frac>_post_ms2.msalign`, or
  `<frac>_ms2.msalign` with `-D`) into `<name>_ms2.msalign` with
  `msalign_frac_merge::mergeFractions`, and each fraction's
  `toppic_record_prsm` (`<frac>_post_ms2.toppic_record_prsm`, or
  `<frac>_ms2.…` with `-D`) into `<name>_ms2.toppic_raw_prsm`; `<frac>` is
  the spectrum file's base, with the `_mono` suffix of a `-A` run and the
  `_calib` suffix of `--mass-calib` (`prsm_post_mass_match::postBaseName`),
  and the MS1 spectra are TopFD's `<frac>_ms1.msalign` without them
  (`topfdBaseName`). The spectra are merged from the
  post-matched files because the report attaches them to the PrSMs, whose
  matches were counted on those spectra. `TopPIC_multi_postprocess` then
  clusters the merged PrSMs into proteoforms and proteins across the
  fractions with `prsm_simple_cluster` (by protein and precursor mass, not
  by feature; the feature files are not merged) and computes one FDR
  (decoys kept), writing
  `<name>_ms2.toppic_record_prsm`, and `TopPIC_report` runs on it. The
  re-clustering is required: `prsm_util::mergePrsmFiles` offsets the
  spectrum and feature ids per fraction but not the proteoform/protein
  cluster ids, so without it the fractions' clusters collide and the
  best-PrSM filters collapse them (79 combined proteoforms instead of 142
  on st_1 + st_2). There is no recount, post mass matching or PTM
  characterization on the merged set. The record file is always kept, so
  `-C` (combine-result-only) does not need `-k`. Combined results are
  always `<name>_ms2_toppic_*`.
- **With `-c` the files are reported after the combined result, with its
  proteoform ids.** `TopPICProgress` skips the file's `TopPIC_report`
  when `combinedOutputName` is set, and `TopPICProgress_multi_file` defers
  the file's `cleanToppicDir` (the report needs the post-matched spectra
  and their index, and the cleanup copies the report's cutoff files to the
  XML outputs). After the combined report,
  `prsm_util::updateProteoClusterIds` (the counterpart of
  `mergePrsmFiles`) gives every PrSM of each file's `toppic_record_prsm`
  the `proteo_cluster_id` of its combined PrSM (spectrum id
  `SpPara::getMaxSpecNumPerFile()` x file index + its own, same precursor
  id, protein, start and end; a PrSM without one, not expected, gets a new
  id above the combined ones with a warning) and rewrites the file; then
  each file gets `geneSpIndex` on its spectra (`<frac>_post_ms2.msalign`,
  or the input file with `-D`), `TopPIC_report` (start time = the combined
  run's, as before) and `cleanToppicDir`. With `-C` the same runs on the
  record files of the earlier run, so their reports are rewritten. Only
  the proteoform ids change: the files' FDRs (the proteoform-level FDR of
  a PrSM is still that of its cluster in the file), `prot_cluster_id`s and
  everything else stay the file's own, so the PrSM and protein selections
  are those of a run without `-c`; the proteoform table picks the best
  PrSM per combined id among the PrSMs that pass, so where the combined
  clustering (`prsm_simple_cluster`) splits or joins the file's clusters
  (`prsm_feature_cluster`) it lists more or fewer proteoforms. TopMG's
  combined runs (`topmg_process`) are unchanged. Results (toppic_test
  experiment 122, 480_F3_02 as two TopFD replicates and as two retention
  time parts at 0-55 / 45-100 min, setting d with `-c`): every row of the
  files' record files, tables and SQLite tables has its combined PrSM's id;
  without the ids the record files, PrSM and protein tables equal main's;
  the combined outputs, `-C` reruns and a run without `-c` equal main's or
  the full run's; the proteoform tables gain 5-6 rows per file where the
  file's feature clustering joined same-protein, same-mass PrSMs of
  barely overlapping ranges (HMGN2 28-56 / 60-90) and lose 1 where it
  compares adjusted precursor masses 1.25 Da apart and the combined one
  the observed ones 1.02 Da apart (PTMA). Combined proteoforms reported
  in both files under one id: replicates 2,002 of 2,002 (main 1,658 of
  1,996), parts 404 of 404 (main 2 of 403). Time unchanged (id update
  ~1 s).
- **The multiple-shift stage runs after the first-pass top selection
  (`TopPIC_identify` → `TopPIC_multi_ptm_search` → `TopPIC_postprocess` →
  `TopPIC_report`, no `TopPIC_post` any more; a combined run's merged set
  goes `TopPIC_multi_postprocess` → `TopPIC_report`).** `TopPIC_multi_ptm_search`
  reads `toppic_pre_raw_prsm` (the first-pass top selection), skips the
  spectra whose best PrSM has an E-value strictly below the skip E-value
  (with `-d` the smallest E-value of a decoy PrSM of the first pass, so
  that every spectrum a decoy could compete with is searched again; 0.01
  without a decoy database; the hidden `--multi-skip-evalue <x>` overrides
  both), searches the rest (`DiagFilterMng::skip_spec_ids_`), merges all
  first-pass PrSMs (`toppic_pre_raw_prsm`, skipped spectra included) with
  the multiple-shift ones and keeps the best PrSM per spectrum
  (`toppic_raw_prsm`; without `-s 2` `TopPIC_postprocess` copies
  `toppic_pre_raw_prsm` to that name). `TopPIC_postprocess` runs the
  clustering, cluster protein correction, recount, post mass matching,
  characterization and FDR once on `toppic_raw_prsm` and writes
  `toppic_record_prsm`;
  `TopPIC_report` runs the decoy removal, the cutoffs and the outputs.
  The PrSM file chain is therefore `toppic_pre_raw_prsm` → `toppic_raw_prsm`
  → `toppic_proteoform_cluster` → `toppic_cluster` → `toppic_cluster_prot`
  → `toppic_recount` [→ the same under the `_post_ms2` name]
  [→ `toppic_cluster_local`] [→ `toppic_cluster_fdr`]
  → `toppic_record_prsm`
  [→ `toppic_prsm_target`] → `toppic_prsm_cutoff` → `toppic_form_cutoff`
  → `toppic_form_cutoff_form` / `toppic_prot_cutoff` → `toppic_prot_cutoff_prot`.
  Do not reintroduce a separate FDR for the multiple-shift PrSMs: with the
  skipped targets missing from its denominator it under-reports (toppic_test
  reports 0023/0024). TopMG's identification is hosted behind the hidden
  `--topmg <disabled|pipeline|topmg-only>` (argument key `topmgMode`,
  default `disabled`; any other value is refused). With `pipeline`,
  `TopPIC_topmg_identify` runs after the multiple-shift stage on the spectra
  the earlier stages left unidentified and hands `toppic_topmg_raw_prsm` to
  `TopPIC_postprocess`; with `topmg-only` (which needs `-b`)
  `TopPICProgress` skips `TopPIC_identify` and
  `TopPIC_multi_ptm_search` altogether, runs `TopPIC_init` (start time,
  argument echo, base data, feature file check) and sets `varPtmTypeNum`
  from the modification file itself, and `TopPIC_topmg_identify` searches
  every spectrum with an empty `raw_prsm_ext`.
- **GUI wiring.** `toppicwindow.ui` has `postMassMatchCheckBox` (checked by
  default from the parser's argument map) on the Basic Parameters tab;
  `gui/util/command.cpp` emits `-D` when the box is **cleared**. The GUI has
  no field for `-I`: `toppicwindow.cpp` erases `postMinPeakNum` from the
  arguments, so toppic uses its default (1).
- **Short options in `toppic_argument.cpp` are nearly exhausted.** Before
  adding one, list the letters in use with
  `grep -oE '"[a-z-]+,[A-Za-z]"' src/console/toppic_argument.cpp`; a clash is
  only reported at run time as "option '-X' is ambiguous" (that is how `-P`
  collided with `--proteoform-ppm-error`). Every option is declared twice
  (the visible `display_desc` and the full parse `desc`); change both.

## ProForma strings (`seq/proteoform_util`, issue #610)

`proteoform_util::getProFormaStr` writes a proteoform in ProForma 2.1
(HUPO-PSI, https://github.com/HUPO-PSI/ProForma; the spec and its compliance
levels are in that repository's `SpecDocument/ProForma_v2.1_Final.pdf`).
`PrsmMatchTableWriter` puts it in the `ProForma` column after `Next amino
acid` (PrSM, proteoform and protein tables of toppic and topmg) and
`PrsmSqlWriter` in the `proforma` column of the SQLite tables; the manual
describes it in `docs/manual/proteoform.md` section 6. Facts to keep
straight (toppic_test experiment 128):

- **PTMs are Unimod name and accession**, `[Oxidation|UNIMOD:35]`,
  `[Cation:Na|UNIMOD:30]` (joined with `|` as synonyms, ProForma 2.1
  section 7.5, a Level 2 feature). The name is `PtmBase::getUnimodName(id)`,
  the abbreviation of the `ptm_base.xml` entry with that Unimod id, never a
  modification file's name (unchecked: `Sodium,...,30` gives
  `[Cation:Na|UNIMOD:30]`); an id without an entry (newer than the file)
  gives the accession alone. A PTM without a Unimod id (a modification
  file's line with UnimodID -1) is its mass with the name as an `INFO` tag
  (section 6.6), `[+12.3457|INFO:MyPTM]`. A name with `|`, `#` or
  unpaired brackets is left out (`isProFormaName`; `Cation:Mg[II]` is
  written). A label (`#s1(0.900)`, `#g1`) goes after the tag's first part
  (`labelTag`): after an `INFO` part it would be read as part of the
  information.
- **`ptm_base.xml` is Unimod**: `No PTM` (TopPIC's only own entry; the
  empty PTM is found by that name, not by mass 0) and every Unimod
  modification (1,560, unimod.obo of 17 Feb 2026) in id order, written by
  toppic_test's `tests/202610/scripts/proforma_output/gen_ptm_base.py`:
  abbreviation = the Unimod (PSI-MS) name, name = the Unimod full name,
  mass = the exact sum of the AME2020 nuclide masses over the Unimod
  composition (bricks expanded). Regenerate it for a new Unimod release;
  do not edit it by hand. Experiment 128 renamed the earlier entries to
  their Unimod names (Carbamidomethylation -> Carbamidomethyl, Oxoalanine
  -> Cys->Oxoalanine, Allysine -> Lys->Allysine, Thiocarboxy ->
  Carboxy->Thiocarboxy, Cholinephosphoryl -> Phosphocholine, TPQ ->
  Quinone; carboxylation and Citrullination were the same Unimod PTMs as
  Carboxy and Deamidated) and removed Decarboxy, which has no Unimod
  entry. The old names are not kept as aliases (the user's choice): a
  FASTA sequence written `C[Carbamidomethylation]` stops with "Unknown
  PTM", a modification file's old name defines a PTM of its own with the
  file's mass, and result XML written with an old name cannot be read
  (result files of other versions are refused anyway, see `topdiff.md`). The
  C57 PTM is found by `getC57AbbrName()` ("Carbamidomethyl").
- **Every Unimod name is built in**, so a modification file's line whose
  name is a Unimod name takes the Unimod entry (`PtmBase::getPtmPtr`,
  documented in `docs/manual/proteoform.md` section 3); when the file's
  mass differs by more than 0.01 Da, a `WARNING:` line is printed once per
  name (on `std::cout`: toppic and topmg set the log level to errors, so
  `LOG_WARN` would not show).
- **Where each mass shift goes**, per alteration (`Alter`), on the
  alteration's residues, or its mass shift's when the alteration has none
  within them (TopMG writes an unexpected shift on one residue with an
  empty alteration range; TopPIC's alterations all have their shift's
  range): a protein
  N-terminal acetylation (PROTEIN_VARIABLE, residue 0, the Acetylation PTM)
  and the truncation N-terminal modification are N-terminal modifications
  `[...]-`; an alteration on one residue follows it (fixed, input, variable,
  or an unexpected shift on one residue); an alteration on several residues
  is a range `(...)[...]`, alterations with the same range sharing it (a
  variable shift's unlocalized PTMs); a local annotation (MIScore, `-B`)
  over several residues is written at its best site with `#s<n>(score)` on
  every site with a positive score (`LocalAnno::getScrVec` is indexed by
  the residue position in the proteoform when its length is the
  proteoform's); an empty or out-of-range alteration is an unknown
  position `[...]?`.
- **Ranges may not overlap or nest** in ProForma (TopPIC's `Proteoform`
  string nests them): the ranges are taken from the left, wider first, and
  one overlapping a taken range becomes a set of positions `#g<n>` with
  its first residue preferred.
- Delta masses: 4 decimals with a sign always (`+0.0000` too).
- **TopMG's variable shifts carry the aligned mass** (an `Oxidation` range
  of 16.2826 Da, a `Phospho;Acetyl;Acetyl;Acetyl` range of 206.0429 Da
  against 206.0980), and its `Proteoform mass` column includes it; the
  ProForma accessions stand for the PTMs' exact masses, so a reader's mass
  of a TopMG string is the theoretical one and differs from that column
  (by up to 0.42 Da on 480_F3_02). TopPIC's variable shifts have the PTM
  masses.
- Checked with pyteomics.proforma (5.0.1) on 480_F3_02 (toppic default
  setting, -B, an N-term variable PTM, topmg, PTMs in the FASTA): every
  string parses, has the proteoform's residues and, for toppic, the
  Proteoform mass within the rounding of the written masses (topmg: after
  the aligned masses above).

## Average masses in msalign files (`ms/env/avg_mass_util`)

`topfd --average-mass` (`-A`, `TopfdPara::output_avg_mass_`) writes the
MS/MS msalign file with average masses: `deconv_ms2_process::deconvMsTwo`
and `topfd_single_process` pass the `DeconvMs` through
`avg_mass_util::getAvgMassMs` just before `MsAlignWriter::writeMs`, which
copies the header (precursor m/z and mass) and the peaks with every
monoisotopic mass replaced by `EnvBase::convertMonoMassToAvgMass` (the
intensity-weighted average of the `theo_patt.txt` distribution, cached in
`Env::avg_mz_`). The SQLite record, the feature files and the MS1 file keep
monoisotopic masses. `toppic --average-mass` (`-A`, as topfd's; `-g`
until 2026-10-12, when `--approximate-spectra` moved from `-A` to `-X`;
argument key `averageMass`) converts each input file in `TopPICProgress_multi_file`,
before any search, with `avg_mass_util::convertToMonoMass` to
`<base>_mono_ms2.msalign` (`EnvBase::convertAvgMassToMonoMass` inverts the
table by fixed-point iteration; the feature file is copied to the new
name) and replaces the spectrum file list with the converted names, so the
per-file runs, the cleanup and the `-c` merge all see monoisotopic files.
Its database is toppic's `<base>_toppic.sqlite`, the same copy of TopFD's
`<base>.sqlite` as without the option (post mass matching section): TopFD's
masses are monoisotopic already (on toppic_test experiment 154's `topfd -D
-A` subset the converted msalign masses equal `ms2_env.mono_mass` within
1e-6 Da, the six decimals). Experiment 154's amendment 1 gave the converted
file a copy of its own, `<base>_mono.sqlite` (main up to then wrote into
TopFD's database with `-A` and wrote no PrSM tables with `-A -D`, which
looked for that file); experiment 156 replaced it with
`<base>_toppic.sqlite`. The results are `<base>_mono_post_ms2_*`
(`<base>_mono_ms2_*` with `-D`; main before fdef25ed wrote
`<base>_post_ms2_*`).
Do not convert masses inside `MsAlignReader` instead: the converted file
is what every later step (and the user) sees. GUI: `toppicwindow.ui` has
`averageMassCheckBox` (Advanced Parameters tab, above the localization
box; the keep-decoy row moved up to make room), read into
`arguments_["averageMass"]`, and `gui/util/command.cpp` emits `-A` when it
is checked.

## Precursor and fragment mass tolerances (`para/sp_para`)

`SpPara` holds two `PeakTolerance` objects, `prec_tolerance_ptr_` and
`frag_tolerance_ptr_`, set by `-E`/`--prec-error-tolerance` and
`-e`/`--frag-error-tolerance` (argument keys `precErrorTolerance` and
`fragErrorTolerance`, both default 10 ppm) in toppic, topindex and topmg and
their GUIs. The index file name (`filter/mng/index_file_name`) carries only
the fragment value, as the integer the searches use (`std::stoi`: `7.5` and
`07` give `7`; until 2026-10-12 the string as given, so `-e 07` missed the
index of `-e 7`): the indexes depend on neither tolerance, and the manuals
say an index is reused when `-f`, `-n`, `-e` and `-d` match. topmg's hidden
`--proteoform-ppm-error` moved from `-E` to `-Q` to free the letter. The
rule for choosing between the two tolerances: the
precursor tolerance covers anything sized on the precursor mass — the
`PeakTolerance::compPrecErrorTole(prec_mass)` checks in the
filters/searches/E-value code (`MsHeader::getFirstPrecErrorTolerance(ppo)`
until the Dalton option below), `refinePrecursorByMatchedIons`, the
localization shift checks
(`LocalMng::prec_tole_ptr_`), **and** the tolerance handed to
`prm_ms_util::getIntMassErrorList` / `getPrmPeakPtrs` and
`comp_prob_value`, which use it only to build the synthetic zero and
precursor PRM peaks (`prm_peak_factory`, sized on the precursor mass; the
real fragment peaks already carry their own tolerances). The fragment
tolerance covers per-peak tolerances set on fragment masses: the PRM /
extend spectrum factories, post mass matching's m/z test, and the
localization M tables (`LocalMng::frag_tole_ptr_`). Neither `SpPara` nor
`PeakTolerance` has XML (de)serialization any more: their `appendXml` and
`XmlDOMElement` constructors were never called and were removed; `SpPara`
is built only from the toppic/topmg arguments in `PrsmPara`.

## Precursor error tolerance in Dalton (`para/peak_tolerance`, `--prec-dalton-error`)

Branch `prec_dalton_error` (toppic_test experiment 133, report 0133, issue
#563; merged into main on 2026-10-07). toppic's hidden `--prec-dalton-error <Da>` (argument
`precDaltonError`, empty by default; checked to be a positive number) sets
`PeakTolerance::setDaltonTolerance` on the precursor tolerance in
`PrsmPara`, and every precursor tolerance is then that value instead of
the precursor mass x the ppm of `-E`. The parameter block prints it in Da.
Facts to keep straight:

- **One rule for a precursor tolerance.** `compPrecErrorTole(prec_mass)` is
  prec_mass x ppo (no minimum) or the Dalton value; it serves the filters'
  and searches' precursor tests, `refinePrecursorByMatchedIons`, the
  E-values' candidate window and isotope offset (`CompPValueArray` keeps the
  tolerance as a member; its functions and `EValueProcessor`'s lost the ppo
  parameter) and the MCMC tolerance. `compStrictErrorTole(mass)` returns the
  Dalton value for every mass but 0: the zero PRM peak keeps the minimum
  tolerance (0.01 Da), which no ppm value ever changed.
  `MsHeader::getFirstPrecErrorTolerance(ppo)` and
  `Precursor::getErrorTolerance(ppo)` are gone, and
  `getFirstPrecMonoMassMinusWaterError` takes the tolerance in Da.
- **Unchanged:** the +-1 Da precursor offsets of the stages, `-e`, `-p` and
  the index. Without the option the 480 run is main's byte for byte (60
  files).
- **The refinement of the shift searches.** `refinePrecursorAndHeaderShift`
  steps over +-2 x tole and pairs fragments within tole + their own
  tolerance. With fewer than two pairs its median index is 0, so the first
  bin wins and the precursor moves by +2 x tole (main 03a0c37a: 264 of
  1,315 one-shift PrSMs of 480 at +2 x 10 ppm; +1.0 Da at 0.5 Da).
    - **Fixed on main by experiment 135** (e720aec2, issue #514): PsAlign and
      PsAlignMultiType call `refinePrecursorByMatchedIons` (the delta in
      [-T, T] that matches the most precursor-dependent ions within the
      fragment tolerance, 0 when none matches). On 480: 3,354 / 1,996 / 440,
      no shift PrSM at +2T. This branch's own fixes (amendment 1: bfa89f11,
      8a5a4582, 75c7ae13) were reverted (3413086c, amendment 2).
    - **TopMG's graph alignment** (topmg and toppic `--topmg`) kept
      `refinePrecursorAndHeaderShift` with the fallback until toppic_test
      experiment 143, which gave it `refinePrecursorByMatchedIons` too and
      removed the old function (section "Precursor refinement of TopMG's
      graph alignment" in `topmg.md`).
    - **In the merge** `refinePrecursorByMatchedIons` takes the precursor
      tolerance instead of its ppo (`compPrecErrorTole`, the same T without
      the option), so `--prec-dalton-error d` sets T = d: the shift searches
      refine the precursor over +-d by the most ions matched with the
      fragment tolerance.
- **Shift E-values depend on the precursor peak width.**
  `CompProbValue::compPrecProbs` multiplies the one-shift precursor
  probability by the precursor peak's width and the i-shift one by the
  width again, so a wider tolerance lowers the E-values of PrSMs with
  unexpected shifts (median log10 change at 0.5 Da: -0.59 one shift, -1.16
  two shifts) and of the decoys among them; zero-shift E-values rise
  (+0.31), from the candidate window.
- **Results on 480** (setting d, 1% FDR; main 3,337 / 2,001 / 437):
  0.05 Da 3,334 / 1,997 / 439, 0.2 Da 3,322 / 1,985 / 443, 0.5 Da 3,316 /
  1,980 / 443; the E-value cutoff falls from 2.86 to 0.80 at 0.5 Da (shift
  decoys). 0.2 and 0.5 Da take 7% longer (zero-shift and variable PTM
  filters). With precursors moved by up to +-0.3 Da, 10 ppm reports 3,174
  PrSMs but only 295 equal to main's (the errors become unexpected
  shifts); 0.5 Da reports 3,057, of them 1,417 equal to main's (1,140 of its
  1,471 zero-shift PrSMs), with an E-value cutoff of 0.042.

## Precursor +-1 Da error per identification stage (`para/sp_para`)

`SpPara` holds one precursor-offset list per stage; each filter/search
expands a spectrum into one spectrum set per offset with
`spectrum_set_factory::geneSpectrumSetPtrVecWithPrecError`. Defaults: zero
shift and one shift `{0, -1 Da, +1 Da}`, variable PTM `{0}`, multi shift
`{0}` (the diagonal filter and PTM search absorb a 1 Da error into the
unknown shift, so they have no switch). `PrsmPara` applies the toppic hidden
options `--disable-zero-shift-prec-error`, `--enable-one-shift-prec-error`
and `--enable-var-ptm-prec-error` (argument keys `zeroShiftPrecOneDaError`,
`oneShiftPrecOneDaError`, `varPtmPrecOneDaError`) through the
`SpPara::set*PrecOneDaError` setters. **toppic's defaults differ from
`SpPara`'s for one stage**: `ToppicArgument::initArguments` sets
`oneShiftPrecOneDaError = "false"` (one-shift stage searched with `{0}`;
`--enable-one-shift-prec-error` restores `{0, -1, +1}`) and
`varPtmPrecOneDaError = "false"` (variable-PTM stage searched with `{0}`,
matching the `SpPara` default; `--enable-var-ptm-prec-error` switches it to
`{0, -1, +1}`, the pre-2026-09 toppic default). Absent keys
keep the `SpPara` defaults, so topmg is unaffected. The fragment-mass
`frag_ext_offsets_` expansion (`{0, -1 Da, +1 Da}` applied by
`extend_ms_factory` to every fragment above `frag_extend_min_mass_`) is separate:
`SpPara::setFragMassExtension(false)` reduces it to `{0}`, reached through
toppic's hidden `--disable-frag-mass-extension` (key `fragMassExtension`,
default `"true"`), and `setFragExtendMinMass` sets the threshold, reached
through the hidden `--frag-extend-min-mass <Da>` (key `fragExtendMinMass`,
default `"5000"`, the `SpPara` default; validated as a non-negative
number).

## E-values of PrSMs found at a +-1 isotope offset (`stat/tdgf/comp_pvalue_array`)

Branch `evalue_adjusted_prec` (toppic_test experiment 117, report 0117;
merged into main on 2026-10-05). `CompPValueArray::process` (the shared mode that every toppic
E-value pass uses: the first pass, the multiple-shift pass, the cluster
protein correction and `--topmg`'s TDGF E-values) groups the PrSMs of a
spectrum by isotope offset. Facts to keep straight:

- **The offset rule** (`getIsotopeOffset`) is the one the candidate number
  already used. A PrSM without unexpected shifts, variable PTMs or a
  truncation N-terminal modification has offset +1 / -1 when its adjusted
  precursor mass is more than the precursor error tolerance above / below
  the spectrum's label-free precursor mass. All other PrSMs have offset 0,
  including the one-shift and variable-PTM offsets of the opt-in
  `--enable-*-prec-error` options, as on main.
- **Offset 0** keeps main's computation: the spectrum set's PRM spectra
  and generating function. Their E-values are unchanged (equal on every
  PrSM of 480_F3_02); a smaller maximum score in the group only changes the
  height of the generating function's table, whose top row saturates.
- **Offset +-1** gets `prm_ms_factory::geneMsSixPtrVec` for the observed
  precursor mass +- 1.00228 Da, the PRM spectra the zero-shift search built
  for that mass, so the generating function's precursor mass is the
  adjusted one. `compMultiExpectedValues` takes the candidate number's
  offset from the PRM spectra's precursor mass instead of the PrSM's
  original precursor, so that the offset is not applied twice. This also
  fixes the separate mode (`compSingleExpectedValue`, unused by toppic),
  which already built the PRM spectra at the adjusted mass and then shifted
  the candidate number once more.
- **Results on 480_F3_02** (`-s 2`, oxidation only):
    - For 94 spectra that experiment 116's oracle searched with the right
      precursor, the new offset E-values equal the direct ones exactly. On
      main the same matches differ by 10^-16.9 to 10^2.6.
    - 1% FDR: 3,336 / 1,996 / 437 against main's 3,337 / 1,997 / 438 on
      TopFD 47b9fc57's output, 3,339 / 2,000 / 437 against 3,340 / 2,001 /
      438 on TopFD ead67ddc's. The E-value pass takes 4-5 s longer (+4%).
    - Decoy E-values move little (log10 ratio within -1.35..1.12); target
      offset PrSMs move by up to 17 orders of magnitude, in both directions.
    - **The best PrSM changes for deamidated spectra.** For an Asn-Gly
      deamidated proteoform whose precursor is right, the unmodified
      zero-shift match through the -1 offset now gets the E-value a right
      precursor gives and wins over the deamidation PrSM (1 shift, (N)[+0.98]),
      with 4-11 fewer matched fragments. Deamidated PrSMs at 1% FDR fall
      from 50 to 45 and from 60 to 53.
    - **Deamidation as a variable PTM repairs only part of it**
      (`Deamidated,0.984016,NQ,any,7` next to oxidation, report 0117
      section 8): 86 / 96 PrSMs with a Deamidated variable PTM against
      main's 91 / 102. In HMGN1 scans 2833, 2838 and 2974 the offset match
      still wins: the generating function on the PRM spectrum of P - 1
      gives its fewer fragments a smaller probability per candidate than the
      one on P gives the variable-PTM PrSM, although the precursor mass is
      4-7 mDa above the deamidated mass and 11-15 mDa below the unmodified
      mass + 1 isotope. Not implemented: choosing by precursor error between
      a PrSM at P and an offset PrSM of the same protein and positions.

## E-values of the variable PTM PrSMs with modified residues (`stat/tdgf/comp_pvalue_array`, `--disable-var-ptm-evalue-mod-residues`)

Branch `var_ptm_evalue_mod_residues` (toppic_test experiment 129, issue
#334), merged into main on 2026-10-06 with the method on by default. The
PrSMs of the variable PTM search get their own generating function; the
hidden toppic option `--disable-var-ptm-evalue-mod-residues` (argument key
`varPtmEvalueModResidues`, default `"true"`; the parameter block prints
"Variable PTM E-values with modified residues: false" only when it is set)
gives them the residues without variable PTMs, as before the merge. Facts
to keep straight:

- **The residues.** Without the method (the other PrSMs, and every PrSM
  before the merge) the generating function draws random sequences from
  the 20 residues with the fixed modifications (C57), each with its
  frequency in the search database (`CountTestNum`), and, for COMPLETE and
  PREFIX PrSMs, the first residues of the protein N-terminal forms. With
  the method, a residue that k - 1 variable PTMs can modify becomes k
  residues, itself and its modified forms, which share its frequency x,
  x / k each (`count_util::geneModResFreqs`): with oxidation on M, M ->
  M and M[Oxidation] at x / 2; with oxidation M, acetyl K, methyl KR and
  phospho STY, K -> K, K[Acetyl], K[Methyl] at x / 3. The PTMs are the
  "any" lines of the `-b` file (`mod_util::readAnywhereModTxt`, set in
  `TdgfMng::var_ptm_mod_ptr_vec_` by `setVarPtmModResidues` in
  `toppic_process.cpp`); N-term lines are not used.
- **Which residue a PTM modifies** follows the variable PTM search
  (`VarPtmAlign::initAllowShiftTable`): the residue must equal the
  modification's original residue (`Residue::isSame`, the amino acid
  without a PTM). A C57 residue takes no variable PTM. The N-terminal
  residue list is split the same way; its residues carry no protein
  N-terminal acetylation (that is an alteration of the proteoform, not of
  the residue), so an acetylated first M splits too, as the search lets it
  be oxidized.
- **Which PrSMs.** `CompPValueArray::useModResidues`: no unexpected shift,
  and variable PTMs or a truncation N-terminal modification, the PrSMs
  whose candidate number already takes the residue sequence mass.
  `process` splits the offset-0 PrSMs of a spectrum into these and the
  others (the +-1 offset PrSMs never are variable PTM PrSMs) and computes
  the generating function once for each group, so a spectrum with both
  computes it twice; the variable PTM group has no unexpected shift and
  needs only layer 0. The other PrSMs' E-values are unchanged: their group
  has a smaller or equal maximum score, which only lowers the saturated
  table height (as in the offset groups above).
- **One table.** `CompProbValue::setResidues` replaces the residue masses,
  frequencies and average residue length before a computation, so one
  CompProbValue (138 MB page table per thread) serves both lists. The
  average residue length of the new list (`count_util::computeAvgLength`)
  is used with it.
- **The candidate number is unchanged**, including the factor
  (variable PTM types + 1) per variable PTM.
- **Passes.** toppic's first E-value pass, the multiple-shift pass and the
  cluster protein correction get the PTMs; TopMG hosted in toppic
  (`--topmg`) does not (its PrSMs come from the graph alignment), and
  neither does topmg.
- **With `--disable-var-ptm-evalue-mod-residues`** toppic gives the files of
  main before the merge (4ecb19ea) byte for byte (56 of 56 files of a 480
  run).
- **Results on 480_F3_02** (`-s 2`, TopFD of main 4ecb19ea, toppic_test
  report 0129):
    - Oxidation only (`var_mods.txt`): 3,337 / 2,001 / 437 against
      4ecb19ea's 3,342 / 2,000 / 437, one decoy oxidation PrSM crossing the cutoff. The
      variable PTM PrSMs' E-values move little (median -0.004 log10 at 5-9
      matched fragments, up to +2.5 for ACBP with 48 fragments).
    - Four PTMs (oxidation M, acetyl K, methyl KR, phospho STY): 3,377 /
      2,013 / 442 against 3,371 / 2,012 / 442, but variable PTM PrSMs 731 ->
      800 and zero-shift PrSMs 1,469 -> 1,425. Their E-values fall (median
      -0.27 log10 at 5-9 fragments, -1.2 at 30 or more, decoys alike), because
      splitting S, T, Y, K, R and M moves 16% of the residue probability away
      from the unmodified residues. The PrSMs with variable PTMs then win ties
      and near-ties against the zero-shift PrSMs of the same spectrum, which
      keep the residues without PTMs. In 19 scans a protein N-terminal acetylation
      becomes an acetylated lysine near the N-terminus, with fewer (12) or
      equal (7) matched fragments. One H2B cluster flips from H2B2F to H2B1C
      + methyl (same mass; the cluster protein correction relabels 18
      scans).
    - The first E-value pass takes +6% (oxidation) and +18% (four PTMs), the
      whole run +1.0% and +1.6%; memory is unchanged.

## Precursor refinement of the variable PTM search (`search/diag/diagonal_util`, `search/varptmsearch/var_ptm_align`)

Branch `var_ptm_refine_prec` (toppic_test experiment 121, report 0121;
merged into main on 2026-10-05). Facts to keep straight:

- **The headers.** `VarPtmSlowMatch::geneVarPtmNTermShiftHeaders` builds one
  diagonal header per variable PTM mass combination, with C-terminal shift =
  observed (label-free) precursor mass - sequence mass - N-terminal shift.
  The C-terminal diagonal is the one whose C-terminal shift is within the
  precursor tolerance, so its C-terminal shift is the precursor error.
- **The refinement.** `VarPtmAlign::geneResult` gets `refine_prec_mass`
  (the proteoform mass plus the variable PTMs of the C-terminal diagonal)
  from `diagonal_util::refineVarPtmPrecursorAndHeaderShift`. That function
  also sets every backtracked header's C-terminal shift from it
  (`DiagHeader::initHeader`, as the slow match built them).
  `diagonal_util::refineVarPtmHeadersBgnEnd` then generates C-terminal ions
  without the precursor error.
- **On main** the headers keep the observed precursor's C-terminal shifts,
  so every y ion of the range refinement carries the precursor error.
    - Scan 5211 of 480_F3_02 is a deamidated proteoform whose precursor
      TopFD's fragment pairs moved by an isotope, 0.0202 Da above the
      proteoform mass. Every y ion below 2,020 Da was lost.
    - The variable PTM shift range spread over the lost ions (report 0120),
      and no ion inside a shift's range is counted.
- **The shift searches** (`PSAlign`, `PSAlignMultiType`, TopMG's graph
  alignment) do their own precursor refinement
  (`refinePrecursorByMatchedIons`, the most precursor-dependent ions within
  the fragment tolerance; until experiments 135 and 143
  `refinePrecursorAndHeaderShift`, the median fragment error). The
  variable PTM search needs none: its proteoform mass is exact once the
  alignment has chosen the variable PTMs.
- **Unchanged:** the dynamic programming (`VarPtmAlign::dp`) and the
  first-pass PrSMs without variable PTMs (31,702 identical on 480_F3_02).
- **Results on 480_F3_02** (`-s 2`, TopFD ead67ddc output):
    - Oxidation: 3,336 / 1,998 / 437 against main's 3,333 / 1,998 / 437.
      Oxidation + deamidation: 3,336 / 1,992 / 436 against 3,331 / 1,989 /
      436. Decoys at fixed E-values are unchanged.
    - First-pass variable PTM PrSMs of the same spectrum and proteoform: 708
      / 743 gain fragments, none lose; E-values better in 250 / 383, worse in
      none. Most gains are at precursor errors of 0.015 Da or more.
    - Reported PrSMs: 15 / 16 of the same proteoform gain fragments (ACBP
      scans 2000-2008 +9 to +13; 5211 38 -> 52, E 9.2e-17 -> 4.2e-33).
      1 / 4 change proteoform (a one-shift PrSM becomes a variable PTM
      PrSM). 3 / 5 PrSMs are new, all borderline (E 0.006-2). None is lost.
- **Remaining (fixed by "Range of a variable PTM shift" below):** the
  range of a variable PTM shift was still the gap between the refined
  diagonals, wider than the modifiable residues (5211's oxidation of M44
  on P43-M44).

## Precursor refinement of the shift searches (`search/diag/diagonal_util`, issue #514)

Branch `c_term_prec_adjust` (toppic_test experiment 135, report 0135; not
merged). Facts to keep straight:

- **Why the refinement matters.** A C-terminal ion of mass m enters the
  alignment as the PRM P - (m - C shift), P the observed precursor. In a
  diagonal anchored at the N-terminus (`n_strict && !c_strict`) it is matched
  with ppo x (m + P), the fragment tolerance plus the precursor's error
  (`prm_ms_factory`, `compDiagPair`); a C-anchored diagonal does the same to
  its N-terminal ions. The PrSM (`Prsm::initScores`) matches the proteoform's
  ions with ppo x m. The refined precursor decides whether the ions the
  alignment matched with the relaxed tolerance (the precursor-dependent
  ions: the C-terminal ions of every diagonal but the last and the
  N-terminal ions of the last) still match.
- **`refinePrecursorByMatchedIons`** (PsAlign, PsAlignMultiType) moves the
  precursor by the delta in [-T, T] (T = P x the precursor ppo, or the Dalton
  value of `--prec-dalton-error`) that matches
  the most precursor-dependent ions within the fragment tolerance (each
  theoretical ion once over the spectra of the group). Within the best
  delta interval it takes the median error of its ions, clipped to the
  interval; of equal intervals the one closest to 0; delta 0 when no ion is
  in reach. The shifts move as before: the last diagonal's N-terminal shift
  and the other diagonals' C-terminal shifts.
- **`refinePrecursorAndHeaderShift`**, main's method until then: the median
  of 0.005 Da bins over +-2T, which can leave the precursor tolerance, and
  +floor(2T / 0.005) x 0.005 Da with fewer than two pairs (report 0133).
  It stayed for TopMG's graph alignment until toppic_test experiment 143
  and is removed ("Precursor refinement of TopMG's graph alignment" in
  `topmg.md`).
- **The variable PTM search needs neither:** its proteoform mass is exact
  (`refineVarPtmPrecursorAndHeaderShift`).
- **Results on 480_F3_02** (`-s 2`, TopFD main 38ba67dc output):
    - First-pass target one-PTM PrSMs: 26,805 C-terminal ions matched with the
      relaxed tolerance, 81.1% of them within the fragment tolerance at P;
      the PrSM keeps 92.8% with main's refinement and 94.9% with this one
      (multiple-PTM: 86.1% and 89.3%). Matched fragments on the same
      alignments +0.9% (decoys +2.4%).
    - Main's fallback hit 251 of 1,281 reported one-shift PrSMs (+2T, about
      +20 ppm in the precursor and the shift mass); none now.
    - With every precursor moved by up to +-5 ppm, 80.1% of 6,779 target
      one-PTM PrSMs get exactly the unmoved run's adjusted precursor (main's
      bins: 0.33 ppm median).
    - 1% FDR: 3,354 PrSMs, 2,001 proteoforms, 438 proteins (main 3,337 / 2,001
      / 437); E-value cutoff 2.49 (2.86); decoys at E <= 1: 26 (19). Same
      running time.
- **Unchanged:** the alignment (it still counts the ions only the relaxed
  window matches), the E-value's generating function (relaxed tolerance on
  the observed precursor).

## The filters' mass list (`ms/factory/prm_ms_util`, `getIntMassErrorList`)

Branch `filter_pairs` (toppic_test experiments 138 and 142, reports 0138 and
0142; merged into main on 2026-10-07 as 5fb46b2c). Facts to keep straight:

- **What the list is.** The PRM peaks of a spectrum (prefix list: ORIGINAL
  m - N shift, REVERSED P - (m - C shift); suffix list: the same minus H2O)
  as integer masses round(x x `filter_scale_`, 100 = 0.01 Da) with the
  tolerance ceil(tol x 100) (ORIGINAL ppo x m, REVERSED ppo x (m + P)), plus
  the zero and precursor entries, sorted. The users add 1 per entry to every
  column in its window: `MassMatch::compScores` (one-shift and variable PTM
  filters), `compDiagMaxScoresImpl` (multi-shift filter) and
  `CompShiftLowMem::findBestShift` (the diagonal finder of `PtmSlowMatch`).
  The zero-shift filter does not use it; TopMG uses it in its one-shift
  filter (and the multi-shift one with `--use-asf-diagonal`).
- **One entry per integer mass and base type.** Peaks of one type with the
  same integer mass give one entry (the largest tolerance); an ORIGINAL and a
  REVERSED peak of the same integer mass are two. Main kept one entry per
  integer mass: at an accurate precursor a b ion and the PRM P - y of its
  complementary y ion differ by (precursor error - fragment error) x the
  masses, a few mDa, and often shared an integer mass, so the right
  proteoform scored one point for the pair (two when a bin boundary fell
  between the peaks). On 480_F3_02 an identified spectrum lost 10.9 entries.
- **Repeated integer masses are fine for the users**: they loop over the
  entries, and `findBestShift` takes its range from the first and the last
  entry (the precursor entry). Main already had them for combined spectra
  (the merging restarts per spectrum).
- **Results** (toppic_test experiment 142, against main 7b5112c2):
  480_F3_02 setting d 3,369 / 2,003 / 441 against 3,354 / 1,996 / 440, more
  targets at every rank of the target-decoy list; the TMT E. coli set of
  experiment 80 444 / 263 / 69 against 454 / 244 / 67 (more targets at most
  ranks; a lower 1% cutoff); TopMG in toppic (`--topmg topmg-only`) 3,025 /
  1,768 / 422 against 3,010 / 1,765 / 420. No change where the precursor and
  fragment errors differ by 10 ppm (P1_high_high_85min_1: nothing merges) or
  on st_1. No time or memory cost. Right proteoforms near a filter's top-N
  cutoff can still drop out when other candidates gain (4 of 23 changed
  spectra on 480).
- **The rule is a hidden option**, `--filter-mass-list
  <merge|fragment|single>` of toppic and topmg (argument key
  `filterMassList`, default `fragment`; branch `filter_mass_list`,
  toppic_test experiment 142, amendment 2; merged into main on 2026-10-08
  as abe33295, after `second_align`, `ps_align_penalty` and `local_align`).
  `SpPara::FilterMassListMode`
  (`kMerge`, `kFragment`, `kSingle`; `setFilterMassListMode` from the
  name, `isFilterMassListModeName` for the parsers' check) is set by
  `PrsmPara` and passed to `getIntMassErrorList` by its four callers, so a
  new caller must pass `getFilterMassListMode()` too. `merge` is main's
  rule up to 7b5112c2 (all peaks in one slot), `fragment` the rule above,
  `single` one entry per peak. Both parsers refuse other values, and the
  parameter block prints `Filter mass list:` only for `merge` and
  `single`, so default output is unchanged. toppic's `--topmg` stages take
  it from toppic's argument map. Checked on st_1 and 480_F3_02 (setting d;
  every output file but the time lines and that parameter line, 64 per
  run): the default and `fragment` write main c7f13d66's and filter_pairs
  57e44368's files, `merge` main 7b5112c2's, `single` filter_all_masses
  f6572a10's; topmg's one-shift filter files on st_1 equal those of the
  three binaries mode by mode.
- **`single` is not better** (toppic_test experiment 142, amendment 1,
  against `fragment`): its extra entries (61 per identified 480_F3_02
  spectrum, +13.6%) are 99% one fragment mass at its other charge states,
  which TopFD lists separately, so a fragment counts once per charge state.
  480_F3_02 3,357 / 2,000 / 437 against 3,369 / 2,003 / 441 (fewer targets
  at most ranks, strong PrSMs lost at the filters' cutoffs); TopMG in
  toppic fewer targets at every rank; topmg 2,017 against 2,052; TMT 525
  against 444 PrSMs, mostly because two high-ranked decoys are not found;
  P1_high_high 708 against 689 (a looser cutoff). The multi-shift filter,
  whose time grows with the square of the list, takes 34% longer on 480 and
  TMT, the runs 5-8%.

## Tag filter (`filter/tag`, `--tag-filter`)

Branch `tag_filter` (toppic_test experiment 150, report 0150; from main
d56ebc44, not merged). A sequence-tag protein filter whose candidates
(simple PrSMs) join the candidate lists of toppic's identification stages.
Off by default; the hidden toppic option `--tag-filter <stages>` (argument
`tagFilterStages`, a comma-separated subset of `var`, `one`, `multi`,
`topmg`, or `all`; any other name is refused) turns it on for the named
stages, with the hidden `--tag-filter-length <n>` (`tagFilterLength`, the
residues of a tag, default 4, at least `TagIndex::K`) and
`--tag-filter-cand-num <n>` (`tagFilterCandNum`, the proteins kept per
spectrum set, default 10). The parameter block prints the three only when
the filter is on. Facts to keep straight:

- **Tags** (`TagFilter::geneTags`). The PRM peaks of a spectrum set
  (`prm_ms_util::getPrmPeakPtrs` on `getMsTwoPtrVec`: each fragment mass
  as the N-terminal PRM, ORIGINAL, and as the precursor minus the mass,
  REVERSED, plus the zero and precursor peaks), per base type, peaks within
  the sum of their strict tolerances of each other merged (one fragment at
  several charge states). Two peaks whose difference is a residue mass (the
  19 letters with the fixed modifications; I is read as L) within the sum of
  their strict tolerances are an edge, one per residue that fits (K and Q
  both at a wide tolerance). A tag is a path of `tag_len_` edges in
  increasing mass, so it reads the protein N to C whichever series it comes
  from; the strings are looked up once each. About 180 tags per 480_F3_02
  spectrum (1,100 without the merging of charge states).
- **The index** (`TagIndex`). The search database (targets and decoys, the
  file the searches read) as one text with I as L, the positions of every
  4-mer in a bucket per 4-mer (as `prsm/search_fasta_match`), and the prefix
  mass before every residue (`residue_util::convertStrToResiduePtrVec` with
  the fixed modifications, so FASTA PTMs count; one double per residue,
  180 MB for the human database). A tag is looked up by its rarest 4-mer
  and verified at each position; every occurrence in every protein is
  returned. Built once per stage, 2-3 s.
- **Score** (`TagFilter::scoreHit`). An occurrence lies on the diagonal
  s = (PRM mass of the tag's first peak) - (prefix mass before the tag);
  the occurrences of a protein sorted by s form a diagonal while
  neighbours are within the sum of their first peaks' N-terminal
  tolerances (`getNStrictCRelaxTolerance`: the precursor's error enters a
  REVERSED peak's PRM). A diagonal's support is the number of distinct
  peaks its tags cover; the score is the best diagonal's support plus the
  second's when that reaches `TagFilterMng::getMinScore()` (tag_len + 3,
  default 7: a ladder of six residues or two tags apart on one diagonal),
  and a protein below the minimum is not a candidate. The first version
  counted the peaks of all occurrences: titin, nesprin and the other
  longest proteins collected hundreds of chance 4-mers and filled every
  list (37% decoys in the var lists, 7 of 174 confident zero-shift
  proteins of the first 500 spectra listed); with the diagonals and the
  minimum 83 of 174 are listed (one stage 88), 0.4% / 33% decoys. Ties go
  to the shorter protein, then the name (a name tie put every `DECOY_`
  before `sp|`). The top `cand_num_` proteins are the candidates.
- **What a candidate carries, by stage** (`TagFilterMng::Stage`; the
  spectrum sets of the stage's precursor error list, `getPrecErrorVec`, so
  the search pairs them by precursor mass as it pairs the filter's):
    - `one`: the diagonals' shifts (the best supported first, at most
      `max_shift_num_` = 10, two within 0.01 Da once) as `n_term_shifts`
      and, as c = P - protein mass - s (`getNTermShiftsFromCTermShifts`
      inverts it), as `c_term_shifts`; the same candidate in the four type
      files (`OnePtmSlowMatch` reads the lists its type needs).
    - `var`: `geneVarPtmPrsms` finds the proteoforms containing the tag
      whose mass matches the precursor: the N-terminus a protein N-terminal
      form of `-n` (`allowMod` on the protein; s = the form's `getProtShift`
      plus a sum of 0 to `-S` variable PTM masses within the first peak's
      tolerance) or a truncation start i >= 1 (s = -prefix[i] + such a sum,
      plus a Position N-term PTM mass with one PTM fewer and i >= 2), the
      C-terminus the end e with prefix[e + 1] - prefix[start] + water +
      the form's `getPepShift` or the N-term PTM + a sum of PTM masses = P
      within `compPrecErrorTole`. One simple PrSM per (start, end, form,
      N-term PTM): the range, `n_term_shifts` {the form's shift} (0 for
      truncations), `c_term_shifts` {0}, the truncation N-terminal
      modification shift when used, in the file of its type (complete /
      prefix for forms, suffix / internal for truncations). `VarPtmSlowMatch`
      exits on a form shift that is not an allowed `ProtMod` and reads
      `getCTermShifts()[0]`, so both lists must be set.
    - `multi` and `topmg`: the protein alone (`PtmSlowMatch` and the graph
      alignment find their own diagonals).
- **Files and joining.** The stage's mass filter runs as before; the tag
  filter (`tag_filter_processor::process`, a thread pool over the spectra,
  `OrderedTaskRunner` so the files are in msalign order) writes
  `toppic_tag_filter_<stage>[_<TYPE>]` (`toppic_topmg_tag_filter` for the
  hosted stage); `tag_filter_merge::process` writes `<filter>_tag[_<TYPE>]`,
  the name the search then reads: the filter's candidates first, unchanged,
  then the new tag candidates of the spectrum (`one`: a protein listed with
  the same precursor variant gets the tag shifts added to its lists;
  `var`: the same protein, range, N-terminal shift and truncation
  modification is not repeated; `multi`: the same protein and precursor
  variant is not repeated). For `topmg` the tag file is one more input of
  the `SimplePrsmStrMerge` that builds `toppic_topmg_graph_filter` (the
  best 20 x 5 instead of 20 x 4), so the graph alignment's score-10 rule
  (`SimplePrsmFilter`) applies to the tag candidates too. The multi stage's
  tag filter skips the identified spectra (`skip_spec_ids_`), the hosted
  stage's scores the spectra TopMG searches (`search_spec_ids_`);
  `cleanToppicDir` deletes the files. The console prints "Tag filtering
  (<stage>) - started / finished", the index size, a summary (spectra,
  tags, occurrences, candidates) and the joined counts.
- **Without the option** the run is main d56ebc44's (the first 500 spectra
  of 480_F3_02, setting d: every file; the toppic manual is unchanged, the
  options being hidden).
- **Results on 480_F3_02** (toppic_test experiment 150, report 0150; the
  branch's TopFD, setting d `-s 2 -b var_mods.txt`, 1% FDR, 16 threads; base
  3,374 PrSMs / 2,004 proteoforms / 440 proteins in 19:28): `--tag-filter
  var` 3,376 / 2,006 / 441 (+44 s; 1,578 proteoform candidates for 1,197
  spectra, 1.0% decoys), `one` 3,386 / 2,012 / 441 (+4:09: 24,774 proteins
  for 2,649 spectra, 33% decoys, the one-shift search 120 s against 51 s,
  and a one-shift decoy at E 4.2e-4 lowered the multiple-shift stage's skip
  threshold from 0.0215, so it searched 1,337 spectra instead of 1,112),
  `multi` 3,382 / 2,013 / 441 (+1:02), `var,one,multi` 3,389 / 2,015 / 443
  (+6:11): 23 spectra gained, 17 of them through proteins only the tag
  filter listed (HMGN3 59-99 with 38 fragments at E 1.3e-12, PFD1 42-66,
  TM263 1-37, RL12 124-165, ...), 8 lost at the cutoff (2.41 against 2.50),
  217 PrSMs of the same protein re-explained (112 better E-values, 57
  worse; RL38 2-70 from E 5.2e-5 to 6.2e-46). The hosted TopMG search
  (`--topmg topmg-only`, oxidation, phospho, N-terminal acetyl): 2,996 /
  1,764 / 407 against 2,979 / 1,751 / 406, 17 gained and none lost, but
  35:04 against 16:40, the graph alignment taking 1,482 s against 420 s on
  the 16,537 tag candidates above its score 10 (whole, long proteins). A
  tag filter stage costs 50 s on 4,085 spectra (about 180 tags, 85 k
  occurrences per spectrum). The tag filter lists the protein of about half
  of base's final PrSMs (53% var, 54% one, 48% multi, 60% TopMG; three
  quarters of those with 30 or more matched fragments, a tenth of those
  with 10-19): tags need five consecutive matched cleavage sites, so the
  filter complements the mass filters on strong spectra whose protein the
  index misses and cannot replace them.

## Second alignment of the shift searches (`search/diag/second_align_util`, `--disable-second-align`, issue #614)

Branch `second_align` (toppic_test experiment 144, report 0144 and its
amendment 1; merged into main on 2026-10-08 in the form below). On by
default; toppic's hidden `--disable-second-align` (argument key
`secondAlign` = "false"; `PrsmPara::useSecondAlign`; the parameter block
prints "Second alignment at the refined precursor: false" only then)
keeps the first alignment's PrSM, main's behaviour before the merge.
Facts to keep straight:

- **Why.** The one-shift and multiple-shift alignments (PsAlign,
  PsAlignMultiType) match the precursor-dependent ions with the fragment
  tolerance plus the precursor's ("Precursor refinement of the shift
  searches" above), and the PrSM's
  annotation matches them with the fragment tolerance at the refined
  precursor, so the dynamic program counted ions the annotation did not
  and could choose a diagonal set on their strength (report 0135 section
  7: "A second alignment on the refined precursor with the strict
  tolerance would make the alignment and the annotation agree").
- **What it does** (`PsAlign::geneResult`, `PsAlignMultiType::geneResult`).
  The first pass is the former `geneResult`, factored into
  `second_align_util::genePrsm` (the precursor refined by the matched
  ions to P', the extended spectra at P', the ranges refined, the PrSM).
  Then, for every call with a first alignment: the PRM peaks are rebuilt
  at P' (`genePrmPeaks`: `prm_ms_factory::geneMsSixPtrVec` with the
  refined precursor, so the REVERSED peaks are P' - (m - C shift), plus
  the zero and precursor peaks sized on P'; every fragment peak's
  N-strict-C-relax and N-relax-C-strict tolerances set to its strict
  tolerance ppo x m); the first alignment's headers are moved to P'
  (`moveHeaders`: an N-anchored header keeps its N-terminal shift, a
  C-anchored one keeps its C-terminal shift, so its N-terminal shift
  moves by P' - P; C-terminal shifts and truncation positions from
  `DiagHeader::initHeader` at P'); the diagonals are rebuilt
  (`geneDiagonals`, headers without a pair dropped); the same dynamic
  program runs on them (PsAlign for the type of the last `compute()`,
  which it now remembers; PsAlignMultiType for the types of the last
  `compute()` together, kept by P' in `second_align_ptrs_` so the types
  whose first results refine to the same mass share one second
  alignment); the second result for the same type and shift number goes
  through `genePrsm` with `keep_prec` (the refinement's delta range is
  [0, 0], so P' is kept and the headers keep their shifts), the range
  refinement and the PrSM. That PrSM is returned; with no second
  alignment of positive score the first PrSM is kept, so no candidate is
  lost. The candidates are ranked by the PrSM's matched fragments as
  before, and the E-values are computed from the PrSM afterwards (the
  generating function is unchanged). The zero-shift search, the variable
  PTM search and TopMG's graph alignment are untouched. The precursor
  never leaves the window of the first refinement (+-T around the
  observed precursor).
- **`refinePrecursorByMatchedIons` takes an optional delta range**
  (`delta_min` / `delta_max`, infinite by default), which cuts its
  [-T, T]; only the second pass uses it, with [0, 0].
- **Research log.** `TOPPIC_SECOND_ALIGN_LOG=<file>` appends one
  tab-separated line per `geneResult` call with a first alignment: stage
  (`one` / `multi`), spectrum id, type, shift number, protein, P, the
  first score, P', the first PrSM (start, end, matched theoretical
  fragments, shifts as left:right:mass), the second score (`-` without a
  second alignment), P' again, the second PrSM, and `second` / `first`
  for the PrSM kept. The alignment's score is the sum of the diagonal
  pairs' scores (one per spectrum and ion type per matched cleavage) with
  the zero and precursor peaks counted with the group's spectrum number
  each, so score - 2 against the PrSM's matched fragments measures the
  agreement (for the multiple-shift search add the large-shift penalty
  back, 3 per shift above 300 Da).
- **Variants tested and dropped** (experiment 144, branch commits
  3dc127c4-26f7dbe3, kept in the branch's history): a second refinement
  by the matched ions after the second alignment (+-T around P', so up to
  2T from the observed precursor: 14% of the one-shift candidates and 68
  reported shift PrSMs of 480_F3_02 ended beyond 10 ppm, and one such
  move lost two ATIF1 oxidation PrSMs through the cluster protein
  correction), the same clamped to the observed precursor's window, and
  a second alignment with the first alignment's relaxed tolerances
  (which left the alignment and the annotation apart: the tolerance is
  the gap, not the precursor).
- **Results on 480_F3_02** (setting d, experiment 131's TopFD output,
  against main 7b5112c2): the second alignment's score is within 1 of the
  PrSM's matched fragments for 98.4% of the 265,195 one-shift and 98.3%
  of the 41,712 multiple-shift candidate alignments (the first alignment:
  81.9% and 69.3%). It changes the positions or shifts of 14% and 19% of
  the candidate PrSMs (the shift positions and diagonal sets; 695
  one-shift candidates gain matched fragments, mostly +1, and 297 lose).
  1% FDR: 3,356 PrSMs / 1,996 proteoforms / 439 proteins against 3,355 /
  1,997 / 438, the same E-value cutoff (2.49) and decoys: within the data
  set's noise, because the candidates that reach 1% FDR are chosen by the
  annotation's fragments and the E-value either way, and the relaxed
  window inflated targets and decoys alike. The one-shift search takes
  38 s instead of 23 s; the multiple-shift search's 256 s are unchanged.
  The pushed main 249433b1 (this, the penalty fix below, `filter_pairs`,
  `topmg_prec_matched_ions` and `local_align`) gives 3,375 / 2,004 / 440
  against the pre-merge main c7f13d66's 3,372 / 2,008 / 439, the same
  cutoff (report 0144 amendment 2).

## Large-shift penalty of the shift alignments (`search/oneptmsearch/ps_align`, `search/ptmsearch/ps_align_multi_type`)

Branch `ps_align_penalty` (toppic_test experiment 145, report 0145; merged
into main on 2026-10-08 as 0bca905d, after `second_align` e61e339b; on the
development subset the merged main with `--disable-second-align` differs
from main c7f13d66 only in the one-shift stage files and what follows them,
as the branch did from its base). `PsAlignPara` fixes `align_large_shift_thresh_` = 300 Da and
`align_large_shift_panelty_` = 3: a shift between two diagonals above the
threshold costs 3 points in the dynamic program. Facts to keep straight:

- **The two programs applied it differently until this branch.**
  `PsAlignMultiType::dp` (the multiple-shift search since 4f845a6f)
  compares the shift predecessor's score minus the penalty with the
  diagonal and truncation alternatives and stores it, so the penalty
  propagates and the end score is 3 lower per large shift. `PsAlign`
  (the one-shift search; the multiple-shift search with
  `--disable-multi-type-align`; the original TopPIC code) used the
  penalty in `getShiftPre` only to choose which diagonal to shift from,
  and `dp` then took the chosen predecessor's score without it: a large
  shift was never discouraged against continuing without a shift, and
  the stored and final scores carried no penalty. `ps_align_multi_type.hpp`
  describes the two as the same recursion; they were not (toppic_test
  report 0144 section 7, found by comparing the alignment scores with
  the PrSMs' matched fragments).
- **This branch makes PsAlign match PsAlignMultiType**: `getShiftPre`
  returns the chosen predecessor's penalized score, and `dp` compares and
  stores it. The predecessor choice is unchanged; what changes is that a
  path through a large shift now competes with the no-shift and
  truncation alternatives 3 points down, and the end score reflects it.
- **Results on 480_F3_02** (toppic_test experiment 145, report 0145;
  setting d, experiment 131's TopFD output, against main 7b5112c2): 3,354
  PrSMs / 1,997 proteoforms / 437 proteins against 3,355 / 1,997 / 438,
  the same E-value cutoff. 825 of the 14,670 best one-shift candidate
  alignments change, 672 of them a shift above 300 Da replaced by a
  smaller one (the same matched fragments in 593, 1-3 fewer in 206, more
  in 26; targets and decoys alike); the multiple-shift stage files are
  unchanged. 31 of 3,353 reported PrSMs get another proteoform, mostly a
  large shift rewritten as a smaller one with the same or better
  fragments and E-value; 2 MED9 PrSMs are lost when the feature's YBOX1
  +14 Da PrSM takes the cluster (the penalty took 2905's -421 Da
  alignment with 19 fragments down to a -49 Da one with 18) and SEM7A at
  the protein cutoff; 1 gained. No time cost.

## Tolerance of the complementary PRM peaks in the E-values (`stat/tdgf/prob_peak`, `--evalue-complement-tolerance`, issue #616)

Branch `evalue_strict_complement` (toppic_test experiment 146, report 0146;
from main 90d744d7; merged into main on 2026-10-08 as a hidden option with
the relaxed tolerance, main's behaviour, as the default). Facts to keep
straight:

- **The three tolerances of a PRM peak** (`prm_ms_factory::addSixMasses`,
  the fragment tolerance ppm): an ORIGINAL peak (a b-ion mass m) has
  strict = N-strict-C-relax = ppm x m and N-relax-C-strict = ppm x (m + P);
  a REVERSED peak (the complementary mass P - m of a y-ion, P the
  precursor the PRM spectrum was built at) has strict = N-relax-C-strict
  = ppm x m and N-strict-C-relax = ppm x (m + P), the precursor's error
  added because P - m moves with P. The searches use the relaxed
  tolerances for the precursor-dependent ions (`compDiagPair`), the
  annotation and the second alignment the strict one at the refined
  precursor.
- **The generating function before the branch** (`ProbPeak`, the "strict"
  flag false in toppic; the separate mode with the PRM spectra at the
  PrSM's precursor, `CompPValueArray::compSingleExpectedValue`, is not
  used) took every peak's N-strict-C-relax tolerance: the strict window
  for an ORIGINAL peak and the relaxed one, about twice as wide, for a
  REVERSED peak, on the six-mass PRM spectra of the spectrum set at the
  observed precursor P (or P +- 1 isotope for the PrSMs found there).
- **The option**: toppic's hidden `--evalue-complement-tolerance
  <relax|strict|separate>` (argument `evalueComplementTolerance`, default
  `relax`; anything else is rejected). With `strict`, `ProbPeak` takes
  the strict tolerance for a REVERSED peak too
  (`PrsmPara::useEvalueStrictComplement`, passed to `CompProbValue` by
  `CompPValueArray`; the cluster protein correction's recomputation goes
  through the same `PrsmPara`); the zero and precursor peaks are ORIGINAL
  and unchanged, so the precursor peak's width, which normalizes the
  shift E-values (`compPrecProbs`), stays the precursor tolerance. The
  parameter block prints "E-value tolerance of complementary ions:
  strict" (or separate) only when set. topmg has no option (its own
  E-values are the MCMC's; `--topmg` in toppic follows the setting).
  Without the option the run is main's before the merge byte for byte
  (development subset: 63 files but the time stamps).
  `separate` (6f612533) selects the code's separate mode instead
  (`PrsmPara::useEvalueSeparate`, `EValueProcessor::process(true)` in
  both E-value stages and the cluster protein correction): one
  generating function per PrSM on PRM spectra rebuilt at the PrSM's
  adjusted precursor, every peak strict, the proteoform mass as the
  precursor. It costs nothing (the E-value stage's time is the table
  setup per spectrum, 176 s on 480 either way).
- **The generating function's resolution** (report 0146 section 7):
  `TdgfMng::convert_ratio_` is 27.4335 units per Da (0.0365 Da; the
  searches use 274.335), residue masses are rounded to integers (proline
  +0.50 unit), and a window is ceil(tolerance x 27.4335) units, so a
  strict window of 10 ppm x m is one unit for every fragment below 3.65
  kDa and a relaxed one, 10 ppm x (m + P), one unit below 3.65 kDa and
  up to six at 20 kDa. The strict setting therefore changes nothing for
  precursors below 3.6 kDa and most for the heaviest, and a PrSM whose
  ions form a long ladder on a proline- or D/E-rich stretch can get a
  probability of 0 ("Invalid e value", the placeholder 1e300) because no
  integer path hits all its one-unit windows (NPM scans 3786 and 3617,
  FNBP4 6285 on 480; the separate mode has them too).
- **Results on 480_F3_02** (setting d, experiment 131's TopFD output,
  against main 8688fd08): strict 3,374 PrSMs / 1,994 proteoforms / 441
  proteins, separate 3,373 / 1,997 / 444, base 3,375 / 2,004 / 440.
  Every E-value falls (targets median x 10^-0.45, decoys x 10^-0.31,
  more for heavier precursors and higher scores) and the cutoffs with
  them (PrSM 0.84 against 2.50); 8 spectra report a one-shift form of
  the same protein instead of main's plain form, the multiple-shift
  stage skips 41 more spectra because its skip threshold (the smallest
  first-pass decoy E-value, experiment 131) falls from 0.021 to 0.009
  and 7 two-shift PrSMs are lost that way (ROA2 5019 with 39 fragments,
  JUPI1 4202 with 40). Not merged; the recommendation is `relax` as the
  default if the option is merged, and a finer generating-function
  resolution before the strict tolerance is tried again.

## Range of a variable PTM shift (`seq/mass_shift`, `search/varptmsearch/var_ptm_align`)

Branch `var_ptm_narrow_range` (toppic_test experiment 122, report 0122;
merged into main on 2026-10-05). Facts to keep straight:

- **The ranges.** `VarPtmAlign::geneShiftVec` takes the variable PTM shifts
  from the gaps between the refined diagonals
  (`diag_header_util::getDiagonalMassChanges`) and narrows each alteration
  to the residues of its gap that can carry its PTM. It now also calls
  `MassShift::narrowRangeToAlters`, which sets the shift's range to the
  union of its alterations' ranges (alterations covering no residue are
  skipped; with none left the range stays as it is). A shift with several
  alterations (a diagonal between two PTMs was empty) covers all of them.
- **All residues of the PTM** (amendment 1). An alteration covers every
  residue that one of the modifications of its PTM modifies:
  `VarPtmSearchMng::ptm_mod_ptr_vec_2d_` keeps them in one list per PTM
  (grouped by PTM, not by mass), `getPtmModPtrVec(mod)` returns it, and
  `geneShiftVec` tests each. `Deamidated,0.984016,NQ,...` gives Deamidated
  on N and on Q; the alignment puts a tied PTM on the first N or Q of its
  stretch (the `>` tie rule of `VarPtmAlign::dp`) and the backtrack takes
  that residue's modification, so testing only it (as before) left the
  other type out of the range (FOXK1 4009 `AVT(Q)[Deamidated]NGK`, now
  `AVT(QN)[Deamidated]GK`). On 480_F3_02 (ox + deamidation) 3,247 first-pass
  ranges and 8 reported ones change; search fragments, E-values and
  identifications do not, and oxidation-only runs are identical.
- **What the range decides.** `Proteoform::getSeqSegmentPtrVec` generates
  no ion inside a shift's range, and the proteoform string and the
  "Variable PTMs" column show it (`P(M)[Oxidation]`, `Oxidation:[44]`
  instead of `(PM)[Oxidation]`, `Oxidation:[43-44]`).
- **No ion of the search spectrum is gained.** The range refinement trims
  each diagonal to its matched ions in the same extended spectrum, so the
  ions the narrowing adds were already unmatched there. On 480_F3_02 the
  first-pass PrSMs keep their counts and E-values (10,791 / 11,854 variable
  PTM PrSMs, a Python copy of `findPairs` predicts every count).
- **Fragments are gained after the search.** The post mass matching looks
  for every unmatched theoretical fragment of the best PrSM in TopFD's raw
  peaks, and the recount uses every mass. With more theoretical ions both
  find more: of the best PrSMs of every spectrum, 203 / 231 gain 482 / 541
  final fragments (88-89% of them masses the post mass matching added);
  none loses. E-values come from the search, so they do not change, except
  for a PrSM that the cluster protein correction gives another spectrum's
  proteoform (the correction recomputes its E-value): 1 of 4,075.
- **Identifications are unchanged** (oxidation 3,336 / 1,998 / 437,
  oxidation + deamidation 3,336 / 1,992 / 436), as are the decoys and the
  run time. Of the reported variable PTM shifts, 438 / 480 (ox) and 543 /
  592 (ox + deamidation) are on one residue (main: 143 and 208); every
  other one covers two or more residues that can carry its PTMs.
- **Determined sites** (branch `var_ptm_split_sites`, toppic_test
  experiment 124, report 0124; merged into main on 2026-10-05). A shift holds several PTMs when the
  diagonal between them had no matched fragment. When each of its PTMs has
  as many candidate residues as copies (and no residue is a candidate of
  two PTMs), the sites are known: `VarPtmAlign::splitDeterminedShift` makes
  one shift per PTM on its residue, with that residue's modification
  (`(MAQM)[Oxidation;Oxidation]` -> `(M)[Oxidation]AQ(M)[Oxidation]`,
  `Oxidation:[33];Oxidation:[36]`). Any other shift is left as it is. The
  fragments between the sites cannot match in the search spectrum (that
  diagonal was empty), so on 480_F3_02 no first-pass count or E-value
  changes (1,589 / 2,451 first-pass PrSMs split); 11 reported PrSMs per
  run change their notation, 3-4 of them gain a final fragment.

## Deamidation or the -1 isotope offset (`prsm/prsm_deamid_selector`)

Branch `deamid_prec_error` (toppic_test experiment 123, report 0123;
merged into main on 2026-10-05). **Off by default**; the hidden toppic option `--deamid-selector`
(key `deamidSelector`, default `"false"`; printed in the parameter block
only when set) turns it on. Then, with a deamidation in the variable PTM
file (a PTM with Unimod id 7 or a monoisotopic mass within 0.001 Da of
0.984016 Da; `getDeamidPtmNames`), `TopPIC_identify` runs
`prsm_deamid_selector::process` between the E-value computation and the top
PrSM selection: it reads `toppic_evalue` and writes `toppic_evalue_deamid`,
which `prsm_top_selector` then reads. Facts to keep straight:

- **Why.** The unmodified proteoform matched at the -1 offset (a zero-shift
  PrSM whose adjusted precursor mass is P - 1.00228) and the deamidated one
  at P (a variable PTM PrSM) explain a spectrum with masses 18.3 mDa apart.
  Since experiment 117 their E-values come from generating functions on the
  PRM spectra of different precursor masses, and the offset match often
  wins with fewer fragments.
- **The pairs.** Per spectrum (all its precursors, as the top PrSM selection
  groups them) and per protein, start and end: the zero-PTM PrSMs
  (`unexpected_ptm_num` and `variable_ptm_num` 0) and the deamidation PrSMs
  (no unexpected shift, a `Variable` alteration whose PTM abbreviation is a
  deamidation, read with `MassShift(XmlDOMElement)`). The best of each kind
  by E-value is compared. The pairs are zero-PTM PrSMs at the -1 offset in
  practice: no other offset gives the two masses within the tolerance.
- **The rule:** matched fragments (`match_theo_frag_num`, the E-value's
  score) 5 or more apart decide. Otherwise the smaller |precursor error|
  decides: the label-free precursor mass + k x 1.00228 - the adjusted
  precursor mass (the proteoform mass in both searches), k = round((adjusted
  - precursor) / 1.00228). A tie keeps the better E-value. Every PrSM of the
  losing kind in the group is removed, so that another deamidation site
  cannot take the loser's place. E-values are not changed.
- **The loser may have been the best PrSM, the winner need not be.** When
  the zero-PTM PrSM wins with a worse E-value than the deamidation PrSM, the
  top PrSM selection can pick a third PrSM (on 480: one-shift PrSMs of the
  same protein).
- **Results on 480_F3_02** (`-s 2`, oxidation and deamidation, against main
  7c8e00e6; TopFD 47b9fc57 / ead67ddc outputs; run with the step on, before
  the option made it opt-in):
    - The step: 539 / 544 pairs, half of them decoys, 620 / 621 PrSMs
      removed, 3.4 s. A Python copy of the rule predicts
      `toppic_evalue_deamid` exactly.
    - 3,335 / 1,994 / 436 and 3,338 / 1,997 / 436 at 1% FDR (main 3,335 /
      1,993 / 436 and 3,338 / 1,996 / 436), the same cutoff and decoys.
    - It recovers report 0117's five HMGN1/HMGN2 deamidations (2488, 2752,
      2833, 2838, 2974) on both inputs and keeps 4679's offset match on D.
    - It loses 7 other deamidations per input: PTMS 33-80 (5,553 Da; five
      scans on S, six on D) and PTMA (12 kDa) by precursor error, LMNA 3734
      (S) by fragments; one-shift PrSMs are reported instead. Deamidated
      PrSMs at 1% FDR 94 -> 92 and 102 -> 100.
- **The cause of the PTMS losses:** TopFD's precursor masses are
  systematically high by about +1.3 to +1.6 ppm (+7 mDa at 5.5 kDa, +17.5
  mDa above 10 kDa; the median of the confident zero-shift PrSMs), so the
  deamidated explanation (+10 to +13 mDa) looks worse than the -1 offset
  (-5 to -8 mDa). Measured from the run's median ppm error instead of zero,
  the deamidation is closer in all five (report 0123's simulation; not
  implemented).

## Truncation N-terminal modifications (`-b`, Position `N-term`)

A variable PTM file line with Position `N-term` is a variable PTM of the
proteoform N-terminus (the alpha-amine of its first residue), kept apart
from the side-chain PTMs of that residue. Facts to keep straight:

- **Representation.** `Proteoform::trunc_n_term_mod_ptr_` (default
  `ModBase::getNTermNoneModPtr()`) records it, and its mass is the residue
  sequence's `n_mod_ptr_` (sequence mass, `BpSpec` prefix masses), as for an
  acetylated protein N-terminal form; it is no mass shift.
  `setTruncNTermModPtr` accepts it only on a proteoform with the NONE
  protein N-terminal form starting at protein residue 3 or later (0-based
  start >= 2; residues 1-2 are the protein N-terminal forms') whose first
  amino acid is the mod's. `proteoform_factory::geneTruncNTermModProteoform`
  makes the variant, `geneSubProteoform` keeps it when the first residue
  stays, the PrSM XML has `<trunc_n_term_mod>` only when set.
- **Output.** `[Name]-` prefix of the proteoform string, outside any
  annotation of residue 1; `Name:[N-term]` first in `Variable PTMs`
  (`getVarPtmStr`); `getVarPtmNum` (`# Variable PTMs`, XML
  `variable_ptm_num`, FDR groups) counts it; one `prsm_mass_shift` row of
  type `N-terminal variable`.
- **Filter.** One more round per N-term mass m (spectrum built as for an
  N-terminal label of m), only on diagonal-index rows of a fitting start
  (suffix and internal types; the N-terminal form index is not scored), up
  to `-S` - 1 other PTMs, each round its own `SimplePrsmMemMerge` group. A
  candidate carries the masses to try in `SimplePrsm::trunc_n_term_mod_
  shifts_` (`{0}` = none, not written to the file then); a candidate of
  several rounds is listed once with the union.
- **Search.** The alignment knows no N-term shift: the processor searches
  the candidate for each listed mass, 0 as it is, otherwise the variants of
  the N-term mods of that mass, with at most `-S` - 1 side-chain PTMs
  (`VarPtmAlign::max_ptm_num_`). So residue 1 can carry an N-terminal and a
  side-chain PTM. When a variant and the candidate as it is score the same
  (one acetylation of K1, on the N-terminus or the side chain), the
  processor's sort (`cmpScoreDecProtIncTruncNTermFirst`) reports the
  N-terminal one; the user chose this over comparing E-values.
- **E-value.** (T + 1)^k x (T_N + 1) for a PrSM with k side-chain PTMs and
  a truncation N-terminal modification, T_N = the number of N-term masses
  (`TdgfMng::n_term_var_ptm_type_num_`, argument `nTermVarPtmTypeNum`); the
  candidate count uses the sequence mass without the N-terminal mod.

## Proteoform classification (`seq/proteoform`, the PrSM tables)

Branch `proteoform_classification` (issue #611; toppic_test experiment 129
and its amendment 1, report 0129; merged into main on 2026-10-06). `Proteoform::getClassification(
prot_hit_num, var_mod_ptrs)` gives the level of the proteoform in the
five-level classification of proteoform identifications (Smith et al., Nat.
Methods 16:939-940, 2019), by the rule the user gave:

|  | one protein | several proteins |
|---|---|---|
| no unexpected shift, variable PTM sites decided | `1` | `2D` |
| no unexpected shift, a variable PTM site undecided | `2A` | `3` |
| unexpected shift(s) | `4` | `5` |

Facts to keep straight:

- **The properties.** An unexpected shift is an alteration of type
  `UNEXPECTED` (`getAlterNum`); protein N-terminal forms, fixed PTMs, the
  truncation N-terminal modification and TopMG's FASTA PTMs (`INPUT`) do
  not count. The protein hits are `SearchFastaMatch::process`'s matches,
  the `# Protein hits` column; 0 counts as one protein. A variable PTM site
  is undecided (`hasUndecidedVarPtmSite`) when a variable mass shift covers
  more than one residue and one of its PTMs does not have exactly as many
  candidate residues as alterations, or two PTMs share one: the candidates
  of a PTM are the residues of the range that one of its modifications in
  `var_mod_ptrs` modifies (Deamidated on N and Q, Phospho on S, T and Y),
  the test of `VarPtmAlign::splitDeterminedShift`. toppic splits a shift
  whose sites are decided (one PTM per residue), so in toppic every
  multi-residue variable shift is undecided; TopMG's shifts span the gaps
  between matched fragments, and most of its ranges hold one candidate
  residue (decided). The proteoform knows neither the database nor the PTM
  file, so both come from the caller.
- **One place computes it.** `prsm_cutoff_selector::process(...,
  var_mod_ptrs)` sets it for every selected PrSM without one, from
  `SearchFastaMatch` on the database file the tables search (built only
  when a PrSM needs it); `TopPIC_report` passes the "anywhere"
  modifications (`mod_util::readAnywhereModTxt`) of `-b` and `-B`,
  `TopMG_post` those of `-i`. The proteoform stores it (`setClassification`,
  `getClassification()` without arguments; "" until set). The proteoform
  and protein cutoffs and the top PrSM filters read the PrSM cutoff's PrSMs
  and keep it, and the tables and the SQLite writer, which read the cutoff
  files, write the stored level (`-` if none).
- **Where it is written.** The last column, `Proteoform classification`, of
  every table of `PrsmMatchTableWriter` (PrSM, proteoform and protein,
  single and with the other matches; toppic, its combined `-c` tables and
  topmg), `proteoform_classification`, the last column of the SQLite tables
  `prsm`, `proteoform` and `protein`, and `<classification>` in the
  `<proteoform>` element of the XML outputs (`_prsm.xml`, `_proteoform.xml`,
  `_protein.xml`), after `<proteo_db_seq>` and the truncation N-terminal
  modification. The column is last so that no existing column moves; the
  rows of a PrSM's other protein matches have no field for it (they end
  before the q-values, 39 of the 45 fields before it). The element is
  written and read only when set, so the PrSM files before the report are
  unchanged.
- **What the rule leaves out.** Levels 2B and 2C never occur. An unexpected
  shift is level 4 whatever its mass, also when it is a known PTM on one
  residue (HMGN1 `(S)[+79.9669]`) or of isotope size. With `-B`, a shift
  that the localization explains by a PTM becomes a `VARIABLE` alteration
  (`stat/local/local_proteoform`), so its proteoform is level 1, 2A, 2D or
  3 (not run). Several proteins can be one gene (LAP2A / LAP2B, both TMPO),
  and that is still 2D. The level is that of the PrSM's proteoform, so the
  PrSMs of one proteoform id can differ (one with an unexpected shift, one
  without), and the proteoform table shows the level of its best PrSM.
- **Results** (toppic_test experiment 129, 480_F3_02, `-f C57 -d -t FDR
  -T FDR`, oxidation): toppic (`-y FDR -s 2`) writes main 4ecb19ea's files
  but the column, the element and the SQLite column (1de20d8b, levels 1 /
  2D / 4 / 5, against main; a2a84e16 against 1de20d8b: the report stage's
  PrSM files equal without `<classification>`, the tables and SQLite
  tables equal but 29 PrSMs, 16 proteoforms and 1 protein moved from 1 to
  2A), in the same time (16:48 against 16:42; the PrSM cutoff 1.4 s
  against 0.3 s, the protein hit index). The level agrees with a Python
  copy of the rule in every PrSM row of toppic's and topmg's tables, and
  the XML element and the SQLite column with the tables. Proteoforms of
  levels 1 / 2A / 2D / 3 / 4 / 5 at 1% FDR: toppic 1,050 / 16 / 70 / 0 /
  812 / 50; topmg (one unexpected shift by default) 1,007 / 15 / 68 / 1 /
  21 / 1, where 210 of its 228 range shifts hold as many candidate
  residues as PTMs (decided).

## PTM localization by a variable PTM search (`stat/local_align`, the default; `--local-original` for `stat/local`)

Branch `local_align` (toppic_test experiment 140, report 0140; from main
6eb40e97, merged into main as 927e2da3 on 2026-10-08). Since 2026-10-08
`LocalAlignProcessor` is **the default** `-B` localization (argument
`localAlign`, "true"); the hidden `--local-original` selects the MIScore
candidate search of `stat/local` (`LocalProcessor`, the method before),
and the hidden `--local-align` (the switch that selected the new method
while it was optional) is still accepted. The parameter block prints
"Localization method: Variable PTM search (local_align)" or "MIScore
candidate search (local, --local-original)". Without `-B` neither runs.
What the method does, and how it differs from `stat/local`:

- **Explanations of a shift.** For an unexpected shift of mass m and the
  precursor tolerance e (`compStrictErrorTole(adjusted precursor mass)`),
  the explanations are the single PTMs X of the `-B` file with
  |mass(X) - m| < e and a modifiable site, and the unordered pairs {X, Y}
  (X = Y allowed) with |mass(X) + mass(Y) - m| < e and two sites in one of
  the two orders. Unlike `local_util::getPtmPtrVecByMass`, a match one
  isotope mass away does **not** count by default (switch `iso`, 0 since
  experiment 140's amendment 1: with it on, a +12.98 Da shift of a PrSM
  with a +1.00 Da second shift became a Methyl at -1 Da, i.e. a 1 Da
  precursor error turned into a PTM). The shift is replaced only when
  there is **exactly one** explanation; with none or several it stays
  (logged "none" / "ambiguous"). Several explanations are first reduced
  by the **tie-break rules** of `local_align_util::breakTies` (switch
  `tie`, 1; experiment 140's amendment 5), applied in this order while
  more than one is left: (1) fewer PTMs (Dimethyl over Methyl+Methyl,
  Trimethyl over Methyl+Dimethyl, a PTM on the N-terminal acetylation
  form over Acetyl + the PTM on the plain form), (2) the N-terminal
  acetylation form over an Acetyl of `-B` on the plain form of the same
  start (terminal candidates); (3) **Acetyl vs Trimethyl** (0.036 Da
  apart, within 10 ppm above 3.6 kDa): when the explanations left are on
  the same proteoform and the same up to Acetyl <-> Trimethyl, both are
  localized and the better match score wins (then the higher confidence,
  then Acetyl; `LocalAlignProcessor::choose`, logged
  `site_scores(Acetyl:31/Trimethyl:33)`) - a first version of amendment
  5 preferred Acetyl outright, which gave calmodulin's K115
  trimethylation (758, 88.5%) an Acetyl at 48.9%. What is left after the
  rules must be one explanation. Without the rules (`tie=0`) Dimethyl is
  never chosen when Methyl is listed too, nor Acetyl when Trimethyl is.
  The PTMs are recognized by the abbreviations `Acetyl` and `Trimethyl`
  of the `-B` file. The log's `rule` column names the rules that decided
  and `all_explanations` the set before them.
- **No filters.** `LocalProcessor` keeps a localization only when the
  matched fragment number stays above the original's minus 5 and above
  0.67 x the original (`DESC_MATCH_LIMIT_`, `desc_ratio_`), when a site
  passes the MIScore threshold `-H`, and (two shifts) when the localized
  PTMs overlap the original shifts; `LocalAlignProcessor` applies none of
  these. When no site reaches `-H`, the PTM is reported on the range from
  the first to the last site with a positive score, with the **best
  site's score** as the confidence (`local_align_util::filterScores`;
  `LocalAnno::getConf`, the `<confidence>` of the PrSM XML). Above the
  threshold the confidence is `scrFilter`'s, the sum of the scores of the
  reported range; the sum over the fallback range would be 1.0 whatever
  the spread, so it was replaced by the best site's score in experiment
  140's amendment 1 (a Methyl over 18 K/R of H2B at 5.5% each: 0.055, not
  1.0). The per-site scores are in the MIScore column either way.
- **The proteoform's termini.** `stat/local` scores the candidate forms
  of `local_proteoform::getAllCandidateForms` (N- and C-terminal
  truncations by up to 5 residues, the N-terminal acetylation forms) and
  takes the best; `local_align` first explains the shift value of the
  proteoform as it is. When that gives no explanation and the shift has
  at most `term_max_matched` (3; 2 until amendment 4) matched cleavage
  sites on its N-terminal (C-terminal) side (`compMatchedSites` of the
  original proteoform), the candidate proteoforms whose start (end) is
  moved by up to `term_range` (5) residues are tried
  (`LocalAlignProcessor::processTerminal`, the forms of
  `local_proteoform::getCandidateForm`, so a new start of 0 or 1 brings
  the N-terminal acetylation variants; the piece keeps at least one
  residue, and a truncation may cut into the shift's range - until
  amendment 4 the shift's residues had to stay inside, which skipped the
  truncations that explain a shift sitting on the terminal residues, 31
  of main's moved-start localizations on 480). The other N-terminal forms
  of the **unmoved** start (the acetylation dropped, or added: main's base
  candidates include them, and +37.95 Da + 42.01 is a Phospho) belong to
  the **standard attempt** (`processPiece`, `n_forms`): when the piece
  holds the N-terminus and the shift has at most `term_max_matched`
  matched sites on its N-terminal side, their explanations join the pool
  of the proteoform's own before the tie-break rules, so that a PTM on the
  acetylation form competes with Acetyl + the PTM on the plain form
  (amendments 4 and 5; until 8476656b they were terminal candidates, tried
  only after an "ambiguous" standard attempt had already stopped). Each
  candidate's unexplained mass is the adjusted
  precursor mass minus the candidate's mass, and the one-explanation rule
  counts over all (candidate, explanation) pairs ("term_localized",
  "term_none", "term_ambiguous" in the log, the explanations tagged
  `@start<n>` / `@end<n>` and the N-terminal form). The result proteoform
  takes the moved terminus (the base form from `getCandidateForm` with the
  kept shifts, the new shifts moved with `shiftPositions`); a two-shift
  PrSM can move its start for the first piece and its end for the second.
  Switch `term` (1; 0 = the proteoform as it is only). Added in
  experiment 140's amendment 3.
- **The site tables** (`compOnePtmSTable`, `compTwoPtmMTable` in
  `local_align_util.cpp`) are `local_util`'s with the known shifts of the
  proteoform (the variable PTMs, and the other unexpected shift of a
  piece) in the theoretical prefix and suffix masses, a cleavage site
  inside a known shift's range giving no fragment (as
  `Proteoform::getSeqSegmentPtrVec` does for the matching);
  `local_util::compOnePtmSTable` ignores the known shifts. The protein
  N-terminal modification of a proteoform read from the PrSM file (the
  acetylation) is the residue sequence's N-terminal mod
  (`ResidueSeq::getNModPtr`), not a residue's mass nor a shift that
  `Proteoform::getMass` adds; `makePiece` adds it to the prefix masses of
  a piece starting at the N-terminus (as `BpSpec` does) and to the N
  offset of any other piece - missed until experiment 140's amendment 7,
  so the N-terminal fragments of acetylated proteoforms never matched
  (13 of 39 acetylated PrSMs localized by both main and the branch had
  different site scores, none of the 194 plain ones). The fragment
  match is a binary search of the sorted extended masses per theoretical
  mass (not `compMTable`'s two-pointer walk, which needs monotonic
  theoretical masses, broken by a negative known shift). The scores are
  `prob_ratio_ ^ (score - best score)` before the normalization, so that
  a score above ~230 does not overflow `pow`.
- **Two PTMs for one shift**: both orders are scored with the DP of
  `LocalProcessor::compTwoPtmScr` and the better order is localized by
  the break-point method of `twoPtmLocalize` (ported, with the site pairs'
  probabilities summed directly). The two VARIABLE shifts take the ranges
  of their sites above the threshold, `[bgn_1, end_1]` and
  `[bgn_2, end_2]`. `local_proteoform::createProteoformPtr` builds its
  alterations with the break point as one end, but `Alter::setLocalAnno`
  resets an alteration's range to its annotation's `[bgn, end + 1)`, so
  main's ranges are the same tight ranges and match the score vectors
  that `Proteoform::getMIScore` and the ProForma writer index from the
  alteration's left position. Experiment 140's amendment 2 first read
  this as a latent out-of-range read and padded the vectors to the
  break-point ranges (4cdb0069), which emptied the MIScore of the second
  PTM of every two-PTM PrSM (13 on 480); the padding was reverted and the
  comment in `createProteoformPtr` states the rule.
- **Two unexpected shifts.** `stat/local` explains the sum of the two
  shifts by one ordered pair over the whole proteoform; `local_align`
  splits the proteoform at a cleavage site b between the shifts (first
  shift's right position <= b <= second's left position): the middle one
  (by index) of the sites in that range matched by a fragment of the
  spectra (`compMatchedSites`, the original proteoform with all its
  shifts), or the middle of the range when none is matched. The piece
  before b holds the first shift, the piece from b the second
  (`makePiece`; a known shift belongs to the piece of its left position
  and is cut at the piece's end), and each piece is explained on its own
  with the residue masses and shifts of the other piece as its fragment
  offsets (`Piece::n_offset`, `c_offset`), so that the pieces are matched
  against the real spectra. A PrSM can come out with one shift explained
  and the other kept (switch `split`, 0 = two-shift PrSMs left alone).
- **The result PrSM** keeps the proteoform's residues, N-terminal form,
  truncation N-terminal modification, start / end and cluster ids; the
  explained shifts are replaced by the VARIABLE shifts, the PrSM takes the
  proteoform with `setProteoformPtr` (matched numbers recounted) and the
  adjusted precursor mass becomes the proteoform's mass, as `stat/local`
  does. A one-shift PrSM whose shift is within [-1, 1] Da + e is left
  alone (`stat/local`'s rule).
- **Alternative explanations** (`Proteoform::alternative_mass_shift_list_`,
  `MassShiftPtrVec2D`; toppic_test experiment 140, amendment 8). When a
  shift has several explanations, the other candidates on the same
  proteoform as the chosen one (the proteoform as it is when none is
  chosen: the shift stays, the alternatives are all its explanations) are
  localized too (`LocalAlignProcessor::choose`, `PieceResult::alternatives`)
  and the result proteoform gets one complete mass shift list per
  alternative (`setAlternatives`): its shifts other than the unexpected
  ones and the variable PTMs the localization added, the alternative's
  variable PTMs, and any other unexpected shift as it is - so for a
  two-shift PrSM each list explains one shift and keeps the other as a
  shift. Candidates on a moved terminus or another N-terminal form than
  the result's cannot be expressed as a mass shift list and are left out.
  PrSM XML: `<alternative_mass_shift_lists>` with one `<mass_shift_list>`
  per alternative under the proteoform (read back by `parseXml`; older
  files have none). PrSM TSV: the column "Alternative PTM explanation"
  after MIScore (`getAlternativePtmStr`: each alternative's variable PTMs
  and unexpected shifts in position order, `Name:[a-b]` / `+mass:[a-b]`
  joined by `;`, the alternatives joined by ` | `). The SQLite PrSM table
  does not have it.
- **Research aids.** `TOPPIC_LOCAL_ALIGN_PARA` ("name=value,..."): `iso`
  (0; 1 = the isotope match of `stat/local`), `max_ptm` (2; 1 = single
  PTMs only), `split` (1), `term` (1), `term_range` (5),
  `term_max_matched` (3), `tie` (1; 0 = no tie-break rules).
  `TOPPIC_LOCAL_ALIGN_LOG=<file>`: one tab-separated line per shift
  handled (spectrum id, PrSM id, protein, shift index and mass, piece with
  the break site and the matched sites of the range, the explanations,
  the decision small / none / ambiguous / localized / no_sites / overlap,
  the sites with the confidence and the score).
- **Results** (toppic_test experiment 140, report 0140; 480_F3_02 in
  setting d plus `-B` with Oxidation M, Acetyl K, Methyl KR, Dimethyl KR,
  Trimethyl K, Phospho STY; the same 3,354 / 1,997 / 438 identifications
  with equal E-values in every run, the stage 2-3 s): the branch replaces
  400 shifts of 311 PrSMs (main 278 of 265): 222 one-shift PrSMs (116
  Phospho, 42 Phospho+Phospho, 36 Methyl ...) and 109 pieces of the 163
  two-shift PrSMs (76 Oxidation; 20 PrSMs with both shifts explained, 69
  with one); 44 one-shift PrSMs are ambiguous (Dimethyl / Methyl+Methyl
  19, Acetyl / Trimethyl 16), 42 localized shifts match one isotope away,
  29 PTMs are on the fallback range. 148 PrSMs are localized by both
  methods, 145 with the same PTM(s) and 141 with the same best site(s).
  Main alone localizes 117 (76 by truncated candidate forms, 12 by
  dropping the N-terminal acetylation form, 29 by choosing among
  ambiguous explanations); the branch alone 163 (63 one-shift PrSMs by
  pairs, 24 by single PTMs that main's filters rejected, 86 pieces of
  two-shift PrSMs).

## Index files version 2, the combined filter pass and the E-value stage (toppic_test experiment 150)

Branch `speedup_four_ptm_s2` (toppic_test experiment 150, report 0150;
from main 915bb387, merged into main on 2026-10-09). `toppic -u 16 -f C57
-d -t FDR -T FDR -y FDR -s 2 -b topmg_var_mods.txt` on 480_F3_02 takes
543 s instead of 816 s (-33%) with every output file byte-identical to
main 915bb387's (the `-k` stage files, the record PrSM file, the PrSM
XML). Facts to keep straight:

- **Index files version 2** (`MassMatch::serializeMassMatchV2`): the
  header `TPMMIDX2` (scale, proteoform / column / row numbers), then the
  seven arrays as zstd -1 frames, the column bounds stored once
  (`col_index_bounds_`, col_num + 1 ints) instead of the begins and ends.
  `deserializeMassMatch` probes the magic and reads a version 1 Boost
  archive or a version 2 file (one `read` into a thread-local buffer,
  `posix_fadvise DONTNEED` after it). `topindex --index-format <1|2>`
  (hidden, default 2). The human reviewed database with decoys takes 93 GB
  instead of 172 GB. Version 1 index directories stay readable; the file
  names (`filter/mng/index_file_name`) are unchanged, so the format is
  known only by reading a file. `MassMatch::compScores` prefetches the
  column bounds and row ids a few masses ahead (no measurable change).
- **One pass for the zero-shift, variable PTM and one-shift filters**
  (`filter/varptm/zero_var_ptm_filter_processor`): with a variable PTM
  file and no approximate spectra, `toppic_process` runs one stage, "Zero
  unexpected shift, variable PTM and one unexpected shift filtering" ("Zero
  unexpected shift and variable PTM filtering" with `-s 0`), that loads
  each block's four zero-shift index files once and scores them with the
  `ZeroPtmFilter`, a `VarPtmFilter` built on the loaded indexes (a new
  constructor) and, with `-s` of at least 1, the `OnePtmFilter`. The
  one-shift filter's index files are the zero-shift ones cut at its column
  limit (`truncateColumns`), so it scores the loaded indexes with its
  windows cut there (`MassMatch::compScores(..., col_limit)`), and where a
  spectrum's mass lists equal the variable PTM filter's (both
  `getIntMassErrorList` at precursor offset 0) and every window lies below
  the limit, it finds its candidates on the variable PTM filter's score
  arrays in place (`VarPtmFilter::computeBestMatch`'s `base_scores_hook`,
  called before the precursor bonus is added). The three stages' output
  files are unchanged (`toppic_zero_filter_*`, `toppic_var_filter_*`,
  `toppic_one_filter_*`), so the searches and the tag filter's hooks read
  them as before; the one-shift index files are not read in this mode
  (topindex still writes them for the separate stages and topmg). Hidden
  `--separate-var-ptm-filter` runs the separate stages; without a variable
  PTM file the zero-shift stage runs alone as before. The zero-shift
  filter computes its N-terminal index scores once per spectrum and copies
  them for the three precursor offsets (the N-terminal mass lists do not
  depend on the precursor); the stage's three spectrum readings run in
  parallel threads. Copying the score arrays between the filters was
  tried and cost more than it saved (cache eviction); read them in place.
- **The multiple-shift search's diagonal finder**
  (`CompShiftLowMem::findBestShift`, `findBestShiftChunked`) counts the
  mass differences in chunks of 65,536 bins from the top of the range
  down, keeps the count events of at least a threshold (4, then 3, then 2)
  in the original loop's order per spectrum mass (a window across a
  chunk's lower edge waits in a pending buffer until the lower chunk has
  emitted its part) and replays them through `checkShift`; when the
  replay does not fill the 20-shift list the next threshold is tried, and
  last the original loop runs. The lists are the original's (an event
  below the final list's smallest count never changes it; checked on a
  dump of 5,046 calls, `TOPPIC_SHIFT_DUMP=<file>`); the search went from
  157 to 35 s, because the original loop touched a random cache line per
  mass pair in a buffer of up to 770 MB.
- **The E-value stages** (`stat/tdgf`): a PrSM with at most
  `comp_evalue_min_match_theo_frag_num_` (4) matched fragments gets the
  maximum E-value without a table and is left out of its spectrum's
  groups; a spectrum whose PrSMs are all such skips the generating
  function. `CountTestNum` is built in a background thread started before
  the filters (`EValueProcessor::prepareCountTestNum`, joined when the
  stage starts). The PrSMs are written in spectrum order through the
  `OrderedTaskRunner` (16 pending tasks per thread; fewer let one slow
  spectrum idle the others) with `PrsmXmlWriter::writeTrimmedLines`, the
  bytes the former per-thread files and merge produced.
- **The generating function in score slices** (`CompProbValue`): the
  table (height = the group's best score + 1, at most 128, x the mass in
  units of 1/27.4 Da) is swept in slices of 8 rows (`kChunkRows`), each a
  sweep over all columns in blocks of 1,350 columns through a circular
  page of 6 blocks, which needs the slice below only at the columns with
  a position score (`saved_rows_`); the look-back window of the residues
  stays in L2, where the whole-height sweep (a page of 5,000 Da x 128
  rows, 138 MB per thread) streamed the table through memory and the stage
  was bound by memory bandwidth. The per-cell operations and their order
  are the old sweep's, but every slice of every layer starts from a
  cleared page, so the columns before column 0 hold zeros (branch
  `evalue_clean_sweep`, toppic_test experiment 152). **The old sweep
  carried a layer's last pass into the next layer**: it cleared the tail
  of its page once per computation, so for precursors whose table wraps
  the page (above about 4.6 kDa) a shift layer's first blocks read,
  through the wrap, the last 5,105 columns (the largest residue mass) the
  previous layer left there: random sequences with layer 0's scores at
  some mass near the end of the table, placed at mass 0 of layer 1, which
  then collect the spectrum's matches a second time. Only the one-shift
  and two-shift probabilities of such precursors were raised (the
  zero-shift layer and the variable PTM PrSMs are layer 0): by a median
  of 10^0.00002 at 5-9 matched fragments and 10^0.0001 at 10-19, but by
  up to 10^8.7 for PrSMs with 20 or more (one-shift targets of the first
  pass above 4.6 kDa); decoys by at most 10^0.02. Main 3b8041dd
  reproduced it (`tail_capture_`) to keep main 915bb387's output; on
  480_F3_02 (setting d with the four variable PTMs) the clean sweep gives
  the same 3,412 PrSMs / 2,016 proteoforms / 447 proteins and cutoffs,
  and 5 spectra report another PrSM of the same protein (a one-shift
  explanation that now beats a zero-shift PrSM at a +-1 isotope offset
  or an isobaric variable PTM, in 4 of them with more matched fragments).
- **Research aids**: `TOPPIC_PROFILE=<file>` (`common/util/
  stage_profile`: thread-seconds per named part and stage; the output
  files do not change) and `TOPPIC_SHIFT_DUMP=<file>` (the diagonal
  finder's inputs).
- **Results** (480_F3_02, toppic_test report 0150, one run at a time):
  main 820 / 812 s -> 542 / 543 s, user time 10,880 -> 7,510 s, system
  time 529 -> 128 s, peak memory 18.9 -> 18.4 GB; 3,412 PrSMs, 2,016
  proteoforms, 447 proteins at 1% FDR in every run. Not kept:
  binary-search shift predecessors in the alignments (the quadratic loop
  was 1% of the search), per-class table heights in the E-value.

## Rescoring of the PrSMs (`stat/rescore`, `--rescore`)

Branch `prsm_rescore` (toppic_test experiment 151, report 0151 and its
amendment 1; from main 3b8041dd, merged into main on 2026-10-09 as the
hidden `--rescore`, off by default). **Off by default**: toppic's
`--rescore` (argument key `rescore`, "false") turns it on, and
`validateArguments` refuses it without a target-decoy search (`-d`) and
together with `--fdr-groups`. From 2026-10-09 (2863fa08) to 2026-10-10 it
was on by default with `-d`, with the hidden `--disable-rescore` to turn it
off; at the user's request it is opt-in again (`rescore_opt_in`, toppic_test
experiment 155, amendment 3), and `--disable-rescore` is gone (a run of it
in the toppic_test reports is a run without `--rescore` now). Since
2026-10-10 (`rescore_visible`, experiment 155, amendment 5) `--rescore` is
a visible option: `toppic -h` lists it after the cutoffs (short form `-q`
since amendment 6, `rescore_short_option`; `toppic_gui` still emits
`--rescore`), and the toppic manual describes it (section 2, steps 13
Rescoring and 14 FDR computation; the `-q`, `--rescore` row of section 4;
an example in section 6). The parameter
block prints "FDR computation: rescored PrSMs (linear SVM, 3-fold
cross-validation; the maximum E-value last), the E-value at a level where
it gives more" only when it runs (the two additions only while their rules,
below, are on; this line and `--fdr-groups`' "FDR computation: separate
PrSM groups" are printed with or without a variable PTM file). It gives
the best PrSM of every spectrum a score of a semi-supervised linear model,
and the FDRs rank by it instead of the E-value. GUI wiring: `toppicwindow.ui` has
`rescoreCheckBox`, "Rescore PrSMs (needs the decoy database)", on the
Advanced Parameters tab under "Keep intermediate files" (500, 248),
unchecked by default from the argument map and read back into
`arguments_["rescore"]`; it cannot be checked while "Decoy database" is
unchecked, "Decoy database" cannot be unchecked while it is checked, and
`checkError` refuses the pair before toppic is started;
`gui/util/command.cpp` emits `--rescore` only when it is checked (the
hidden `--disable-rescore-fallback` and `--disable-rescore-evalue-floor`
have no GUI control). Facts to keep straight:

- **Where it runs.** `TopPIC_postprocess` keeps its order (clustering,
  cluster protein correction, recount, post mass matching, PTM
  characterization) and runs `prsm_rescore::process` on the result
  (`toppic_rescore`) just before the FDR; `TopPIC_multi_postprocess` (`-c`)
  does the same on the merged PrSMs and spectra (not tested). The E-values
  are not changed; `prsm_fdr::process(..., use_rescore)` sorts by
  `PrsmStr::cmpRescoreDecProtInc` and its three FDR functions compare the
  key `rescoreKey` (the score negated; `eValueKey` by default), so the
  cutoffs (`-t FDR`, `-T FDR`, `-y FDR`) apply to the rescored FDRs.
  `-t EVALUE` cutoffs still compare E-values.
- **The score** (`stat/rescore/rescore_model`): Percolator's algorithm
  (Kall et al., Nat. Methods 2007). The spectra are split into 3 folds at
  random (`std::mt19937` with seed 1 on the spectrum ids in increasing
  order, `rng() % 3`); for each fold a linear SVM
  (`stat/rescore/linear_svm`: L2 penalty, squared hinge loss, C = 0.1,
  "balanced" class weights, dual coordinate descent, the bias as a
  penalized feature; scikit-learn's `LinearSVC` defaults) is trained on the
  other two folds 5 times, each time with the targets at 1% FDR of the
  current score (TopPIC's ratio, `passFdr`) as the positives and every
  decoy as a negative, the first iteration ranking by -log10 E-value; the
  features are standardized on the training folds. The held-out fold's
  scores are put on one scale (the training folds' 1% threshold -> 0, their
  median decoy -> -1). A fold with fewer than 10 positives gives every PrSM
  of the search its -log10 E-value (the FDR is then the E-value's). The
  searches (`<tool>` toppic, topmg) are rescored apart, as the FDR counts
  them apart. The `<rescore>` element of the PrSM XML (after
  `<protein_fdr>`, written and read only when set: `Prsm::getRescore`,
  `PrsmStr::getRescore`, NaN without) carries the score to the FDR and
  the outputs; the TSV tables do not show it.
- **The features** (`prsm_rescore::compFeatures`, names in
  `getFeatureNames`): -log10 E-value (clipped to [-5, 100]), log10 test
  number, -log10 one-protein probability, matched fragments and masses (the
  record's, after the recount and the post mass matching), matched masses
  per spectrum mass and fragments per residue, log10 masses of the
  spectrum, 0 or 2 unexpected shifts, variable PTMs, log10(1 + total
  |shift|), a shift within 0.02 Da of a common modification mass, NME,
  acetylated N-terminal form, start <= residue 1, end at the protein's
  last residue, log10 length and mass, precursor charge, log10 intensity,
  intensity / the largest precursor's, precursor number, log10 feature
  intensity, ECScore; and from the b / y (c / z) pairs of
  `peak_ion_pair_util::genePeakIonPairs` on the spectra read with the
  EnvCNN cutoff of `-F` (the masses the search used): matched intensity
  fraction, median and median absolute deviation of the ppm errors (pairs
  at an isotope offset left out), |median|, isotope-offset share,
  cleavage-site coverage, the longest run of covered sites, complementary
  sites. One more: |RT - predicted RT|, the prediction a ridge regression
  (alpha 1) on the amino acid composition, log10 length and log10 mass,
  fitted in each fold on the training folds' targets at 1% FDR by E-value.
  The console prints, per search, the PrSMs, decoys and the targets at 1%
  FDR by E-value and rescored, and the 8 largest standardized weights.
- **The FDR caveat** (toppic_test report 0151's follow-up): with this
  order the cluster protein correction has already given some spectra the
  proteoform of their cluster's best PrSM (on 480_F3_02 331 records, 140 of
  them decoys before the correction, 192 at E = 1e300), so the decoy count
  does not check those records, and a rescored PrSM among them is accepted
  on the cluster's identification (feature-level propagation) rather than
  on its own spectrum. Rescoring the PrSMs before the correction keeps the
  decoys valid (+2.4% PrSMs instead of +6.3% in the Python prototype).
- **Two rounds of FDRs** (branch `rescore_fdr_fallback`, toppic_test
  experiment 153, report 0153; from main 2e97163f; merged into main on
  2026-10-10 with `rescore_evalue_floor` as c122997d, on main e5c15398).
  `prsm_fdr::process` with `use_rescore` computes the three FDRs by the
  E-value (on copies of the PrSMs, sorted from the input order as the
  E-value computation sorts them, so the E-value round equals a
  computation without `--rescore`) and by the rescore (on the PrSMs, as
  before), counts the targets within each level's cutoff (`FdrCutoffs`:
  the `-t` / `-T` / `-y` value when the type is FDR, else 1%;
  `getFdrCutoffs` in `toppic_process`) and keeps a level's rescored FDRs
  unless the E-value's give more targets (a tie keeps the rescore; each
  search chooses on its own; the console prints "FDR computation -
  <tool>: targets within the cutoffs by E-value / rescore: ..."). Only the
  PrSMs whose E-value FDR differs are rewritten, so where the rescore wins
  at every level the files are those of a single rescored round. The
  hidden `--disable-rescore-fallback` (argument `rescoreFallback`) keeps
  the rescored FDRs at every level; the parameter block's FDR line ends
  with ", the E-value at a level where it gives more" when the rule is on.
  It is used by `TopPIC_postprocess`, `TopPIC_multi_postprocess` and the
  second precursor merge. Results: 480_F3_02 setting d keeps the rescore at
  every level (3,573 / 2,041 / 447 against 3,374 / 2,004 / 440 by E-value)
  and writes main's record, XML and table files but the parameter line;
  the Bruker run of experiment 152 (b_p1) 1,276 / 82 / 71 (PrSMs by the
  rescore, proteoforms and proteins by the E-value; main 1,276 / 64 / 59,
  `--disable-rescore` 806 / 83 / 72); the 500-spectrum subset of 480 with
  `--mass-calib` (6 decoys) 474 / 345 / 152 by the E-value at every level
  (main 26 / 22 / 12). Without the rescoring the files are main's.
- **The rescore's PrSM-level FDR is not valid on the Bruker run**
  (toppic_test experiment 153, entrapment search: human SwissProt plus 0.5
  x shuffled `ENTRAP_` proteins): at 1% FDR the rescore accepts 1,357
  PrSMs with 59 entrapment PrSMs (estimated FDP 13.0%), the E-value 777
  with 3 (1.2%); 47 of the 59 have E = 1e300 (records the cluster protein
  correction gave a cluster's proteoform and that then failed the E-value),
  and so do 465 of the 1,357. The model clips -log10 E at -5, so such a
  record is ranked by its other features (feature intensity and the like,
  shared with its cluster); a decoy of this kind (rescore 1.06) is also
  what stops the rescored proteoform list after 64 targets. The two rounds
  fix the proteoform and protein levels but cannot see this at the PrSM
  level, where the rescore gives more targets. On 480_F3_02 the rescored
  FDP was 0.84% (report 0151).
- **The lowest rescore for the PrSMs with the maximum E-value** (branch
  `rescore_evalue_floor`, toppic_test experiment 155, report 0155; from
  `rescore_fdr_fallback` e10d3fd0; merged into main on 2026-10-10 as
  c122997d, on main e5c15398: the merge check of experiment 155,
  amendment 2). `prsm_rescore::process`
  takes an `EValueFloor` (`getRescoreEValueFloor` in `toppic_process`, for
  the three callers): the PrSMs with E >= `ExpectedValue::getMaxDouble()`
  (1e300) are left out of the training and get the lowest double,
  `std::numeric_limits<double>::lowest()` (-DBL_MAX; the user's choice,
  amendment 1, fc2e8b46; f3d913c6 gave 1 below the lowest rescore of the
  search's other PrSMs, the same order), so the FDRs rank them last, as the
  E-value does; a search without a model keeps the -log10 E-value. Every
  use of a rescore is a comparison or a negation (`rescoreKey`), so the
  value is safe; `Prsm::appendXml` writes a rescore of magnitude 1e100 or
  more with `max_digits10` digits (`-1.7976931348623157e+308`, read back
  exactly by `std::stod`), not with `str_util::toString`, whose fixed
  notation writes 309 digits. The console prints "<n> PrSMs with the
  maximum E-value (<d> decoys) get the rescore -DBL_MAX (the lowest double)
  and are left out of the training", and the parameter
  block's FDR line "(linear SVM, 3-fold cross-validation; the maximum
  E-value last)". Hidden `--disable-rescore-evalue-floor` (argument
  `rescoreEvalueFloor`) turns it off (the files of `rescore_fdr_fallback`,
  checked on the Bruker run and the 500-spectrum subset); research
  switches in `TOPPIC_RESCORE_PARA`: `floor_train=1` (the model trained on
  every PrSM as before, only their scores replaced) and `floor_value=<v>`.
  With -DBL_MAX every output equals f3d913c6's but the `<rescore>` values of
  those PrSMs (toppic_test experiment 155, amendment 1).
  Not 0: the rescore's scale puts the training folds' 1% threshold at 0 and
  the median decoy at -1, so a 0 sits at the acceptance threshold and,
  with the tie rule of the FDR (tied decoys first), stops every list there
  (the Bruker run's 2% and 5% lists equal its 1% list). Results at 1% FDR,
  entrapment FDP (human + 0.5 x shuffled `ENTRAP_`) in brackets: 480_F3_02
  setting d 3,481 / 2,046 / 443 (main's rescore 3,573 / 2,041 / 447,
  E-value 3,374 / 2,004 / 440); on the entrapment database 3,475 PrSMs
  (0.9%), 2,032 proteoforms (0.9%), 442 proteins (0), E-value 3,353 (1.0%) /
  1,981 (1.1%) / 439 (2.7%), main's rescore 3,556 (0.8%) / 2,020 (0.4%) / 442
  (0): valid, gain +3.6% PrSMs, all gained PrSMs with E > 1 (the E = 1e300
  ones gone). The Bruker run: 905 PrSMs (rescore) / 82 / 71 (E-value); on
  the entrapment database the floored rescore wins at every level, 918
  PrSMs (3.9%, against 13.0% without the floor and the E-value's 1.2%), 92
  proteoforms with 4 entrapment ones, 80 proteins with 4: still not valid.
  Its accepted entrapment PrSMs have E 0.4-346 and sit on intense features
  (the model's largest weight is the feature intensity once the E = 1e300
  records are out of the training): among the records with E < 1e300, junk
  on features above 1e7 is entrapment and decoy 10 : 10, elsewhere about
  1 : 4 (the expected 1 : 3). Trained on every PrSM (`floor_train=1`) 903
  PrSMs (5.3%), 87 proteoforms with 1; with 0 (`floor_value=0`) 895 (4.4%).
  Few decoys: st_1 (63 decoys) main's rescore keeps 16 PrSMs against the
  E-value's 177, the floored rescore 64, the two rounds the E-value's 177 /
  56 / 35; the 500-spectrum subset of 480 the E-value at every level
  (with `--mass-calib` 481 / 343 / 151, the PrSMs by the rescore). No time
  cost.
- **Results on 480_F3_02** (setting d of toppic_test report 0062, TopFD
  main 915bb387 output of experiment 151; toppic_test experiment 151,
  amendment 1): 3,573 PrSMs / 2,041 proteoforms / 447 proteins at 1% FDR
  against main 3b8041dd's 3,374 / 2,004 / 440 (the toppic_test Python
  prototype, `rescore.py` with the same features: 3,588 +- 4 / 2,045 / 443
  over 5 seeds); 242 spectra gained and 43 lost, 98 of the accepted PrSMs
  with E = 1e300 (cluster-corrected records, the caveat above). The step
  takes a few seconds (12:26 against 12:53 for the run, the same memory).
  Without the option the branch writes main's 62 files (56 equal, 6 but for
  the time lines). On the human + 0.5x shuffled entrapment
  database of experiment 151 (`ENTRAP_` proteins, indexed by this
  topindex in 3:19, 140 GB) the estimated false discovery proportion
  (entrapment PrSMs x (1 + 1 / 0.5) / accepted PrSMs) is 0.84%, 1.66% and
  2.86% at 1%, 2% and 5% FDR (3,557 PrSMs at 1%), against the E-value's
  0.98%, 1.41% and 3.56% (3,353); with 8-20 entrapment PrSMs the estimate
  has a standard error of about 0.3 percentage points at 1%.

## The second precursor search (`prsm/prsm_second_prec`, `--second-precursor`)

Branches `second_precursor` (toppic_test experiment 151, report 0151 and its
amendment 2; from main d31939d9) and `second_precursor_id` (amendment 3: the
units keyed by spectrum id and precursor id, the `Precursor ID` column; from
main 2863fa08). Off by default; toppic's hidden
`--second-precursor` (argument key `secondPrecursor`, "false"; the parameter
block prints "Second precursor search: yes" only then) searches the spectra
with two or more precursors a second time with the second precursor, and
reports up to one PrSM per precursor of a spectrum. Refused with `-c` and
with `-r` other than 1. Facts to keep straight:

- **Why.** TopFD writes every MS1 feature it assigns to an MS/MS spectrum's
  isolation window as a precursor (`PRECURSOR_MZ=a:b:c`, by intensity;
  `deconv_ms2_process.cpp` step 3; 2,440 of the 4,323 spectra of
  480_F3_02 have two or more), and the search uses the first
  (`MsHeader::getFirstPrec*`, `spectrum_set_factory`). A co-isolated
  proteoform's fragments stay unexplained. TopMPI (Wang, Xiong and Liu,
  Proteomics 2025) searches the second precursor after removing the first
  proteoform's matched masses; this is that search on TopFD's precursor
  list.
- **The passes** (`TopPIC_second_precursor` in `console/toppic_process`,
  after `TopPIC_postprocess` and before `TopPIC_report`; `TopPICProgress`'s
  identification stages are `TopPIC_identifyAll`, which both passes call).
  The first pass is a normal run up to its record file. Then
  `prsm_second_prec::writeResidualSpectra` writes `<base>_resid_ms2.msalign`
  from the spectrum file the first pass searched: the spectra with two or
  more precursors, their precursor list with the second precursor first
  (the ids of the scan's list and the feature ids kept: `MsAlignWriter`
  adds a `PRECURSOR_ID=1:0:...` line to a spectrum whose precursors are not
  at the positions of their ids, `MsAlignReader` takes the ids from it, and
  without the line numbers the precursors by position, so TopFD's files
  and every other msalign file are unchanged; the second pass's PrSMs thus
  carry `precursor_id` 1), without the masses that the
  spectrum's accepted first-pass PrSM (a target under the `-t` cutoff)
  explains: within the fragment tolerance of one of the proteoform's
  fragment masses (`theo_peak_util::geneProteoformTheoPeak`, so the shifts'
  ranges give no ions and the N-terminal label counts), as it is or one
  isotope (`mass_constant::getIsotopeMass`) lower or higher. The feature
  file is copied to `<base>_resid_ms2.feature`. The second pass runs the
  same identification stages and post-processing on that file with its own
  copy of the arguments (`rescore` "false": the merged set is rescored
  instead), so its files are `<base>_resid_ms2.toppic_*` and, post-matched,
  `<base>_resid_post_ms2.*`: its post mass matching reads and appends to
  the first pass's database, toppic's `<base>_toppic.sqlite` (the copied
  arguments carry `toppicDbFileName`), and `postBaseName` keeps `_resid` in
  the post-matched name.
- **The units.** A PrSM is the unit (spectrum id, precursor id), the
  precursor id the position of its precursor in the scan's list (0 the
  first pass's, 1 the second pass's); the spectrum ids are the scan's own
  (until 2026-10-09 the second pass's spectra were appended to the
  spectrum file with their ids raised by `SpPara::getMaxSpecNumPerFile()`,
  10,000,000, and a raised id marked a second unit). The first pass's
  record is kept as `toppic_first_record_prsm`, and `mergePrsms` writes
  the units to `toppic_second_merged`: every first-pass PrSM, and every
  second-pass PrSM whose protein (decoy prefix aside) and precursor feature
  no first-pass PrSM of the spectrum has; one that repeats a first-pass
  PrSM's protein or feature is kept in its place only with a better
  E-value (a truncated or modified form of the first proteoform is not a
  second unit: the residual fragments can still be its); in spectrum id,
  then precursor id order. The units are clustered as `TopPIC_postprocess`
  clusters (`prsm_feature_cluster`, or `prsm_simple_cluster` with `-x`;
  `toppic_second_cluster`, `toppic_second_prot`), without a cluster
  protein correction, recount or post mass matching (each pass did its
  own), rescored (with `--rescore`;
  `toppic_second_rescore`) and given one FDR over all units
  (`toppic_second_fdr`; `--fdr-groups` applies), which becomes the file's
  `toppic_record_prsm`. `TopPIC_report` then runs as usual.
- **The spectra of a unit.** `prsm_reader_util::readPrsmsWithSpectra`
  attaches the spectra of a PrSM's spectrum id, every PrSM of the id from
  the one spectrum file; with its `second_sp_file_name` (the residual
  post-matched spectra, `<base>_resid_post_ms2.msalign`, whose first
  precursor is the scan's second) a PrSM gets, from either file, the
  spectrum whose first precursor has its precursor id, as
  `prsm_util::isMatchMs` matches. The rescoring of the merged set reads
  with the residual file, so a second unit's fragment features are
  computed on the masses its search used, as its matched numbers were; the
  report's tables read the one output file, so a second unit is shown on
  the scan: `# Masses` is the scan's, the matched numbers the search's,
  and `Charge` (`MsHeader::getPrecPtrById`, the first precursor when none
  has the id) is the unit's precursor's. The output spectrum file is the
  first pass's post-matched file with the masses the second pass's post
  mass matching added (`mergePostMatchedMasses`: a residual mass the scan
  does not have within the strict fragment tolerance, appended with the
  next peak ids; not with `-D`), one spectrum per scan with the scan's
  precursor list, so its spectrum ids are TopFD's and the SQLite `prsm`
  rows join `ms2_spectrum`.
- **The outputs name the unit.** `Precursor ID` after `Spectrum ID` in the
  PrSM, proteoform and protein tables of `PrsmMatchTableWriter` (toppic
  and topmg; 0 in every run without the option; the rows of a PrSM's other
  protein matches have an empty field), `precursor_id` after `spectrum_id`
  in the SQLite tables `prsm`, `proteoform` and `protein`, and
  `<precursor_id>` in the PrSM XML, which it always had.
- **Cleanup.** `cleanToppicDir(..., copy_outputs = false)` deletes the
  residual pass's intermediates (no XML outputs to copy), and without `-k`
  the residual spectrum, feature and index files too; the merge's files are
  in its lists. The residual pass's record file stays, as record files do.
- **Not supported:** combined runs (`-c`: the merge works per file and the
  fractions' ids would collide), several spectra per set; `--topmg` modes
  are run in both passes but were not tested.
- **The FDR caveat.** The decoys of the merged set are each pass's own, so
  the target-decoy count checks a second unit against the second pass's
  decoys on the same residual spectrum. What it cannot check is a second
  PrSM explained by the first proteoform's fragments; the mass removal is
  what limits that (toppic_test report 0151 section 4.2: 1% of the second
  PrSMs share most of their matched masses with the first against 23%
  without the removal). A spectrum's two units are not independent
  measurements of two proteoforms when the precursors are one species'
  isotopes or adducts; TopFD's features separate those by mass.
- **Results on 480_F3_02** (setting d of toppic_test report 0062, TopFD main
  915bb387 output; toppic_test experiment 151, amendment 2): 2,440 spectra get a second
  precursor; 61,231 masses of 1,934 identified spectra are removed; the
  second pass finds 2,357 PrSMs (targets and decoys), of which 2,103 join
  the first pass's 4,076 units, 40 replace a first-pass PrSM of the same
  protein or feature and 214 are dropped for one. 4,100 PrSMs / 2,382
  proteoforms / 529 proteins at 1% FDR against main d31939d9's 3,374 /
  2,004 / 440 (825 PrSMs of the second pass), the prototype's combination
  of two runs (toppic_test report 0151, `resid_prec.py` + `combine.py`)
  giving 4,107 / 2,394 / 528 with 4,098 accepted units in common; with
  `--rescore` (the first pass rescored too, so the residual spectra follow
  the rescored acceptance) 4,507 / 2,559 / 565 (prototype 4,517 / 2,565 /
  554). 25.5 min against 12.5 (the second pass is a full search of 2,440
  spectra), 7-8.6 GB. Without the option the 62 output files equal main's
  (6 but for the time lines). On the human + 0.5x shuffled entrapment
  database of experiment 151 (one run with both options, 35.7 min, 10.6 GB)
  the estimated false discovery proportion (entrapment PrSMs x (1 + 1 / 0.5)
  / accepted PrSMs) by E-value is 0.96%, 2.01% and 4.92% at 1%, 2% and 5%
  FDR (4,044 PrSMs at 1%, 786 of the second pass, 8 of the 13 entrapment
  PrSMs of the second pass), and rescored 0.40%, 1.26% and 3.72% (4,459),
  the prototype's 0.96% and 0.80%; with 6-28 entrapment PrSMs the estimate
  has a standard error of about 0.3 percentage points at 1%.
- **Results with the units keyed by (spectrum id, precursor id)** (branch
  `second_precursor_id`, toppic_test experiment 151, amendment 3; the same
  TopFD output): the same 4,100 / 2,382 / 529 accepted units with
  `--disable-rescore` (825 of the second pass), rescored 4,519 / 2,563 /
  568 against the raised-id run's 4,507 / 2,559 / 565 (4,493 units in
  common; the folds are drawn over the spectrum ids, which a scan's two
  units now share); 8,274 masses of the second pass's post mass matching
  added to 1,831 scans; 25-26 min, 8.0-8.3 GB. Without the option the run
  writes main 2863fa08's files but the `Precursor ID` column: on the first
  500 spectra 57 of 63 files are identical, the 6 tables equal without the
  column and the time lines, the SQLite tables equal without
  `precursor_id` (0 in every row).

## Two-pass m/z recalibration (`ms/calib`, `topcalib`, `toppic --mass-calib`)

Branch `mass_calib` (toppic_test experiment 152, report 0152 and its
amendments 2 and 3; from main 2863fa08, merged into main on 2026-10-09 as
140d9479: on the 500-spectrum subset of 480_F3_02 in setting d the merged
main writes main 6949cad7's files without the option and the branch's model,
spectra and PrSM files with `--mass-calib`). `topcalib` is the C++
form of the experiment's `calib.py` and `recal_msalign.py`
(`tests/202610/scripts/two_pass_recal/` of toppic_test; the manual is
`docs/manual/topcalib_manual.md`), and toppic's hidden `--mass-calib` runs
the method inside a search. Facts to keep straight:

- **The method.** A zero-shift search (`toppic -s 0`, `-D` is enough) gives
  the calibrants: PrSMs with E-value <= 1e-5, no unexpected shift, no
  variable PTM (`mass_calib_fit::readCalibrants`). Precursor errors
  (precursor mass - proteoform mass, ppm of the proteoform mass) of the
  PrSMs whose difference rounds to 0 isotopes of 1.00235 Da; fragment
  errors of the spectrum's masses within 20 ppm of the proteoform's N- and
  C-terminal ion masses (the activation's ion types, the nearest one). Two
  models fitted apart by robust least squares (`fitModel`, `robustFit`: 10
  rounds of weights sqrt(min(1, 1.5 s / |r|)), s = 1.4826 median |r|,
  weighted normal equations): MS1 e = c0 + c_rt (rt - rt0) / 60 on the
  precursors, MS/MS e = c0 + c_mz (m/z - 800) / 1000 + c_rt (rt - rt0) /
  60 on the fragments; rt in minutes, rt0 the median of all calibrants.
  The two levels are fitted apart because their offsets differ: on the
  Bruker maXis run of experiments 151 and 152 MS1 -3.16 ppm, MS/MS -4.17
  ppm.
- **The shortcut** (`msalign_mass_calib::process`): the correction is
  applied to TopFD's msalign masses, not to the raw peaks: m/z' = m/z (1 -
  e 1e-6), so a neutral mass M of charge z becomes M - e 1e-6 (M + z
  proton). Precursors with the MS1 model at their m/z, fragment masses with
  the MS/MS model at (M + z proton) / z, both at the spectrum's
  `RETENTION_TIME`; six decimals; every other line copied. Experiment 152
  found it equal to recalibrating the raw peaks and rerunning TopFD at
  `-e 5` (shortcut against raw peaks: Bruker 851 / 850 PrSMs, 480 3,308 /
  3,310, the same targets at E <= 1e-5). TopFD's SQLite database is not
  changed; toppic `--mass-calib` writes a corrected copy of it (below).
- **The model file** (`MassCalibModel::write`/`read`) has calib.py's keys
  (`mz0`, `rt0_min`, `ms1` / `ms2` with `c0`, `c_mz`, `c_rt_per_hour`, the
  calibrant numbers, the robust SDs); it is written by hand with `%.17g`
  (Boost's property tree writes numbers as strings) and read with Boost's
  JSON parser, so either tool reads the other's file. `apply` needs only
  `mz0`, `rt0_min`, `ms1`, `ms2`.
- **Masses.** topcalib takes the residue, PTM and ion type masses from
  `res/base_data` (AME2020 masses) and the proton from `mass_constant`;
  calib.py has its own table (M 4.4e-7 Da lighter, from the sulphur mass).
  The fragment errors therefore differ by about 1e-4 ppm, and so do the
  MS/MS coefficients (480: c0 0.47366 against 0.47349 ppm); with
  TopPIC's masses put into calib.py the two fits agree to 1e-14. The
  proteoform string is parsed for single-residue fixed modifications
  (`(C)[Carbamidomethyl]`, found with `PtmBase::containsAbbrName`) and an
  `[Acetyl]-` prefix; any other PrSM gives no fragment calibrants.
- **Results** (toppic_test experiment 152, amendment 2; 9b4adfda):
  `topcalib apply` with calib.py's models writes recal_msalign.py's files
  byte for byte (480_F3_02 4,323 spectra, Bruker 4,337); `topcalib fit` on
  the zero-shift searches has calib.py's calibrant numbers (480 1,252
  precursors, 51,894 fragments of 1,381 PrSMs; Bruker 224 / 3,657 / 238)
  and MS1 coefficients to 1e-14; the corrected files of the C++ chain
  differ from the Python chain's by at most 5e-6 Da (0.006 ppm) in 40% /
  53% of the fragment lines, precursor lines identical. toppic main
  915bb387 at `-e 5` on the C++-corrected files: 480_F3_02 setting d 3,493
  PrSMs / 2,090 proteoforms / 448 proteins against the Python chain's
  3,492 / 2,089 / 448, Bruker 848 / 115 / 79 against 851 / 117 / 79 (the
  same targets at E <= 1e-5; the differences are spectra near the cutoff).
  fit 0.24 s (calib.py 1.7 s), apply 0.7 s (recal_msalign.py 1.1 s).
- **Built into toppic: the hidden `--mass-calib`** (argument `massCalib`;
  toppic_test experiment 152, amendment 3). `TopPICProgress_multi_file`
  runs `TopPIC_mass_calib` for every file, after the `-A` conversion and
  before the searches:
    - **The calibration search** is `TopPIC_identify` on a copy of the
      arguments: `shiftNumber` 0, no variable PTM file, no tag filter,
      TopMG disabled, no deamidation selector, and precursor and fragment
      tolerances of at least `min_tole` (10 ppm), so that an offset of a few
      ppm does not cost the calibrants their fragments and the index of the
      default tolerance serves (the index file names carry the fragment
      tolerance as the arguments write it, "10"). Its parameter block is
      printed as for any search; its files are deleted by `cleanToppicDir`
      with `copy_outputs` false (kept with `-k`).
    - **The calibrants** (`prsm/prsm_mass_calib::readCalibrants`, on
      `toppic_pre_raw_prsm` read with its spectra): target PrSMs (no
      `DECOY_`) with E-value <= 1e-5, no unexpected shift and no variable
      PTM; the precursor error of `getOriPrecMass()` against the proteoform
      mass, the fragment errors of the deconvoluted masses against
      `theo_peak_util::geneProteoformTheoPeak` (the activation's ion types,
      `SpPara::getMinMass`, the N-terminal label), both through
      `mass_calib_fit::addPrecCalibrant` / `addFragCalibrants`, which
      topcalib's table reader uses too (the `ms` layer cannot include
      `prsm`, hence the split). With fewer than 3 precursor or 4 fragment
      calibrants a warning is printed and the file is searched uncorrected.
    - **The corrected spectra**: the model `<base>_mass_calib.json`, the
      spectra `<base>_calib_ms2.msalign` (`msalign_mass_calib::process`) and
      a copy of the feature file; the file list then names the corrected
      file (with `-C` the earlier run's, when it exists). Both stay after the
      run, as `-A`'s `_mono` files do.
    - **The corrected database** (`ms/calib/sqlite_mass_calib::process`,
      toppic_test experiment 154; until then the post mass matching
      corrected the peaks it read in memory and TopFD's database stayed
      uncorrected, while the envelopes it appended were at the corrected
      m/z): toppic's database `<base>_toppic.sqlite` (experiment 156; in
      experiment 154 `<base>_calib.sqlite`, and `<base>_mono_calib.sqlite`
      from `<base>_mono.sqlite` after `-A`), written after the corrected
      spectra from TopFD's `<base>.sqlite` (also after `-A`), which is only
      read (the md5 stays) and which
      the user chose not to correct in place, so that nothing has to be
      undone for a later run without the option. The copy goes through
      SQLite's backup API (`sql/sql_db_copy`, shared with `-A`'s copy;
      TopFD writes its database in WAL mode, and pages
      may still be in the `-wal` file), toppic's own tables of an earlier
      run are dropped (no VACUUM: it may renumber the rows of the tables
      without an INTEGER PRIMARY KEY), and every m/z and deconvoluted mass
      is corrected with the model at its own spectrum's retention time
      (the spectrum table's `retention_time` in seconds / 60; SQL functions
      over `MassCalibModel::correctMz` / `correctMass`): the MS1 model for
      `ms1_peak.mz`, `ms1_env.mono_mass` / `ref_mass`, `ms1_env_peak.mz` and
      the 3D tables of `topfd -D` (`PEAKS<k>.MZ` at `RETENTIONTIME` / 60,000,
      each `CONFIG` row's `MZMIN` / `MZMAX` the range of its corrected layer;
      the k-th row is layer k), the MS/MS model for `ms2_peak.mz`,
      `ms2_env.mono_mass` / `ref_mass` and `ms2_env_peak.mz`. Not changed
      (the user's decision): the isolation windows (`ms2_spectrum.target_mz`,
      `begin_mz`, `end_mz`, the instrument's settings), nor intensities,
      retention times, ids and scores. A table `mass_calib` holds the model
      (`MassCalibModel::toJson`, the text of the model file), the source
      database (TopFD's), the version and the time. The file is written as
      `<base>_toppic.sqlite.tmp` with no journal and renamed when complete,
      and every run writes it again from TopFD's database. The envelopes that
      toppic versions before experiment 156 appended to TopFD's database
      cannot be told apart from TopFD's own (`ms2_spectrum` has no envelope
      count) and are copied and corrected with them. Without TopFD's
      database (TopFD `-N`, possible with `-D`) a note is printed and there
      is no corrected database.
    - **The search uses it.** The post mass matching
      (`prsm_post_mass_match::process`) reads the corrected database's peaks
      as they are and appends its envelopes there, and `prsm_sql_output`
      writes the PrSM tables there.
    - **Names**: the corrected spectra are `<base>_calib_ms2.msalign`
      (after `-A` `<base>_mono_calib_ms2.msalign`; their residual spectra of
      `--second-precursor` `..._calib_resid_ms2.msalign`), and
      `postBaseName` gives `<base>_calib_post_ms2`
      (`<base>_calib_resid_post_ms2`, `<base>_mono_calib_post_ms2`), so the
      result files of a run with the option are `<base>_calib_post_ms2_*`
      (with `-D` `<base>_calib_ms2_*`) beside those of a run without it
      (`<base>_post_ms2_*`) in the same directory; their database is the
      one `<base>_toppic.sqlite` of the last run (experiment 156; in
      experiment 154 `<base>_calib.sqlite` and `<base>_mono_calib.sqlite`).
      `-c` takes the fractions' record names from `postBaseName`. A combined run merges the fractions' MS1 spectra
      from TopFD's `<frac>_ms1.msalign` (`prsm_post_mass_match::
      topfdBaseName`, the spectrum file's base without `_ms2`, `_resid`,
      `_calib` and `_mono`; `msalign_frac_merge::mergeFractions`'
      `ms1_base_names`): the merge took them from the post-matched spectra's
      base, which keeps `_calib` now, and with `-D` main already looked for
      `<frac>_calib_ms1.msalign` (`-c -D --mass-calib` stopped there) and
      `<frac>_mono_ms1.msalign` (`-A`). The parameter block prints "Mass
      calibration: two-pass (zero-shift calibration search)" only with the
      option.
    - **Research switches** in `TOPPIC_MASS_CALIB_PARA` ("name=value,..."):
      `max_evalue` (1e-5), `tol_ppm` (20), `min_tole` (10). (`post_match`,
      the in-memory correction's switch, went with it in experiment 154.)
    - **Results of the corrected database** (toppic_test experiment 154,
      report 0154; branch `mass_calib_sqlite` 3affd343 and 52a8b87a against
      main 506f18b4, merged into main on 2026-10-10 as fdef25ed with
      amendment 1; setting d with `--mass-calib -e 5 --disable-rescore`):
      every corrected value equals a Python recomputation from TopFD's
      database and the run's model bit for bit (480_F3_02 27.7 M values in
      26.2 M rows, Bruker 18.3 M, a `topfd -D -A` database with its 10 layer
      tables 57.3 M), every
      other value and TopFD's schema unchanged, TopFD's database's md5
      unchanged by the run; the corrected `ms2_env` masses are the corrected
      msalign masses within 1.2e-6 Da (six decimals, the msalign file's
      rounded retention time). 480_F3_02 3,488 / 2,096 / 452 and Bruker 852 /
      116 / 79 PrSMs / proteoforms / proteins, as main; the same tables, XML
      files and SQLite PrSM tables (main's in TopFD's database), the same
      22,342 / 11,533 post-matched masses, of which 34 / 142 differ from
      main's by 1e-6 Da (main corrected the peaks with the msalign file's
      rounded retention time). A rerun writes the same files and database
      tables; a run without the option in the same folder writes main's
      files beside them, and a later `--mass-calib` run copies the
      envelopes that run appended to TopFD's database (23,585 / 21,978) into
      its corrected database without changing its files. Writing the
      corrected database: 3-5 s (480: 875 MB, 26 M rows; Bruker MS/MS only),
      26 s with topfd `-D`'s layer tables (2.5 GB); run time and peak memory
      as main's. Without the option main's files (500-spectrum subset).
    - Without the option the run writes main 2863fa08's files (a
      500-spectrum subset of 480_F3_02, setting d: every file but the
      tables' time lines, every SQLite table).
    - **Results** (toppic_test experiment 152, amendment 3; e69a52c4,
      `-u 16`): 480_F3_02 setting d with `--mass-calib -e 5
      --disable-rescore` 3,492 PrSMs / 2,089 proteoforms / 448 proteins,
      as the Python chain of experiment 152 (topcalib's corrected file
      3,493 / 2,090 / 448, no option at 10 ppm 3,374 / 2,004 / 440); with the
      default rescoring 3,667 / 2,103 / 462 against 3,573 / 2,041 / 447.
      Bruker (report 0151's line-filtered TopFD output, `-d -t FDR -T FDR -y
      FDR`, `--disable-rescore`) 852 / 116 / 79 against topcalib's 848 / 115
      / 79 and 806 / 83 / 72 without the option. The built-in models equal
      topcalib's to 2e-3 ppm in c0 (the same 1,252 / 224 precursor
      calibrants; 51,905 / 3,661 fragment calibrants against 51,894 /
      3,657, the difference not examined; the Bruker MS/MS slope 0.315
      against 0.333 ppm per 1,000 m/z). The post mass matching on corrected
      peaks adds 11,533 masses on Bruker (offset -4 ppm) against 10,186 on
      uncorrected ones, and 503 reported PrSMs gain matched fragments (203
      lose, +934 in all; the identifications are equal); on 480 (+0.5 ppm)
      22,388 against 22,412, no change. The calibration search takes 1.6
      min (480) and 2.6 min (Bruker); with `-e 5` the 480 run takes as long
      as without the option at 10 ppm (11.3 against 11.2 min).
