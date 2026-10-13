# TopPIC notes: after the search

Facts to keep straight when changing what toppic does after the
search: the order of the post-search stages, post mass matching and
toppic's database `<base>_toppic.sqlite`, combined runs (`-c`), the PTM
localization, the PrSM tables, ProForma strings and the proteoform
classification (topmg writes them too), average-mass input, the second
precursor search and the two-pass m/z recalibration (also topcalib)
(`src/prsm`, `src/stat/local*`, `src/ms/calib`, `console/toppic_process`).
One of three TopPIC files (the others: `toppic_search.md` and `toppic_evalue.md`), split on 2026-10-12
from `toppic.md`, which had moved out of the repository's `CLAUDE.md`;
the code-wide conventions and the source layout are there.

## Post mass matching (`prsm/prsm_post_mass_match`, toppic only)

After the search, the clustering, the cluster protein correction and the
match recount, `console/toppic_process` runs
`prsm_post_mass_match::process`, which matches the theoretical fragment
masses of each PrSM that no deconvoluted mass matched against the
**centroided MS/MS peaks of TopFD's SQLite database, read only** (since
toppic_test experiment 168 in every mode; from toppic's copy
`<base>_toppic.sqlite` before it; MSPathFinderT's envelope test), scores
the added masses with EnvCNN, appends their envelopes to toppic's database
`<base>_toppic.sqlite` (not with `--no-sql`), and writes
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
  `-z`, possible with `-D`) nothing is written and the PrSM tables are not
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
- **`--no-sql` writes no SQLite file** (branch `no_sql_output`, merged into
  main on 2026-10-12 as 4cdaa7a1, whose tree is the tested 94c36445;
  toppic_test experiment 164; visible in toppic and topmg; `-z` in toppic,
  topmg and topfd since branch `no_sql_short_option` (2026-10-13, experiment 164
  amendment 2; merged into main with `topdia_no_sql` as daa9098e; topfd's was `-N`, which is `--n-terminal-label` in toppic,
  and topfd_gui passes `-z`)). toppic
  (`arguments["noSql"]`): `TopPICProgress_multi_file` still deletes an
  earlier `<base>_toppic.sqlite` but writes none (`writeToppicDb`
  skipped); `TopPIC_mass_calib` writes no corrected copy and passes the
  model file instead (`arguments["massCalibModelFileName"]`, the
  `<base>_mass_calib.json` it writes with %.17g; "" without a model, e.g.
  too few calibrants), per file through `calib_model_lst`; the post mass
  matching calls `prsm_post_mass_match::processReadOnly` on TopFD's
  database (`arguments["topfdDbFileName"]`, now set for every file before
  `TopPICProgress`), and `TopPIC_report` skips `prsm_sql_output::write`
  ("The PrSM tables are not written to an SQLite database (--no-sql).").
  `processReadOnly` opens the database as `file:<path>?immutable=1`
  (`SQLITE_OPEN_READONLY | SQLITE_OPEN_URI`; no locks, no `-wal`/`-shm`
  files next to TopFD's WAL database, which `mode=ro` would leave; the file
  must not change during the run), and `LoadedPeakReader` reads the peaks
  of the PrSMs' scans in one pass over `ms2_peak` (the spectrum of a scan is
  the first in id order, as `RawPeakReader`'s query; the peaks are
  stable-sorted by m/z), because the per-spectrum query needs the index on
  `ms2_peak.spec_id` that TopFD does not write and a read-only connection
  cannot create; with a model, every m/z is corrected with
  `MassCalibModel::correctMz(2, mz, retention_time / 60)`, as
  `sqlite_mass_calib` corrects the copy. No envelope is written. The
  parameter block prints "SQLite database output: No" only when the option
  is set. topmg (`topmg_process`) skips its `prsm_sql_output::write` into
  TopFD's database. The GUIs have no control for it. Checked (report 0164)
  on st_1 against main 01d2db5a: without the option every result file and
  `<base>_toppic.sqlite` (all tables; `mass_calib` but its time) equal
  main's for the default, `-D`, `--mass-calib`, `--second-precursor`, `-c`
  and `-c -C`, topmg (`-u 1`) too, with TopFD's database after topmg equal;
  with it every result file equals the run without it (apart from the
  parameter line), the folder holds TopFD's database alone with its md5
  unchanged (no copy, no -wal / -shm), an earlier copy is deleted, and
  without TopFD's database the run stops with the usual error (runs with
  `-D`). 480_F3_02 setting a: the 15 result files of main (3,207 / 1,918 /
  436), 947 MB of `<base>_toppic.sqlite` not written; post mass matching
  13.1 s against 16.4 s, and the 5.4 s table write and the 0.9 s copy
  gone.
- **toppic's database holds TopFD's MS/MS envelope tables only**
  (branch `toppic_db_ms2_env` from main 04b96954, toppic_test experiment
  168, the user's request: "copy only ms2 envelope and ms2 envelope peaks
  ... in default setting. Copy all peaks only if recalibration is turned
  on. For average mass, copy only ms2 envelope ..."). `writeToppicDb`
  writes `<base>_toppic.sqlite` with TopFD's `ms2_env` and `ms2_env_peak`
  (`sql_db_copy::copyTables`: a new file under the `.tmp` name, TopFD's
  database attached read only as immutable (`sql_util::immutableUri`, moved
  there from the post mass matching), each table created with TopFD's own
  CREATE statement from `src.sqlite_master` and filled with `INSERT ...
  SELECT ... ORDER BY rowid`, closed in WAL mode and renamed as
  `closeCopy`); `sql_db_copy::copy` (the backup API copy of every table)
  is gone. The post mass matching appends its envelopes to it and the
  report writes the PrSM tables into it, as before. Only `--mass-calib`
  with a model writes every table (`sqlite_mass_calib`'s corrected copy,
  unchanged); with too few calibrants the database is the envelope copy
  (the peaks are TopFD's). `-A` gets the envelope copy (TopFD's database
  holds monoisotopic masses), `-A --mass-calib` the corrected copy; `-D`,
  `--second-precursor` and the fractions of `-c` the envelope copy; `-C`
  and `--no-sql` as before. A reader of toppic's database takes every
  table it lacks (the spectra, the MS1 and MS/MS peaks, the MS1 envelopes,
  `CONFIG` and the 3D levels, `ms_info`) from TopFD's database: TopMSV
  (topmsv_private branch `two_sqlite_files`) opens toppic's file and
  attaches TopFD's, so an unqualified table name resolves to toppic's file
  first. The post mass matching (`prsm_post_mass_match::process`, one
  function since then: `processReadOnly` and `RawPeakReader`, which
  created the index `ms2_peak_spec_id_index` in the copy, are gone) reads
  the peaks of the PrSMs' scans from TopFD's database in every mode, read
  only as immutable in one pass (`LoadedPeakReader`), corrected in memory
  with the model after `--mass-calib` (`TopPIC_mass_calib` passes the
  model file with or without `--no-sql`), and then opens toppic's database
  (its third argument, "" with `--no-sql`) to append the envelopes.
- **A missing database is an error, not a skip.** TopFD's database of a
  spectrum file is `prsm_post_mass_match::sqlFileName(sp_file)`
  (`<base>.sqlite` for `<base>_ms2.msalign`; use it, do not re-derive the
  name). It is checked in three places with the same message ("run TopFD
  without disabling its SQLite database output (without -z) ..., or run
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
Its database is toppic's `<base>_toppic.sqlite`, with the same tables of
TopFD's `<base>.sqlite` as without the option (post mass matching section):
TopFD's masses are monoisotopic already (on toppic_test experiment 154's `topfd -D
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
  `<base>_resid_post_ms2.*`: its post mass matching reads TopFD's database
  and appends to the first pass's database, toppic's `<base>_toppic.sqlite`
  (the copied arguments carry `topfdDbFileName` and `toppicDbFileName`),
  and `postBaseName` keeps `_resid` in the post-matched name.
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
      the 3D levels of `topfd -D`: since toppic_test experiment 165 the
      `LEVEL` column of `ms1_peak`, corrected with it, and each `CONFIG`
      row's `MZMIN` / `MZMAX` set to the range of its corrected level (the
      k-th row is level k, the rows with `LEVEL` >= k; one `GROUP BY LEVEL`
      pass over `ms1_peak`); in a file of an earlier TopFD the tables
      `PEAKS<k>` (`MZ` corrected at `RETENTIONTIME` / 60,000, the k-th row
      the range of `PEAKS<k>`), the MS/MS model for `ms2_peak.mz`,
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
      database (TopFD `-z`, possible with `-D`) a note is printed and there
      is no corrected database.
    - **The search uses it.** The post mass matching
      (`prsm_post_mass_match::process`) appends its envelopes there, and
      `prsm_sql_output` writes the PrSM tables there. The post mass
      matching reads TopFD's peaks and corrects them in memory with the
      model file (since toppic_test experiment 168; up to then it read the
      corrected database's peaks, which are the same values, and with
      `--no-sql` it did as now).
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
