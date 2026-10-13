# TopFD and TopDIA notes

Facts to keep straight when changing the spectral deconvolution and the
feature detection of topfd and topdia (`src/topfd`, `src/topdia`,
`src/ms/env`, `src/ms/mzml`, `src/common/emass`, the ONNX models in
`res/`) and their GUIs. TopDIA runs the TopFD layer, so its notes are
here too.
Moved out of the repository's `CLAUDE.md` on 2026-10-12; the code-wide
conventions, the dependencies and the source layout are there.

## Quantized ONNX models (`res/envcnn_models`, `res/ecscore_models`)

The models the tools load, `envcnn_two_block.onnx` and
`ecscore_seven_attr.onnx` (file names fixed in `onnx_env_cnn::initModel` and
`onnx_ecscore::initModel`), are dynamically quantized (int8 weights,
activations quantized at run time) forms of the float PyTorch exports kept
next to them as `envcnn_two_block_full.onnx` and
`ecscore_seven_attr_full.onnx`, which no code loads (toppic_test experiment
131, issue #595). Both are made by `quantize.py` of that experiment (Python
onnxruntime 1.24.1; the files keep opset 11 and IR version 6 and load in the
vendored runtime 1.23.2). They are used by TopFD and TopDIA (both models) and
by TopPIC's post mass matching (EnvCNN).

- **Gemm -> MatMul + Add first.** The export writes every fully connected
  layer as `Gemm(A, W, b, transB = 1)`, which onnxruntime's
  `quantize_dynamic` does not quantize; the rewrite `MatMul(A, W^T) + b`
  gives the same logits.
- **EnvCNN (`fc_row`).** The three fully connected layers (fc.0 9,600 ->
  2,048 holds 19.7 M of the 21.9 M parameters) get int8 weights with one
  scale per output channel, and their inputs are quantized to uint8 with
  **one scale per row** (envelope): `max(row) / 255`, valid because each is
  a ReLU output (>= 0), built from ReduceMax, Max, Div, Round, Cast,
  MatMulInteger and Mul. `quantize_dynamic` itself uses one activation scale
  per call: the deconvolution scores up to 256 envelopes per call, and its
  input has a few rows with values up to 825 (experimental intensities
  divided by the theoretical envelope's highest peak; a neighbouring intense
  peak), which coarsened the whole batch: scores moved by 0.06 on average,
  4.7% of the fragment masses crossed toppic's 0.4 (`-F`), and 142 PrSMs of
  480_F3_02 changed their proteoform. With one scale per row a score does not
  depend on the batch (mean change 0.002). The convolutions stay float:
  ConvInteger was slower than the float convolutions and less accurate. The
  file is 22 MB instead of 87 MB.
- **ECScore (`fc`).** `quantize_dynamic` with int8 weights, one scale per
  tensor (113 k parameters, 122 kB). TopFD calls it once per candidate
  feature (about 1 s per run in total either way); its scores change by
  0.027 on average (its seven inputs have different ranges but share one
  activation scale).
- **Results** on 480_F3_02 (toppic_test experiment 131, report 0131; main
  7eb5dbe9 with the models swapped in `res/`, i7-11700 with AVX-512 VNNI,
  `topfd -u 16`, three runs each): TopFD 283 s against 408 s with the float
  models (MS1 deconvolution 144 against 216 s, feature detection 64 against
  75 s, MS/MS deconvolution 74 against 115 s; EnvCNN was 84-88% of the
  deconvolution time and its int8 calls are 1.6-1.8x faster). 42,447
  features against 45,360 (the quantized ECScore; the EnvCNN alone gives
  45,367), 320 of the 4,318 MS/MS precursors changed. toppic (setting d of
  toppic_test report 0062): 3,346-3,347 PrSMs, 1,996-1,998 proteoforms and
  436 proteins against 3,336-3,342 / 1,998-2,000 / 437; 33 PrSMs with
  another proteoform than in a float run (1-3 between float runs). The int8
  speed-up depends on the CPU's int8 dot-product instructions; it was not
  timed on CPUs without AVX-512 VNNI or on ARM.
- **In the tree** (branch `quantized_onnx_models` on main 72534ded, merged
  into main on 2026-10-06; toppic_test experiment 131, amendment 1; one run
  each): TopFD 278.6 s
  against 402.8 s with the float `*_full` models; its outputs are as close to
  the experiment's arm B as two runs of B are to each other. toppic: 3,347
  PrSMs, 1,997 proteoforms, 437 proteins against main's 3,338 / 2,001 / 436.
  TopPIC's post mass matching with the int8 EnvCNN (5.7 s of EnvCNN against
  16.9 s) gives the same tables and XML files as with the float model on the
  same TopFD output; only 3,772 of the 724,767 EnvCNN scores in
  `_post_ms2.msalign` change (up to 0.12, 41 across 0.4). TopDIA was not
  tested.
- **The mass-difference input of unmatched peaks** (`onnx_env_cnn::
  extractFeature`, and ECScore's EnvCNN feature in `comp_env_cnn_score`): a
  theoretical peak without an experimental peak within 0.02 m/z gets -0.02
  or +0.02, alternating over the peaks (even index -0.02). It was a random
  sign from `std::rand() % 2`, one stream for the process, so an envelope's
  score depended on how many envelopes had been built before it and, with
  several threads, on their scheduling: the "EnvCNN jitter" of two runs of
  one binary (toppic_test report 0036 named the cause). toppic_test
  experiment 139: a PTCR spectrum gave 95 or 94 species at `topfd -u 1`
  depending on whether an HCD scan was deconvoluted before it in the same
  thread (scores 0.562 against 0.561 moved an envelope across 0.2). With
  the alternating sign two `topfd -u 16` runs on 480_F3_02 write the same
  `_ms1.msalign`, `_ms2.msalign` and feature files byte for byte and the
  same SQLite rows (in another insertion order; each thread buffers its
  records); two runs of main 6eb40e97 differ in 771 of 2,803 MS1 spectra
  and 193 MS1 features. Against main the output lies within main's own
  variation (814 MS1 spectra, 131 MS1 features), and toppic in setting d
  of toppic_test report 0062 gives 3,355 PrSMs, 1,998 proteoforms, 438
  proteins, main's first run's 3,355 scans and proteins (its second run
  3,354 / 1,997 / 438). TopPIC's post mass matching uses the same input
  builder, so its EnvCNN scores are repeatable now too (not measured).

## Envelope refinement and the core-peak intensity fit (`ms/env/match_env_refine`)

After the DP stage, `DeconvSingleSp::postprocess` calls
`match_env_refine::mzRefine`, which regenerates each envelope's theoretical
distribution from its reference peak and rescales it to the experimental
peaks. It is two-stage, and the split is deliberate (issue 532):

1. **Monoisotopic choice** (current / previous / next isotope) uses the
   whole-envelope distance, unchanged from the original code. A core-only
   choice was tried and matched no more theoretical b/y fragments of horse
   myoglobin on a test MS/MS spectrum, so do not "simplify" stage 2's core
   into stage 1.
2. **Intensity ratio** of the chosen distribution is refit on the **core
   peaks only**: theoretical intensity >= `EnvPara::refine_core_ratio_`
   (0.5) x the reference peak, and present in the experimental envelope;
   fewer than 3 such peaks -> whole envelope (small envelopes unchanged).
   `compDistWithNorm`/`compDist` take an index list (empty = all peaks).
   Fitting every peak let tails inflated by overlapping neighbouring
   envelopes drag the ratio up, leaving theoretical apexes far above the
   observed ones (1.58x on the 16896 Da envelope near m/z 891 of
   `~/code/mms/data/myoglobin_peaks.txt`).

The theoretical intensities set here are what reach the outputs: the
`msalign` intensity per mass is the theoretical envelope's intensity sum,
the `.env` file's `THEO_INTE` column, and the SQLite `ms1_env`/`ms2_env`
rows. `assignIntensity` (called before `mzRefine` and again after sorting)
only splits an experimental peak among envelopes that share it; it never
changes theoretical intensities. The greedy MS1 and MS/MS searches (TopFD's
defaults) do not run this refinement. When touching any of this, validate
with `topfd -T -j --deconv-ms2 msdeconv` on the myoglobin peak list (without
the option `-T` runs the greedy search): compare the deconvoluted masses
against the sequence's b/y ions (46 exact / 14 off-by-one-isotope is the
reference), and count envelopes whose theoretical apex exceeds the observed
apex by > 20% (77 of 611 before the core fit, 1 after).

## Deconvolution methods (`--deconv-ms1`, `--deconv-ms2`)

topfd and topdia choose the deconvolution of each spectrum level with two
visible options: `--deconv-ms1 <greedy|msdeconv>` and `--deconv-ms2
<greedy|msdeconv|flashdeconv>` (`greedy`: the greedy searches of the next
sections; `msdeconv`: MS-Deconv's dynamic programming, `DeconvSingleSp`;
`flashdeconv`: the FLASHDeconv port below). Both tools also take the
short names `-p` and `-q` (2026-10-08; in topfd briefly `-x` and `-y`); in
topdia `--ms1-intensity-correlation-cutoff` moved from `-p` to `-y` for
them and `--ms2-intensity-correlation-cutoff` from `-P` to `-Y` (topdia_gui's
command emits `-y` and `-Y`). Facts to keep straight:

- **`TopfdPara` stores the names** (`ms1_deconv_method_`,
  `ms2_deconv_method_`; `getMs1DeconvMethods()` / `getMs2DeconvMethods()`
  list the accepted values, `setMs1DeconvMethod` / `setMs2DeconvMethod`
  return false for any other, which both parsers report as an error). The
  drivers ask `isMs1GreedyDeconv()`, `isMs2GreedyDeconv()` and
  `isUseFlashDeconv()`, so a new method is a new name plus its branch in
  `deconv_ms1_process::geneTask`, `deconv_ms2_process::deconvMsTwo`,
  `topfd_single_process` or the `topfd_process` / `topdia_process` choice
  of the MS/MS process class.
- **Defaults differ by tool**: `TopfdPara`'s own default, TopFD's, is
  `greedy` for both; TopDIA's `getTopfdParaPtrForTopdia` sets `msdeconv` for
  both (the greedy searches were measured on DDA data only), so a default
  TopDIA run is unchanged. TopDIA's MS/MS step (`ms2_raw.msalign`, run with
  `setMissingLevelOne(true)`) can use `FlashDeconvMs2Process`, which handles
  the missing-MS1 mode as `DeconvMs2Process` does.
- **`-T` refuses `flashdeconv`** (the FLASHDeconv process reads mzML files
  only); before, `-F` was silently ignored there.
- **Parameter block**: "MS1 deconvolution: Greedy" and "MS/MS
  deconvolution: Greedy" / "FLASHDeconv" are printed only for those methods,
  so `msdeconv` output (and TopDIA's default output) keeps the earlier block
  with the DP's ranking window and harmonic check lines.
- **Removed options** (unrecognised now): `-F`/`--ms2-flashdeconv`, the
  hidden `--ms1-greedy-deconv` / `--ms2-greedy-deconv`, and the short-lived
  `--deconv-ms1-msdeconv` / `--deconv-ms2-msdeconv` of commit b3d71867.
  The letter `-F` is free in topfd.
- **The envelope ranking score** (2026-10-12): `-n`/`--env-rank
  <envcnn|msdeconv>` (default `envcnn`; `--deconv-score` for a few hours)
  replaces the flag `-n`/`--msdeconv` (unrecognised now) in topfd and
  topdia. `TopfdPara` stores the name
  (`deconv_score_`; `getDeconvScores()` lists the values, `setDeconvScore`
  returns false for any other, which both parsers report), and
  `isSortUseMsDeconv()` is `deconv_score_ == "msdeconv"`, so the
  deconvolutions and the parameter block's "Fragment envelope ranking" line
  are unchanged. The GUIs keep their checkboxes (topfd_gui "Sort using
  MS-Deconv score", topdia_gui "Use MS-Deconv score", `msDeconvCheckBox`),
  and `gui/util/command.cpp` always emits `--env-rank <score>`. The ranking
  only decides which envelopes `--frag-num-filtering` (`-j` in topfd and
  topdia; MS/MS envelopes only: `DeconvSingleSp`, which also deconvolutes
  the MS1 spectra with `--deconv-ms1 msdeconv`, applied it to MS1 envelopes
  too until 2026-10-12) keeps; the EnvCNN cutoffs and the msalign
  score column use the EnvCNN score with either value.
- **The EnvCNN cutoffs** (2026-10-12): `-f`/`--envcnn-cutoff-ms1`
  (`TopfdPara::ms1_env_cnn_score_cutoff_`, default 0.2) and
  `-g`/`--envcnn-cutoff-ms2` (`ms2_env_cnn_score_cutoff_`, default 0; topfd's
  `-v`/`--env-cnn-cutoff` before; topdia had none) of topfd and topdia, with the greedy and MS-Deconv deconvolutions
  only. MS1: the greedy search drops its envelopes below the larger of the
  cutoff and `GreedyPara::min_envcnn` (now 0; its default 0.2 became the
  option's), `DeconvSingleSp::postprocess` those of `deconvMsOne`'s
  whole-spectrum search (`ms_level` 1; `--deconv-ms1 msdeconv` had no
  EnvCNN filter before, so its default output lost the envelopes below
  0.2); the precursor-window envelopes (`addPrecWinEnvs`, `deconvMsOne`'s
  precursor-window step) are not filtered. MS/MS as before
  (`DeconvSingleSp`, `greedy_ms2_deconv`; MS3 and `-T` too; not
  FLASHDeconv; PTCR keeps `TOPFD_PTCR_PARA`'s `min_envcnn`). The parameter
  block prints "MS1 envelope EnvCNN score cutoff" in the MS1 section.
  TopDIA takes the same defaults, so its default MS1 output (MS-Deconv)
  lost the envelopes below 0.2 too. The letters of topfd's
  `--frag-num-filtering` (`-g` -> `-v` -> `-j`; `-F` for a few hours) and
  `--disable-additional-feature-search` (`-f` -> `-d`) moved for them.
  topdia's `-d`/`--final-filtering` became `--frag-num-filtering`, and
  both tools take `-j` for it (2026-10-12; topdia's `-v` is
  `--pseudo-cutoff`, and topdia has no `-d` now); the commands of
  topfd_gui and topdia_gui emit `-j`. `topfd -h` and `topdia -h` list `-o
  -w -p -q -f -g -n -j`.
- **GUI wiring**: `topfddialog.ui` and `topdiadialog.ui` have
  `ms1DeconvComboBox` and `ms2DeconvComboBox`, whose item texts are the
  option values (as `activationComboBox`'s are), set from the dialog's
  default `TopfdPara` (TopFD: a default-constructed one; TopDIA:
  `Argument::getTopfdParaPtrForTopdia`) and read back with
  `setMs1DeconvMethod(currentText())`; `gui/util/command.cpp` always emits
  `--deconv-ms1 <m> --deconv-ms2 <m>`. In topfd_gui they are the middle
  column's rows 3 and 4, "MS1 deconvolution" and "MS2 deconvolution" (was
  "MS/MS deconvolution"; they came in as the third column's rows 4 and 5,
  the window grew by 40 pixels, and the former "Use FLASHDeconv for MS/MS"
  checkbox is gone). Since 2026-10-12 topfd_gui's parameters are three
  columns of five rows, 281 pixels wide at x = 30, 340 and 650 (the check
  box columns of "Additional settings" start at the same x): maximum
  charge, maximum mass, m/z error, "Default prec window" (was "Default
  precursor window (m/z)"), fragmentation; MS1 and MS2 S/N ratio, the two
  deconvolution drop-downs (115 pixels: "flashdeconv" needs 111), thread
  number; "MS1 EnvCNN cutoff" and "MS2 EnvCNN cutoff" (`ms1EnvCnnCutoffEdit`,
  `ms2EnvCnnCutoffEdit`, set from `TopfdPara`'s defaults 0.2 and 0, read
  back with `setMs1EnvCnnScoreCutoff` / `setMs2EnvCnnScoreCutoff`; the
  command emits `-f` and `-g`), ECScore cutoff, feature split ratio,
  feature min scan num. In topdia_gui the middle
  column's rows 3 and 4, above the MS1 and MS2 ECScore cutoffs ("Use
  single scan noise level" moved into the Additional settings box). The
  forms place widgets at fixed coordinates, so moving one means editing the
  geometries of the widgets below it. topdia_gui's window is 1103 pixels
  wide (was 965; left column rows 305, right column rows 385, deconvolution
  drop-downs at most 120; the correlation cutoff labels say "MS1 envelope"
  and "MS2 envelope", without "isotopic") so that every label and
  drop-down fits the desktop font Ubuntu Sans 11, not only Qt's offscreen
  Sans Serif 9. To check a form, build a small program on the generated
  `build/<tool>_gui_autogen/include/gui/<tool>/ui_<tool>dialog.h`, run it
  with `QT_QPA_PLATFORM=offscreen` and the font set by
  `QApplication::setFont`, and compare each `QLabel`'s
  `sizeHint().width()` with its `width()` (and `grab().save()` a
  screenshot).

## FLASHDeconv MS/MS deconvolution (`topfd/deconv/flashdeconv_ms2_process`)

`FlashDeconvMs2Process` is a second `DeconvMs2Process`: the same mzML group
reading, one task per MS/MS spectrum on the thread pool, per-thread msalign
files merged at the end and buffered SQLite records, but the spectrum is
deconvoluted with the FLASHDeconv algorithm (OpenMS `SpectralDeconvolution`,
`PeakGroup`, `PeakGroupScoring`, `FLASHHelperClasses`) instead of
`DeconvSingleSp`. `--deconv-ms2 flashdeconv` selects it in topfd and topdia
(`TopfdPara::isUseFlashDeconv`, `topfd_process.cpp`, `topdia_process.cpp`;
previous section). Facts to keep straight:

- **Everything lives in the cpp**, in the namespace `flashdeconv_ms2_process`:
  `Averagine` (the theoretical distributions), `LogMzPeak`, `PeakGroup` (one
  deconvoluted mass with its scores) and `FlashDeconvSingleSp` (one spectrum).
  Each function's comment names the OpenMS function it ports; keep that when
  editing so the two can be compared.
- **FLASHDeconv's constants, not TopFD's.** The tolerance is 10 ppm
  (`kMs2TolerancePpm`), the isotope spacing 1.002371 Da (`kIsoDaDistance`, not
  1.00335), the log m/z bin is tolerance / 2.5, the thresholds are cosine 0.85,
  charge SNR 0.25, QScore 0.2, charge score 0.7, and the QScore weights are
  the fixed ones of `PeakGroupScoring`. TopFD's `-e` (m/z error) is used only
  to match the theoretical envelopes to the peaks for the output.
- **The averagine table is EnvBase's** (`res/base_data/theo_patt.txt`, one
  entry per 25 Da from 50 Da to `-m`), prepared as `PrecalculatedAveragine`
  prepares CoarseIsotopePatternGenerator's: the weak end isotopes are trimmed
  until 0.01% of the power is gone, the intensities are L2-normalised (so
  `getCosine` needs no norm of the theoretical side), and the left / right
  isotope counts from the apex are a running maximum over the masses. It is
  built once per run in `process()` and shared read-only by the threads.
- **Output.** Every surviving peak group becomes a `MatchEnv` of its
  representative (highest per-charge SNR) charge: the EnvBase envelope of the
  monoisotopic mass, scaled so its intensity sum is the group's intensity over
  all charges, **trimmed to the noise floor exactly as `env_detect::compFullEnv`
  trims the standard pipeline's own envelope** (`getSubEnv(1.0, bound_inte,
  INT_MAX, INT_MAX)`; `bound_inte`/`min_inte` are `deconvMsTwo`'s own estimate
  of `DeconvData::initMinInte`, via `baseline_util::getBaseLine` on the
  spectrum, gated by `-estimate-min-inte`/`isEstimateMinInte()`) before
  matching to the peaks by `ExpEnv`. Without this trim every envelope carried
  the averagine table's full isotope tail (EnvBase's entries run to 100+
  peaks for large masses) into the SQLite output, 3–4x more than
  `DeconvMs2Process` writes for the same data (toppic_test experiment,
  2026-09). The QScore is stored as the EnvCNN score (what
  `match_env_util::getDeconvMsPtr` reports and sorts by) and the isotope
  cosine as the MsDeconv score; no EnvCNN model runs. Note FLASHDeconv's own
  msalign writer reports the maximum charge of the range, this path the
  representative charge.
- **Not ported:** the MS1 side and the precursor peak group search (the MS/MS
  max mass and charge come from TopFD's precursor features, as in
  `DeconvMs2Process`), decoy / FDR runs, target masses, and the tolerance
  auto-estimation.

## Greedy MS1 deconvolution (`topfd/deconv/greedy_ms1_deconv`, TopFD's default)

`deconv_ms1_process::deconvMsOneGreedy` is TopFD's MS1 deconvolution
(`TopfdPara::isMs1GreedyDeconv`, true by default; the parameter block then
prints "MS1 deconvolution: Greedy"). It replaces the precursor-window step
and `DeconvSingleSp` (candidate envelopes, filters, MS-Deconv's DP) of
`deconvMsOne` for MS1 spectra with `greedy_ms1_deconv::deconv`, a
THRASH-like subtractive search in which a peak's charge is decided with all
the charge states of its mass; the output (msalign lines per charge state,
EnvCNN scores, SQLite rows) is written by the same `writeMsOne` as
`deconvMsOne`'s. `deconvMsOne`, TopFD's default up to main 89df2a8c, runs
with `--deconv-ms1 msdeconv`, TopDIA's default (see "Deconvolution
methods" above for the option, the defaults and the GUIs). Facts to keep
straight:

- **One spectrum** (`GreedySearch::run`): the most intense peak of the
  residual spectrum (residual >= noise x MS1 S/N, the thresholds of
  `DeconvData::initMinInte`) is the seed; for every charge z with a peak one
  isotope spacing from the seed (pre-check), the EnvBase distribution's apex
  is aligned with the seed, scaled to the residual intensities of its core
  isotopes (>= 0.5 x apex, least squares), and extended to the same mass at
  z-1, z-2, ... and z+1, z+2, ... while those envelopes fit (isotope cosine
  >= 0.6, the apex matched, `min_match_peak_num_`; one failing charge in a
  row allowed). The charge whose series has the highest score wins; an
  envelope's score is the explained intensity (sum of min(observed,
  theoretical)) minus the theoretical intensity predicted where there is no
  peak, times the squared isotope cosine (`score_mode` 1). Then the
  monoisotopic offset (+-2 isotopes) is chosen by the cosine of the
  averagine distribution with the isotope intensities summed over the
  series' charges (`mono_by` 1, FLASHDeconv's
  `getIsotopeCosineAndIsoOffset`), and the mass's envelopes at all its
  charges are subtracted: a matched peak above 1.5x its theoretical
  intensity keeps the difference, others are used up, and the tail
  isotopes up to where the observed series stops or rises again are removed
  whole (`dyn_bound`).
- **Guards against artefacts**: an envelope within 3.5 Da of a mass already
  found at the same charge is not used (`dup_window`, as main's 12 Da
  `filterEnvByMz`; without it about 20% of a dense spectrum's masses have
  a same-charge near-duplicate, mostly from peaks the earlier envelope never
  matched rather than from its residual, and the feature detection merges
  them, so the identifications do not change); a seed charge z whose
  isotopes off the grid of every 2nd/3rd isotope are much weaker than those
  on it is rejected as a harmonic (`harm_ratio` 0.4); a mass found at one
  charge only needs an isotope cosine of 0.8 (`single_min_cos`); envelopes
  with an EnvCNN score below 0.2 are not written (`-f`/`--envcnn-cutoff-ms1`,
  until 2026-10-12 the switch `min_envcnn`; a third of the envelopes, which
  only seeded extra MS1 features).
- **Shared peaks** (toppic_test experiments 157 and 158, reports 0157 and
  0158; branch `greedy_shared_peaks`, merged into main on 2026-10-10 as
  2df670d5, then branch `shared_peaks_same_charge`, merged into main on
  2026-10-11 as 8269c2c2, so the same-charge rule with the default -1 is
  main's from then on). Before a mass's envelopes are added to the
  result, `GreedySearch::removeSharedEnvs`
  builds each one as it will be reported (`getMatchEnv`: the peaks within
  TopFD's m/z error with an intensity of at least the noise) and drops it
  when it shares more than k spectrum peaks (the count of
  `CoTable::cntSharedPeakNum`) with any one envelope of the same charge
  already in the spectrum's result (charge-series order, the seed charge
  first). Envelopes of other charges do not count (the same-charge rule of
  `shared_peaks_same_charge`; 2df670d5 counted every charge). Envelopes of
  one charge share peaks only when their isotope grids line up: on
  480_F3_02 without a limit the MS1 pairs of one charge are almost all 4
  or more isotope spacings apart (7,028 of 7,030; mostly 4-7 Da at charges
  4-6, a second mass on the tail of a large envelope). A dropped envelope
  is neither reported nor subtracted, as one that fails `single_min_cos`;
  a mass with none left is dropped. k is the hidden
  `--ms1-max-shared-peak-num` (`TopfdPara::getMs1MaxSharedPeakNum`;
  default -1, no limit: the test is off and the search is ecf67378's, since
  no k raised the identifications on 480_F3_02; 2df670d5 had the default
  1); the parameter block prints "MS1 maximum shared peak number (same
  charge)" under "MS1 deconvolution: Greedy" only when the limit is set,
  so that a default run writes the block of earlier versions. Not
  tested: the precursor-window envelopes of `addPrecWinEnvs` (an envelope
  of a window without one may share any number of peaks with the search's
  envelopes) and the PTCR search (`deconvSeries` passes -1). On 480_F3_02
  ecf67378's MS1 envelopes form
  124,239 pairs that share peaks (7,030 of one charge). Any charge (2df670d5,
  k = 0 / 1 / 2 at both levels): the MS1 masses fall from 845,164 to
  742,233 / 821,822 / 837,855, the MS/MS precursors from 9,228 to 8,893 /
  9,087 / 9,181, and toppic (setting d) finds 3,340 / 3,379 / 3,372 PrSMs,
  1,969 / 2,004 / 1,995 proteoforms and 433 / 445 / 440 proteins against
  ecf67378's 3,374 / 2,004 / 440. Same charge (k = 1 / 2 / 3 at both
  levels): MS1 masses 843,382 / 844,278 / 844,602, precursors 9,160 /
  9,164 / 9,231, toppic 3,381 / 3,375 / 3,373 PrSMs, 2,011 / 2,006 /
  2,004 proteoforms, 442 / 444 / 441 proteins (any charge at k = 3:
  3,370 / 1,997 / 439). TopFD is repeatable on 480_F3_02 (a second run
  writes the same msalign and feature files).
- **Precursor windows without an envelope** (postprocessing step,
  `addPrecWinEnvs`, switch `prec_win_post`, on by default; toppic_test
  experiment 113): `deconvMsOneGreedy` passes the precursor windows of the
  MS/MS spectra that follow the MS1 spectrum (`PrecWin`, not when
  `prec_win` runs the precursor-window step first). After the EnvCNN
  filter, a window in which no reported envelope has its reference (apex)
  m/z, the test the MS1 feature detection uses to seed its search for
  windows without a feature (`env_coll_detect`, `SeedEnv::getReferMz`),
  gets the envelope of its most intense peak (`GreedySearch::fitWindow`,
  on the original intensities): one peak in the window is a charge 1
  peak; with more, every charge is tried as for a seed but without the
  acceptance thresholds (the harmonic test stays) and the best series
  score wins; the monoisotopic offset is kept only when the apex stays in
  the window. One envelope per window (the window's charge), with its
  EnvCNN score but not filtered by it; a window without a peak gets
  nothing. `TOPFD_GREEDY_PREC_WIN_LOG=<file>` writes one line per window.
  On 480_F3_02 (toppic_test experiment 113, report 0113) 80 of 4,323
  windows had no envelope (40 deconvoluted at charge 1); 20 more MS/MS
  spectra got a precursor (4,318 against 4,298), none identified (their
  MS/MS spectra have a median of 4 fragment masses), and one harmonic
  precursor was corrected (3,487 Da at charge 3 -> 6,974 Da at charge 6, a
  better PTMS proteoform). PrSMs 3,331 and 3,338 against main's 3,337 and
  3,338 (the run-to-run variation); same TopFD time. A new charge 1 feature
  can cover a later window at the same m/z, so that the feature search for
  windows without a feature skips it (scans 5806 and 5815 lost a 10,943 Da
  precursor this way; neither is identified in any run).
- **Research switches**: the defaults are in `GreedyPara` (cpp); the
  environment variable `TOPFD_GREEDY_PARA` ("name=value,...") overrides
  them, so that variants can be run without a rebuild (toppic_test
  experiment 104 lists the runs). `TOPFD_GREEDY_DEBUG_SCAN=<scan>` writes
  the decisions of one MS1 scan to `greedy_scan<scan>.txt` in the working
  directory.
- **Cost**: the search itself is cheap (about 70 ms per 480_F3_02 MS1
  spectrum); almost all of the extra MS1 time against main is EnvCNN on
  the envelopes before the EnvCNN filter (`-f`; about twice main's
  number).
- **Results** (toppic_test experiment 104, report 0104): on 480_F3_02
  with main's TopPIC in the benchmark setting, 3,181-3,182 PrSMs, 1,887
  proteoforms and 425 proteins against main's 3,141-3,153 / 1,870 / 420;
  the gain is corrected precursor isotopes and charges (and fewer decoys
  on the spectra whose precursor changed). TopFD takes about the same time
  (334-336 s against 285-333 s) and memory. With both greedy searches
  (toppic_test experiment 109, rounds 2 and 3): 3,316-3,318 / 1,981 / 433.
  The search gives fewer MS/MS spectra a precursor on some data: 612 of
  1,381 against MS-Deconv's 855 on st_1.mzML (MS/MS-matched features 894
  against 1,158), 4,289 against 4,303 on 480_F3_02.
- **`--deconv-ms1 msdeconv` writes the earlier default's output** (with
  `-f 0` since 2026-10-12, when the MS1 EnvCNN cutoff 0.2 began to apply to
  it): on st_1
  all 1,418 MS1 scans' masses and charges, both feature files and the
  parameter block equal main 89df2a8c's (the last main with the DP as the
  default), the MS/MS spectra up to the EnvCNN jitter of two runs of one
  binary. The greedy search itself is not exactly repeatable (about 5% of
  st_1's MS1 scans differ between two runs).
- **Not done**: the precursor-window step of `deconvMsOne` (switch
  `prec_win` runs it first; it did not change the identifications) and
  FAIMS/`-T` specifics (the MS1 path is the same). The MS/MS spectra have
  their own greedy search (next section).

## Greedy MS/MS deconvolution (`topfd/deconv/greedy_ms2_deconv`, TopFD's default)

`greedy_ms2_deconv::deconv` is TopFD's MS/MS deconvolution
(`TopfdPara::isMs2GreedyDeconv`, true by default; `deconv_ms2_process::
deconvMsTwo` and `topfd -T` call it; the parameter block then prints
"MS/MS deconvolution: Greedy" instead of the DP's ranking window and
harmonic check lines). `DeconvSingleSp`, TopFD's default up to main
89df2a8c, runs with `--deconv-ms2 msdeconv`, TopDIA's default (see
"Deconvolution methods" above). `-g`, `-v` and `-n` apply to both methods
(`greedy_ms2_deconv::deconv` ends with the same filters). It started as a
copy of the greedy MS1 search (previous section) and keeps its structure;
the output (msalign lines, EnvCNN scores, the `-F`-relevant score, SQLite
rows) is written as `DeconvSingleSp`'s. Facts to keep straight:

- **Bounds**: fragment masses below the precursors' largest mass and
  charges up to their largest charge (`deconvMsTwo`'s `max_mass` /
  `max_charge`), the MS/MS signal-to-noise ratio for the seed threshold.
- **What differs from the MS1 search** (each a `GreedyMs2Para` switch,
  toppic_test experiment 106):
    - every charge competes with the best of its isotope offsets
      (`charge_offsets`): below about 5,000 Da the averagine apex and the
      monoisotopic peak are of similar height, and the apex alignment of the
      right charge was often one isotope off and lost to a wrong charge;
    - the offset is chosen on the seed charge's envelope alone
      (`offset_by_seed`), and a neighbouring charge joins the series only
      when its own best offset agrees (`nb_consistent`); a neighbouring
      charge was often the envelope of another mass one isotope away;
    - the search matches peaks within 5 ppm (`ppm_tol`), not TopFD's 0.02
      m/z: on 480_F3_02 the peaks of confident envelopes lie within 0.6 ppm
      (median), and the wide window let high-charge envelopes collect the
      peaks of other species; a peak only within 0.02 m/z (two envelopes'
      peaks merged into one centroid) is neither matched nor counted as
      missing (`near_no_pen`). The reported envelopes are matched with
      TopFD's tolerance, as EnvCNN expects;
    - the same-charge duplicate test only blocks masses an integer number
      of isotope spacings from an earlier one (`dup_iso`): distinct
      masses of one charge can lie closer than 2 Da;
    - no harmonic test and a single-charge cosine of 0.7 (most fragments
      appear at one charge);
    - the charge pre-check also accepts a neighbour peak at the averagine
      distribution's own isotope spacing (`precheck` 2): the M+1 peak of
      averagine is 1.00286 Da from the monoisotopic peak at every mass,
      while the isotope mass constant (`mass_constant::getIsotopeMass`,
      1.00228 Da, the apex spacing of human averagine envelopes above about
      5 kDa; 1.00235 before toppic_test experiment 111) fits the apex
      (toppic_test experiments 107 and 108); at m/z 200 and charge 1 the
      difference is 2.9 ppm, more than half of the 5 ppm window.
- **Shared peaks**: the MS1 search's test (previous section, the
  same-charge rule) with the hidden `--ms2-max-shared-peak-num` (default
  -1, no limit, as for MS1; also MS3 spectra and `-T`); not for the alternative
  masses of `alt_ratio` and the refits of `refit_rounds` (both off). On
  480_F3_02 ecf67378's MS/MS envelopes form 126,804 pairs that share peaks
  (99,798 one peak); 20,448 are of one charge: 1-3 isotopes apart (13,454,
  within 0.02 Da or looser), 4 or more apart (6,362, mostly 4-5 Da) and
  632 near-identical masses. Of the pairs sharing one peak 10% are
  harmonics (charges 2:1 or 3:1). Any charge (2df670d5, MS1 unlimited):
  k = 0 / 1 / 2 removes 95,834 / 29,389 / 8,905 of the 701,494 masses
  (median EnvCNN 0.18 / 0.14 / 0.15) and adds 25,671 / 15,505 / 6,410
  (median 0.05 / 0.03 / 0.02: seeds the removed envelopes no longer
  subtract); toppic finds 3,351 / 3,377 / 3,375 PrSMs and 1,981 / 2,005 /
  1,999 proteoforms. Both levels at k = 1: any charge removes 31,025 MS/MS
  masses (4.4%), same charge 11,940 (1.7%; median EnvCNN 0.36); the PrSMs
  unchanged against no limit match 0.7% (any charge) and 0.2% (same
  charge) fewer fragments, and the same-charge PrSMs match 0.5% more
  fragments than the any-charge ones. k = 1 is the DP's
  coexistence rule (`CoTable::checkCoexist`, any charge); any charge at
  k = 0 costs identifications, every other setting is within 7 PrSMs of no
  limit.
- **Research switches**: `TOPFD_GREEDY_MS2_PARA` ("name=value,...")
  overrides the defaults in `GreedyMs2Para`; `TOPFD_GREEDY_MS2_DEBUG_SCAN=
  <scan>` writes one scan's decisions to `greedy_ms2_scan<scan>.txt`.
  Switches kept off by default are listed there with their experiment 106
  results (`refit_rounds`, `alt_ratio`, `mono_by`, `score_mode` 2, ...).
- **Results** (toppic_test experiment 106, report 0106): on 480_F3_02
  with main's TopPIC in the benchmark setting, 3,296-3,298 PrSMs, 1,961-
  1,964 proteoforms and 427 proteins against main's 3,141 / 1,870 / 420 on
  the same MS1 output (+199 to +228 targets at E <= 1e-5 to 0.01, no more
  decoys); with the greedy MS1 search too 3,315 / 1,980 / 433. The gain is
  15% more b/y fragments of the identified proteoforms at main's precision
  (mostly 1,000-5,000 Da: right charges and isotopes, and envelopes the DP
  did not report). The MS/MS stage takes about 100 s against main's 135 s
  at 16 threads.
- **`--deconv-ms2 msdeconv` writes the earlier default's output**: on
  st_1, with `--deconv-ms1 msdeconv` too, the MS/MS spectra equal main
  89df2a8c's default output up to the EnvCNN jitter, and `topfd -T` on the
  myoglobin peak list writes main's msalign byte for byte (681 masses; the
  greedy search 2,476). With the MS1 deconvolution fixed, the greedy MS/MS
  search is repeatable up to that jitter. `--deconv-ms2 flashdeconv` writes
  main 89df2a8c's `-F` MS/MS spectra unchanged (st_1, all 1,381).

## Precursor isotope correction from fragment pairs (`topfd/deconv/prec_pair_refine`, `--skip-ms2-prec-pair-refine`)

On by default in topfd; the hidden `--skip-ms2-prec-pair-refine` turns it
off (branch `high_res_prec_pair_refine`, toppic_test experiment 116, report
0116; tested as the hidden `--ms2-prec-pair-refine`, which the merge into
main replaced; `--skip-ms2-prec-pair-refine` was visible in `-h` until
2026-10-08, when it was moved to the hidden options).
b + y = the proteoform mass, so fragment pair sums show the mass of the
molecule that was fragmented. Facts to keep straight:

- **Two functions.** `prec_pair_refine::refine` is the low-resolution
  method (`-R`, branch `low_res_greedy_ms1` f902aa09, copied unchanged; main
  has no `-R`, so nothing calls it here; a later merge of a low-resolution
  branch has to rename its `--ms2-prec-pair-refine` too). `correctIsotopes`
  is the high-resolution one that TopFD runs, in `deconvMsTwo` after the
  precursors are built from the features.
- **Where it runs.** `TopfdPara::ms2_prec_pair_refine_` is true by default.
  It runs in `DeconvMs2Process` only (greedy and msdeconv MS/MS
  deconvolution, not FLASHDeconv or `-T`) and only for MS/MS spectra with an
  MS1 spectrum: without one (`-o`) the MS1 check cannot run, and with the
  pairs alone the rule moved 47 right precursors of the known-mass spectra
  away and fixed 18 (report 0116 section 3).
  TopDIA's argument parser turns it off (DIA MS/MS scans hold the fragments
  of every co-isolated precursor; not tested). The parameter line
  `Precursor isotope correction from fragment pairs: Yes` is printed only
  when the correction runs, so the other output files keep the parameter
  block of earlier versions.
- **The rule (`correctIsotopes`).** For every precursor of a spectrum (TopPIC
  searches only the first), it counts the pair sums of the fragment masses
  with EnvCNN >= 0.4 within 10 ppm of P - 1, P and P + 1 isotope (1.00228 Da).
  The precursor moves to the best offset when it has >= 3 pairs, more than at P
  and >= 2 x the count at P. It does not move when the averagine fit of the
  precursor in the group's MS1 spectrum prefers P by 0.05 or more (cosines on
  one grid of isotope positions, +-2 isotopes around the envelope, summed
  over the charges z-1..z+1).
    - The moved precursor keeps its charge, feature id and intensity; its m/z
      moves by j x 1.00228 / z.
    - `DeconvMs2Process::process` passes the group's MS1 spectrum to the
      MS/MS task only when the correction is on.
- **Research switches** in `TOPFD_PREC_REFINE_PARA`: `iso_range` (1),
  `iso_tol_ppm` (10), `iso_min_envcnn` (0.4), `iso_min_pairs` (3),
  `iso_ratio` (2), `iso_ms1_check` (1), `iso_ms1_margin` (0.05),
  `iso_ms1_tol_ppm` (10), `iso_ms1_min_rel` (0.05). The low-resolution
  `window`, `tol_ppm`, `min_pairs` and `promote` are unchanged.
- **Results on 480_F3_02** (`-s 2`, oxidation only):
    - It moves 96 first precursors. At 1% FDR TopPIC finds 3,341 PrSMs /
      2,003 proteoforms / 438 proteins, against main's 3,337 / 1,997 / 438.
    - The moved spectra get better PrSMs: of 86 identified by either run,
      56 match more fragments and 10 fewer.
    - 12 of them are Asn-Gly deamidations (HMGN1, parathymosin, hnRNP A1-like
      2) that main reports without the deamidation. A deamidated proteoform is
      0.018 Da short of one isotope heavier, and MS1 cannot separate its
      envelope from the unmodified one.
- **Why the rule is conservative.** TopFD's fragment masses carry isotope
  errors, mostly one isotope light: they are right in 98% of the matches
  below 1 kDa but 78% at 3-4 kDa and 55% at 6-10 kDa. A pair sum has the
  isotope error of its larger fragment. With pairs alone (no MS1 check,
  moves of +-2) 240 precursors move and many PrSMs get worse.
- **TopFD's own precursor errors** are rare: 97 of 1,866 known-mass spectra,
  half of them histone H2B at 13.8 kDa (one isotope light). The pairs cannot
  reach those: complementary pairs are missing or flat there.

## PTCR spectra (`topfd/deconv/ptcr_ms2_deconv`)

Branch `ptcr_ms2_deconv` (toppic_test experiment 136, report 0136; merged
into main on 2026-10-07 as 5f19b7cb: with the branch's float EnvCNN / ECScore
models the merged topfd writes the branch's PTCR, fragment PTCR and LC-MS
outputs, and main's quantized models change only EnvCNN scores there, no
mass line; 480_F3_02 within the variation of two runs of main faa8e9c3;
report 0136 section 9). A PTCR spectrum is an MS/MS spectrum of proton transfer charge reduction of the
isolated ions without dissociation: no fragments, only the isolated intact
species at lower charges. Facts to keep straight:

- **Reading.** `PwMsReader::parseActivation` decides on the activation
  terms of all precursors of a spectrum (an MSn scan has one per stage;
  msconvert writes the first stage first). `getActivationName` names the
  dissociation (proton transfer is none); `hasPtcrTerm` finds MS:1003249
  (proton transfer charge reduction) or MS:1003248 (proton transfer
  reaction). A proton transfer term without a dissociation: a PTCR
  spectrum (`MsHeader::isPtcr`, an activation object of its own,
  `getPtcrActivationPtr`, placeholder ion types B / Y). `PTCR` is
  deliberately not in `activation_base.xml`, the `-a` lists or the GUIs: no
  search reads these spectra. A dissociation (of the file or of `-a`) and a
  proton transfer term: a fragment PTCR spectrum (`isFragPtcr`, below), of
  the dissociation's activation, when the proton transfer stage comes after
  the last dissociation stage or in its stage; so `-a HCD` on a PTCR scan
  now makes a fragment PTCR spectrum, not main's fragment spectrum. **The
  order of the stages** decides (toppic_test experiment 139):
  `getScanStages` lists them in acquisition order (the filter string's
  stages, each with the file's precursor of the same target, 0.01 m/z;
  without a filter string the precursors in file order, msconvert's stage
  order), and when a dissociation stage comes after the last proton
  transfer stage (`getReducedPrecStage`) the spectrum is the dissociation
  of a charge-reduced precursor (`isReducedPrec`, below), not a fragment
  PTCR spectrum. Main 6eb40e97 read the MS3 scans "X@ptr5 Y@hcd80" of
  PRIDE PXD074730 as fragment PTCR spectra with the window of the proton
  transfer stage and, without a charge, `-o`'s placeholder precursor. The hidden `topfd
  --ptcr` sets the activation to `PTCR` for every MS/MS spectrum
  (`TopfdPara::isPtcrForced`, which the activation check accepts), for
  files whose converter wrote no proton transfer term. `parsePrecursor`
  keeps the selected ion's m/z and charge of the first precursor
  (`MsHeader::getSelectedMz` / `getSelectedCharge`, -1 / 0 without) and the
  isolation window of the proton transfer stage (`getPtcrWinBegin` /
  `getPtcrWinEnd`: the last precursor with a proton transfer term, else the
  precursor window).
- **MSn scans** (n > 2) with a proton transfer stage are MS/MS spectra
  (`MzmlMsGroupReader`'s `isMsTwo`, in `-o` and with MS1 spectra, and the
  count of `readProfile`): grouped with the last MS1 scan before them
  (`parsePrecursor` ignores their precursor reference, the MS(n-1) scan).
  Other MSn scans (a dissociation of fragments, whose precursor is not the
  window's proteoform) are skipped, as main skipped every MSn scan.
  **Missing stages** (`getMissingStages`): ThermoRawFileParser 2.0 writes
  one precursor for the MS3 scan "ms3 1402.5000@cid60.00
  1100.0000@ptr10.00", the proton transfer stage (window 1100 +- 900, the
  PTR term, the scan's charge state), and loses the CID. When a scan has
  fewer precursors than stages and its filter string (MS:1000512) lists
  them, the stages' terms (cid, hcd, etd, ptr) join the activation terms,
  and when the file's precursor is not the first stage, the precursor
  window and the selected ion m/z are the first stage's target with the
  default window (`-w`); the charge stays the file's. The same holds for
  **MS2 scans without a precursorList** (toppic_test experiment 139):
  ThermoRawFileParser writes none for a scan whose precursor scan is not in
  the file (a scan cut out of a longer run, as all MS2 and MS3 files of
  PXD074730), and main stopped on them ("Unable to read the activation
  method"); the filter string's stage now gives the target with the
  default window, the selected m/z (no charge) and the activation.
- **Where it runs.** `deconv_ms2_process::deconvMsTwo` sends a PTCR
  spectrum to `ptcr_ms2_deconv::deconv` and buffers the result
  (`PtcrBuffer`, under a lock); `DeconvMs2Process::process` writes
  `<base>_ptcr.msalign` (msalign format, `ACTIVATION=PTCR`, one line per
  species: mass, intensity, highest charge, best EnvCNN score; `-A` average
  masses) and `<base>_ptcr.tsv` (one line per species with its charges, the
  intensity and EnvCNN score per charge, the isolation charges, the
  precursor flag) after the pool, sorted by spectrum id, only when there are
  PTCR spectra. They are not in `_ms2.msalign` or the SQLite database. In a
  run with MS1 spectra, a species whose mass matches an MS1 feature of the
  spectrum's window (`feature_ppm` 10 ppm, also one isotope apart) gets its
  id and mass (the TSV columns `Feature ID`, `Feature mass`), and the
  precursor species' feature gives the precursor its id and intensity.
  Not in `FlashDeconvMs2Process` (`--deconv-ms2 flashdeconv`) or `-T`.
- **The search** is the greedy MS1 search through
  `greedy_ms1_deconv::deconvSeries`, a second entry point of
  `greedy_ms1_deconv` (no precursor-window step) that takes a seed S/N
  (the MS/MS ratio `-s`) and a `ChargeBound` (`fitEnv` rejects an excluded
  (mass, charge): seeds, series and offsets) and returns the envelopes
  grouped by mass (EnvCNN scored, below 0.2 dropped, intensities shared
  among the kept ones). The PTCR bound: charge z of mass M only when some
  z0 >= z puts the envelope of M (its isotopes >= 10% of the apex in
  TopFD's averagine table) in the isolation window, since proton transfer
  only lowers charges. Without an isolation window there is no bound.
- **A species** is one mass's envelopes merged: the summed intensity, the
  charges, the mass refined from the matched isotope peaks of all charges
  (the intensity-weighted median of (observed - theoretical m/z) x z; the
  envelope's own mass rests on the seed peak). **The precursor** is a
  species with an envelope at the selected ion's charge whose monoisotopic
  m/z at that charge is the selected ion's m/z (`sel_ppm`, 10 ppm), else
  the most intense species with an envelope at that charge whose m/z
  range, from its monoisotopic peak to its last isotope of >= 10% of the
  apex, covers the selected ion's m/z, else the most intense one that could
  have been isolated at that charge, else the most intense one, else the
  selected ion itself. Thermo converters write the instrument's
  monoisotopic m/z (the "Monoisotopic M/Z" trailer) as the selected ion
  when it has one, else the isolation target. The cover test began at the
  first isotope of >= 10% (toppic_test experiment 136), at 17 kDa several
  isotopes above the monoisotopic peak, so a monoisotopic selected m/z lay
  in the envelope of the 16 Da lighter species only: on PXD074730's
  Lysozyme_MS2_PTCR (selected 941.9354, 18+: 16,936.70 Da) main 6eb40e97
  took the 16,920.73 Da species (6.9 M) instead of the 16,936.69 Da one
  (80.7 M; experiment 139).
- **Fragment PTCR spectra** (`isFragPtcr`: a dissociation, then proton
  transfer charge reduction of the fragments) are deconvoluted as other
  MS/MS spectra of their activation (`deconvMsTwo`, `--deconv-ms2`) and
  written to `_ms2.msalign` and the SQLite database, so TopPIC searches
  them with the dissociation's ion types. The selected ion bounds the
  fragments (`ptcr_ms2_deconv::getSelectedIonPrec`: its charge, and its
  mass at the upper end of the precursor window; not capped by `-c` /
  `-m`), with MS1 features the larger of these and the features' bounds (a
  feature of the window need not be the isolated ion: the AAV2 MS1 scans
  are not isotope-resolved and their features are noise, which as the only
  bound left 470 masses per spectrum instead of 4,100-5,500). Without MS1
  features of the window the selected ion (its m/z and charge as the file
  gives them, not monoisotopic) is the precursor. The console counts them
  ("Fragment PTCR spectra: N"). Switch `frag`: 1 this (default), 0 as other
  MS/MS spectra (`-o`: the dummy precursor of `-c` / `-m`), 2
  `deconvFragments`, the charge-series search of the intact spectra with
  the PTCR bound, z0 <= the selected charge and masses below the
  precursor's (also its surviving precursor ladder as the precursor), one
  mass per fragment species, `<base>_ptcr_fragment.tsv`. The MS/MS search
  matched more theoretical fragments of AAV2's capsid proteins than the
  series search (below): with the wide PTCR windows of these data
  (1100 +- 900, 650 +- 255 m/z) the band bound excludes almost nothing
  above 1 kDa, and the MS/MS search's fragment-specific choices (5 ppm,
  per-charge isotope offsets, no harmonic test, the single-charge cosine
  0.7) find more low-charge fragments.
- **The dissociation of a charge-reduced precursor** (`isReducedPrec`;
  toppic_test experiment 139): an MSn scan whose dissociation stage comes
  after its proton transfer stage, e.g. the MS3 scan "X@ptr5 Y@hcd80" of
  the cycle MS1, PTCR MS2 of X, HCD MS3 of one charge-reduced ion at Y
  (PXD074730, iodoTMT quantification of intact proteins: the MS3 scans
  hold the reporter ions). It is read as an MS/MS spectrum of the
  dissociation's activation (`isMsTwo`, `isReadMsTwo`) and written to
  `_ms2.msalign` (`LEVEL=3`) and the SQLite database; its precursor window
  is the first stage's, the MS1 window of the isolated ions, matched to
  the MS1 features as for every MS/MS spectrum. The header keeps the
  dissociation stage's target and window (`getReducedPrecTargetMz`,
  `getReducedPrecWinBegin` / `End`: the file's precursor of that target,
  else the filter string's target with the default window `-w`) and the
  scan its last precursor references (`getParentScan`). Its precursor
  (`deconv_ms2_process::deconvMsTwo`, steps 0b and 1c): (a) from the
  parent PTCR spectrum (`findPtcrParent`: the referenced scan, else the
  last PTCR spectrum read before it with the same isolation target, 0.01
  m/z, among the last 64), deconvoluted again in the MS3's task
  (`ptcr_ms2_deconv::getReducedPrec`; `getPtcrSpec` reads the parent's
  header and does not change it, so both tasks can run at once): the most
  intense species with an envelope at one of its charges, else at any
  charge up to its isolation charge, whose m/z range (monoisotopic peak to
  the last isotope >= 10%) contains the dissociation stage's target, at
  that charge, with its intensity and MS1 feature; (b) without a parent,
  the MS1 features of the window with an envelope at a charge up to theirs
  containing that target (`getReducedPrecFeatures`), all of them when none
  has one; (c) without MS1 features the selected ion with its charge
  (`getSelectedIonPrec`); (d) else `-o`'s placeholder. The precursor's mass
  and charge bound the fragments. The console reports "Dissociation
  spectra of charge-reduced precursors: N, M with the precursor from a
  PTCR spectrum". The MS3 scan does not tell the charge of the ion at Y;
  the parent's charge ladders do. The precursor mass is the parent
  species' own, so it shares the species' isotope error (below).
- **Research switches** in `TOPFD_PTCR_PARA` ("name=value,..."): `band`
  (1), `win_margin` (0 m/z, added to both window ends), `iso_min_rel`
  (0.1), `min_envcnn` (0.2), `sn_ratio` (0 = `-s`), `refine` (1), `frag`
  (1, above), `frag_min_envcnn` (0, frag 2), `frag_prec_margin` (20 Da,
  frag 2: species this close to the precursor mass or above are not
  fragments), `feature_ppm` (10), `sel_ppm` (10, the monoisotopic
  selected ion). The greedy search's own switches
  (`TOPFD_GREEDY_PARA`) apply too.
  `TOPFD_PTCR_DEBUG_SCAN=<scan>` writes one scan's decisions to
  `ptcr_scan<scan>.txt`.
- **`-o` lost the first MS/MS spectrum of every file without MS1
  spectra**: `MzmlMsGroupReader`'s constructor reads one spectrum to check
  that the file is not empty and rewound only when MS1 spectra are present
  (fixed in this branch). TopDIA's MS/MS step runs in that mode on files
  that start with an MS1 spectrum, which the skip consumed harmlessly.
- **Results** (toppic_test experiment 136; data/ptcr_ms2, two single PTCR
  scans of an Orbitrap Eclipse, isolated at 16+ and 17+ with +-5 m/z, of a
  16.94 kDa proteoform and its co-isolated variants; ThermoRawFileParser
  drops the precursor of these scans, cut out of longer runs, and the
  experiment's add_precursor.py writes it back): the dominant proteoform at
  16,936.657 / 16,936.663 Da in the two scans (0.3 ppm; 0.9 ppm without the
  refinement), 11+-16+ and 11+-17+, as their precursor; 77 / 95 species,
  no 1-isotope twins; species with 75% / 72% of the intensity found in both
  scans at the same mass (median 0.47 ppm apart), 17% / 16% at an isotope
  shift. The co-isolated variants lie a near-whole number of Daltons apart,
  so their isotope peaks merge into one centroid (5 ppm at every charge),
  and a minor species' monoisotopic mass is not determined to the isotope
  (a joint least-squares fit of all species does not decide it either).
  Without the bound 81 / 102 more species (9% / 10% of the intensity),
  2 of them at the same mass in the other scan. One species whose envelope
  ends 0.8 and 0.5 m/z below the two windows is reproducible (16,835.55 Da,
  1% of the intensity): the window edges are not sharp, but margins of
  0.5-3 m/z added mostly irreproducible species, so the margin stays 0. main 38ba67dc writes 420 /
  487 envelope lines (one per charge state) with a dummy precursor (`-o`:
  charge 30, mass 50,000 Da), 21 / 33 masses above 10 kDa with a
  1-isotope twin, 72% / 77% of the lines (15% / 19% of the intensity) of
  masses outside every band, charges up to 30. On 480_F3_02 (no PTCR
  spectra) the branch's TopFD output equals main's within the variation of
  two main runs (MS/MS spectra identical; 702 against 727 MS1 spectra with
  a differing mass; same time and memory).
- **Results of fragment PTCR and LC-MS runs** (toppic_test experiment 136,
  amendment 1; MassIVE MSV000095045, AAV2 capsid on an Orbitrap Ascend: 36
  MS3 scans "ms3 1288|1300|1402.5@cid60 1100@ptr10" and 36 "825|830@etd3
  650@ptr10|8", matched against the b/y or c/z-dot ions of VP1, VP2 and VP3
  at 10 ppm with reversed-sequence decoys): the MS/MS search (frag 1) 1,125
  CID and 334 ETD target matches above the decoys' (sums over the
  spectra), the series search (frag 2) 937 and 206, frag 0 with `-o`'s
  bounds 1,113 / 320 and with `-c 64 -m 90000` 1,122 / 335. Main reads no
  MS3 scan (0 MS/MS scans). ThermoRawFileParser's own files (first stage
  from the filter string) give the masses of the files with a precursor per
  stage in 69 of 72 spectra (3 ETD spectra differ near the bound: the
  default window of 3 m/z against +-2.5). A synthetic LC-MS run (the
  AAV2_MS1 run's 2,882 MS1 scans with the 72 MS3 scans and 3 MS1 scans
  copied as intact PTCR scans): every MS3 scan in the group of the MS1 scan
  before it, the masses of the `-o` runs in all 72; 62 precursors from
  features (noise: 43 of charge 1), 10 the selected ion; an intact copy's
  species linked to its window's feature. The intact PTCR outputs of
  data/ptcr_ms2 are unchanged (but the two TSV columns); 480_F3_02 within
  the variation of two main runs.
- **Results on PRIDE PXD074730** (toppic_test experiment 139, report 0139;
  Kaulich and Tholey 2026, iodoTMT-labeled intact lysozyme, BSA and
  beta-lactoglobulin by direct infusion on an Orbitrap Eclipse: MS1, PTCR
  MS2 of X, HCD MS3 of one charge-reduced ion at Y with the reporter ions;
  19 single scans cut out of runs, ThermoRawFileParser 2.0.0-dev, 18 without
  a precursorList). Main 6eb40e97 stopped on the 11 MS2 files without a
  precursorList, took the 16,920.73 Da species as the precursor of
  `Lysozyme_MS2_PTCR` (selected 941.9354 at 18+ = 16,936.69 Da, the
  strongest species), and read the 7 MS3 scans as fragment PTCR spectra
  with the placeholder precursor. The branch reads all 19; with each PTCR
  scan and its MS3 scans in one file, the 6 MS3 precursors come from the
  PTCR species at Y (lysozyme 16,936.66-16,936.69 Da at 15+ to 18+, within
  0.6 ppm of the instrument's monoisotopic m/z of Y in the MS3 scan's
  trailer; beta-lactoglobulin 19,915.54 Da at 18+, one isotope above the
  instrument's mass and its MS1 feature), also from the stage precursors
  alone (no filter string) and with a spectrumRef; with MS1 scans (one
  run's 396 MS1 scans with the PTCR scan and two MS3 scans inserted) each
  MS3 gets its own proteoform and feature, where main gave both scans both
  features. Fragment PTCR on AAV2 and the 4 MS1-only files unchanged
  (masses; EnvCNN scores by up to 0.02); on 480_F3_02 the reading fixes
  alone (7e5b5a64) are within the variation of two main runs, and with the
  EnvCNN input (above) toppic's identifications equal main's.

## MS3 spectral deconvolution (`topfd/deconv/deconv_ms3_process`, `--disable-ms3`)

Branch `ms3_deconv` (toppic_test experiment 146, report 0146; from main
205dc03c, merged into main as ceab7362 on 2026-10-08). An MS3 scan is the
dissociation of an MS2 fragment ion (a
tryptic digest on an Orbitrap Fusion Lumos: CID MS2 of the peptide, CID
MS3 of up to three of its fragments). TopFD read only the MS3 scans with a
proton transfer stage (PTCR, above); the others were skipped. Since this
branch they are deconvoluted by a stage after the MS/MS deconvolution,
`DeconvMs3Process`, on by default when the file has such scans
(`MzmlProfile::getMs3Cnt`, counted by `PwMsReader::readProfile`; the
parameter block prints "Number of MS3 scans" and "MS3 spectral
deconvolution: on / off"; hidden `topfd --disable-ms3`; not run for
FAIMS data). Facts to keep straight:

- **Reading.** The stage reads the mzML file again with `PwMsReader`
  (which reads every scan; the MS/MS stage's `MzmlMsGroupReader` is what
  skips MS3 scans) and keeps the level-3 scans that are not PTCR,
  fragment PTCR or charge-reduced-precursor spectra. The parent MS2 scan
  of an MSn scan (`MsHeader::getParentScan`) is the `spectrumRef` of the
  precursor element whose isolation target is the **last** stage of the
  filter string (`getFilterStages`): ThermoRawFileParser 2.0 writes the
  MS3 stage first and the MS1 stage second, msconvert the MS1 stage
  first; without a filter string the last precursor (msconvert's order).
  The MS1 scan is the cycle's (`MS_ONE_SCAN`, `prev_ms1_scan_id_`). The
  isolation target and window are the first precursor element's
  (`parsePrecursor`), which is the MS3 stage for ThermoRawFileParser
  files - with msconvert files it would be the MS1 stage (not tested).
  The activation is the stages' (CID).
- **The precursor** of an MS3 spectrum is an MS2 fragment: the species
  of the parent MS2 scan's deconvolution (`<base>_ms2.msalign`, read
  back by scan) whose isotopic envelope (the monoisotopic m/z and the
  `1 + mass / 800` isotope peaks above it) has a peak in the MS3
  isolation window, the most intense when several (`findPrecSpecies`);
  its mass, charge, m/z and intensity are the precursor
  (`PRECURSOR_FEATURE_ID` -1: a fragment has no MS1 feature). Without
  one - in the test data the MS3 scans of the noise MS2 scans at the
  gradient's ends, whose MS2 precursor has the dummy charge 1 - the
  selected ion m/z of the file at every charge from 1 to the MS2
  precursor's charge, written as multiple precursors (the file's charge
  state of an MS3 precursor is not usable: 10 for a third of the scans).
  On the test data the chosen species is a b or y ion of the identified
  parent peptide in two thirds of the cases and one isotope off (the
  MS2 deconvolution's monoisotopic error: -1 Da for 2-6 kDa fragments at
  charge 3-6, +1 Da below 3 kDa) in a quarter, so the species' mass is
  **corrected by one isotope from the parent MS2 scan's centroid peaks**
  (`correctIsotope`; the stage keeps the MS2 scans' peaks as it reads
  them): when the peak one isotope above the reported monoisotopic m/z
  is more than `TOPFD_MS3_ISO_K` (4) times the theoretical ratio
  (`EnvBase::getEnvByMonoMass`) times the reported peak, the mass goes up
  one isotope; else when a peak one isotope below holds at least
  `TOPFD_MS3_ISO_C` (0.1) of the reported peak, it goes down one;
  `TOPFD_MS3_ISO_CORRECT=0` turns it off and the `_ms3.tsv` column
  "Isotope correction" records the shift. A cosine of the whole
  theoretical envelope against the peaks was tried first and did not
  separate the cases (experiment 146, amendment 1). TopPIC's zero-shift
  search tolerates the +-1 Da precursor error anyway.
- **The deconvolution** is the MS/MS one (`greedy_ms2_deconv::deconv`,
  or `DeconvSingleSp` at level 2 with `--deconv-ms2 msdeconv`) bounded by
  the precursor's mass (+1 Da) and charge (the MS2 precursor's in the
  fallback case), single-threaded (8,000 MS3 scans of 75 peaks take
  about 80 s).
- **Output.** `<base>_ms3.msalign` in the msalign format with `LEVEL=3`,
  `MS_ONE_ID` / `MS_ONE_SCAN`, and the new lines `MS_TWO_ID` /
  `MS_TWO_SCAN` (the parent MS2 spectrum's id in `_ms2.msalign` and its
  scan) and `MS_TWO_PRECURSOR_MZ` / `_CHARGE` / `_MASS` (the parent's
  first precursor; `MsHeader::setMsTwoId`, `setMsTwoPrec`; written by
  `MsAlignWriter` for level 3, ignored by `MsAlignReader`), the
  precursor lines and the masses; `<base>_ms3.tsv`, one row per MS3 scan
  (scan, MS2 scan and id, MS1 scan and id, retention time, the MS2
  precursor, the isolation target and window, the precursor source
  `ms2_species` / `selected_ion`, the precursor(s), peaks, masses). The
  MS1 / MS2 msalign and feature files and the SQLite database are
  unchanged (checked file for file against main). TopPIC searches an MS3
  msalign file with `-x` (no feature file) and `--disable-post-match`
  (the SQLite database has no MS3 peaks); a y-ion precursor is a
  both-truncated proteoform with no shift, a b-ion one carries -18.011
  Da as its unexpected shift.

## Activation of MS/MS spectra (`ms/mzml/pw_ms_reader`, `res/base_data/activation_base.xml`)

With `-a FILE` (the default) `PwMsReader::parseActivation` takes each
MS/MS spectrum's activation from the activation cvParams of its first
precursor, through `getActivationName` (branch `ethcd_activation`, issue
#530, toppic_test experiment 132). Facts to keep straight:

- **The set of terms decides, not their order.** MS:1002631 (EThcD), or ETD
  (MS:1000598) with beam-type CID (MS:1000422, pwiz `MS_HCD`) or supplemental
  beam-type CID (MS:1002678) -> `ETHCD`; MS:1003182 (ETciD), or ETD with CID
  (MS:1000133) or supplemental CID (MS:1002679) -> `ETCID` (EThcD is tested
  first, so ETD with both kinds of CID is `ETHCD`); otherwise `ETD`, `HCD`,
  `CID`, `MPD` (MS:1000435) in that order; no known term -> the "Unable to
  read the activation method" error, as before. The reader used to keep the
  first of CID/HCD/ETD/MPD in file order, and the converters order the terms
  differently: msconvert writes CID, ETD, ECD, PQD, HCD, then the supplemental
  term (EThcD: ETD, beam-type CID, supplemental beam-type CID, read as ETD;
  ETciD: CID, ETD, supplemental CID, read as CID), ThermoRawFileParser the
  first reaction, then the supplemental one (EThcD: ETD, supplemental
  beam-type CID). Neither writes MS:1002631 or MS:1003182, so do not wait for
  those terms alone. The vendored pwiz (CV 4.1.182) knows all of them;
  updating it changes nothing here.
- **`ETHCD` and `ETCID` have the ion types of ETD** (C / Z_DOT in
  `activation_base.xml`): TopFD writes `ACTIVATION=ETHCD` or `ETCID` in the
  msalign file and C / Z_DOT in the SQLite table `ms2_spectrum`, and toppic
  and topmg read the name back (`ActivationBase::getActivationPtrByName`)
  and search EThcD and ETciD spectra as ETD ones. Searching b / y (or c / y
  for ETciD's CID) ions with c / z-dot would need several ion types per
  `Activation` (one N- and one C-terminal ion type now; used by
  `prsm/theo_peak_util`, the PRM spectra's `getN_BYShift`/`getC_BYShift`,
  the zero-shift fast search and the MCMC E-values).
- **An activation name is listed in several places**:
  `activation_base.xml`, the `-a` help text and validation of
  `topfd_argument.cpp`, `topdia_argument.cpp`, `toppic_argument.cpp` and
  `topmg_argument.cpp`, the `activationComboBox` items of the five GUI forms
  (`topfddialog.ui`, `topdiadialog.ui`, `toppicwindow.ui`, `topmgwindow.ui`,
  `topmergedialog.ui`; the item text is the option value), the `-a` rows of
  the four manuals and `res/README.md`. Change them together.
- **Results** (toppic_test experiment 132, amendment 1): on
  PXD019792_191018_Mmazei_R2_3_ETHCD_5sec (Orbitrap Fusion Lumos, msconvert
  3.0.24284) all 2,190 MS/MS spectra are ETHCD (main: ETD); the MS/MS
  masses, charges and headers and the feature files equal main's (the MS1
  msalign differs in 7 of 1,590 spectra, main against a second main run in
  9: EnvCNN and the greedy search are not exactly repeatable). toppic on
  main's TopFD output gives main's files; with `ACTIVATION=ETHCD`, or `-a
  ETHCD`, the same but the Fragmentation column (346 PrSMs, 206
  proteoforms, 66 proteins at 1% FDR against the M. mazei proteome). Test
  spectra in 12 activation cases (msconvert's and ThermoRawFileParser's
  EThcD and ETciD, the two terms alone, EThcD with HCD first, ETD, HCD,
  CID, CID with HCD, MPD) all get the expected activation and ion types;
  the other 15 mzML files of toppic_test's `data/` (34,143 MS/MS spectra,
  HCD or CID) keep theirs. Amendment 2 (`ETCID`): the three ETciD cases
  (msconvert's, ThermoRawFileParser's, the term alone) give ETCID and ETD
  with beam-type CID and CID gives ETHCD; toppic with the label ETCID (from
  the msalign or `-a ETCID`) writes the ETD run's files but the
  Fragmentation column. None of the local files has ETciD spectra.

## Isotope distributions (`common/emass`)

`common/emass` is a C++ port of toppic_test's `emass.py` (experiments 107,
108 and 111): emass's convolutions of isotope patterns (abundance and
abundance-weighted mass offset per number of extra neutrons, trailing peaks
<= 1e-100 pruned, binary powers for whole atoms) and J.C.P. Miller's
recurrence for fractional atom numbers (`BelowOne::kMixture` takes fewer
than one atom as a mixture with and without it). `isotope_data` holds the
isotopes of C, H, N, O and S (AME2020 masses, IUPAC representative
abundances) and the averagine compositions: `kHumanAveragine` (the UniProt
human reference proteome, experiment 108; the default),
`kHumanNoSulphurAveragine` (the same proteome without cysteine and
methionine, so C, H, N and O only; experiment 112) and `kSenkoAveragine`. `averagine_table::writeAveragineTable` writes TopFD's
`res/base_data/theo_patt.txt` and reproduces the file of experiment 111
byte for byte (toppic_test experiment 112, in 2.5 s against `emass.py`'s
14 s). Facts to keep straight:

- **Bit-identical to `emass.py`**, which experiment 112 checked for every
  rounding and fractional atoms on 14 masses from 100 Da to 100 kDa: the
  monoisotopic masses and the averagine unit mass are summed with
  Neumaier's compensation as Python's `sum()` does (Python >= 3.12), the
  convolutions sum each output in the order of the longer series as
  `numpy.convolve` does, and `averagineFormula` rounds halves to even as
  Python's `round()`. The elements are always convolved in the order C, H,
  N, O, S; `emass.py`'s H-filled formulas list H last, which moves those
  distributions by about 1e-16.
- **It does not replace `EnvBase`.** TopFD still reads the table from
  `theo_patt.txt` (`ms/env/env_base`); `common` cannot depend on `ms`, so
  the table's 11,000 entries and 10 Da interval are repeated in
  `averagine_table.cpp`. Keep the two in step.
- **`theo_patt_gen`** (`src/console/theo_patt_gen.cpp`) is the command-line
  front end: `theo_patt_gen [-o file] [-c human|human_no_sulphur|senko]
  [-s k]` (`--output`, `--composition`, `--sulphur-num`; defaults
  `theo_patt.txt`, `human`, 0). It is a
  developer tool, declared
  `EXCLUDE_FROM_ALL`: `make` does not build it (`make theo_patt_gen`
  does), it is deliberately absent from `TOPPIC_EXECUTABLES` (CMake
  documents installing an `EXCLUDE_FROM_ALL` target as undefined behaviour)
  and it has no manual. Regenerate the resource with
  `bin/theo_patt_gen -o res/base_data/theo_patt.txt`. `--sulphur-num k`
  (0-5; k > 0 only with `human_no_sulphur`) gives every formula exactly k
  S (`averagineFormulaWithSulphur`): the composition is scaled to the entry
  mass minus k 32S masses (at least 0, so the lightest entries are S_k
  alone) and k S are added.

## The MS1 peaks of the database and the 3D levels (`sql/ms_sql_writer`, `--sql-3d`)

Since toppic_test experiment 165 the MS1 peaks and the levels of the 3D
view (`topfd -D`) are one table:

    ms1_peak(LEVEL, spec_id, peak_id, mz, intensity) WITHOUT ROWID,
      PRIMARY KEY(LEVEL, spec_id, peak_id)

- **Before**, `ms1_peak` was a rowid table without `LEVEL`, and `-D` added
  `PEAKS0..n` (`MZ`, `INTENSITY`, `RETENTIONTIME` in ms, `COLOR`), each
  with an index on `(RETENTIONTIME, MZ)`. `PEAKS0` repeated every MS1 peak
  of intensity > 0, and each further table kept the most intense peak of
  each grid block of the one below, so the levels are nested. On
  480_F3_02 the tables were 1,641 of the file's 2,516 MB.
- **`LEVEL`** is the highest level that keeps the peak; level k is the
  rows with `LEVEL >= k`. Without `-D` every peak has `LEVEL` 0. A peak of
  intensity <= 0 is not in the map and has -1 (a guard: the mzML reader
  drops such peaks, `PwMsReader::parsePeaks`, so TopFD writes no -1 today).
  `CONFIG` is unchanged (row k + 1: the ranges and the `COUNT` of level k).
- **Not stored**: a peak's retention time is its spectrum's
  (`static_cast<int>(rt * 1000)` ms, as before), and readers compute the
  colour as `setColors` did, from `CONFIG`'s first `INTMIN` / `INTMAX`
  (0 when `INTMIN <= 0` or `INTMAX <= INTMIN`, else
  `clamp(int(7 (ln I - ln INTMIN) / (ln INTMAX - ln INTMIN)), 0, 6)`).
- **Queries**: a spectrum's peaks are `LEVEL IN (-1, 0, ..., n) AND spec_id
  = ? ORDER BY peak_id`, one key range per `LEVEL` (no index on
  `spec_id`). A 3D window of level k is `LEVEL IN (k, ..., n)` over the
  window's spectra: MS1 spectrum ids follow the retention time, so a
  retention-time range is a `spec_id` range. TopMSV tells the layouts
  apart by the `LEVEL` column (topmsv_private branch `merged_ms1_peak`).
- **Who writes it**: without `-D`, `MzmlMsSqlWriter` as before, with
  `LEVEL` 0 or -1. With `-D`, `deconv_ms1_process` keeps every MS1 peak with
  its spectrum's id and its index (the key `MzmlMsSqlWriter` would write),
  `MzmlMsSqlWriter` writes no `ms1_peak` rows, and `MsSqlWriter::write`
  builds the levels, writes `CONFIG` and inserts `ms1_peak` in key order,
  in one transaction after the deconvolution.
- **The levels are those of before only if their input order is.**
  `downsample` and `addRow` use `std::sort`, which is not stable, so the
  input order decides between peaks of equal retention time or m/z; a
  block keeps the first of two equally intense peaks. `write` therefore
  builds the levels from the peaks of intensity > 0 in the order
  `deconv_ms1_process` collected them, sorts each level by
  `(rt, mz)` (`sortByRtMz`, the order the `PEAKS<k>` tables were written
  in) before the next is down-sampled, and records each peak's level
  through `MsSqlPeak::pos`, its input position (the levels are copies).
  Level 0 is kept to the end and put back in input order in place
  (`restoreInputOrder`), so the rows of a level go in key order.
- **toppic `--mass-calib`** corrects `ms1_peak.mz` in place (`mz` is not a
  key column) and sets each `CONFIG` row's `MZMIN` / `MZMAX` from
  `ms1_peak` (`sqlite_mass_calib`; see `.claude/tools/toppic_output.md`).
- **Checks (toppic_test experiment 165, report 0165)** against main
  01d2db5a on st_1 and 480_F3_02, with and without `-D`: msalign and
  feature files equal (but the time line), every other table the same
  rows, `ms1_peak` row for row the conversion of main's file (levels from
  `PEAKS<k>`, matched by retention time and m/z), `PEAKS0..n` rebuilt from
  it with `COLOR` equal to main's, toppic's results equal.
- **Size and time on 480_F3_02 `-D`**: 894 MB instead of 2,516 (`ms1_peak`
  420 MB instead of 401 + 1,641); the 3D step 11 s instead of 16; peak
  memory unchanged (feature detection, 5.2 GiB; the 3D step 3.8 instead of
  3.6 GiB: level 0 is kept, and a peak carries its key). Without `-D` the
  file is 881 MB instead of 875.
- **Page fill**: SQLite balances the pages of a `WITHOUT ROWID` table (an
  index B-tree) on every split, also on appends in key order, so
  `ms1_peak` is written about 88% full; only a rowid table's appends fill
  their pages (writing the MS1 records in spectrum order does not help).
  A `VACUUM` packs it: 843 MB with `-D` and 839 without, in 3-4 s on
  480_F3_02. So TopFD and TopDIA end each database with
  `TopfdPara::closeSqlDb()`: `VACUUM` (its copy in a temporary file,
  `temp_store = FILE`, not in memory as `MzmlMsSqlWriter`'s setting would
  put it; a failure is a warning and leaves the file as written), then
  close. topfd `-T` does not (one spectrum); its database is closed by
  `createSqlDb` or the destructor as before.
