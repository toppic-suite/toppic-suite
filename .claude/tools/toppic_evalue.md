# TopPIC notes: E-values and FDRs

Facts to keep straight when changing how toppic scores its PrSMs: the
generating-function E-values (`src/stat/tdgf`), the deamidation
selector, and the rescoring of the PrSMs before the FDRs
(`src/stat/rescore`, `prsm_fdr`).
One of three TopPIC files (the others: `toppic_search.md` and `toppic_output.md`), split on 2026-10-12
from `toppic.md`, which had moved out of the repository's `CLAUDE.md`;
the code-wide conventions and the source layout are there.

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
