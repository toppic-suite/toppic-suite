# TopPIC notes: tolerances, filters and searches

Facts to keep straight when changing how toppic finds its candidate
PrSMs: the precursor and fragment tolerances, the filters and their
index files (written by topindex), the tag filter, the zero-shift,
one-shift, multiple-shift and variable PTM searches and their
precursor refinements (`src/para`, `src/filter`, `src/search`,
`src/ms/factory`). topmg uses the same tolerances and filters' mass list.
One of three TopPIC files (the others: `toppic_evalue.md` and `toppic_output.md`), split on 2026-10-12
from `toppic.md`, which had moved out of the repository's `CLAUDE.md`;
the code-wide conventions and the source layout are there.

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
