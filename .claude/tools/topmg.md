# TopMG notes

Facts to keep straight when changing topmg's graph alignment
(`src/search/graph`, `src/search/graphalign`, `src/stat/mcmc`,
`console/topmg_*`) and topmg_gui. toppic hosts the same identification
with `--topmg`. The notes on what topmg shares with toppic are in the
TopPIC files: the filters' mass list and the tolerances in
`toppic_search.md`, the PrSM tables, ProForma strings and the proteoform
classification in `toppic_output.md`.
Moved out of the repository's `CLAUDE.md` on 2026-10-12; the code-wide
conventions, the dependencies and the source layout are there.

## Graph alignment tolerance (`search/graphalign`, TopMG)

TopMG's graph alignment matches the mass differences (distances) of the
spectrum graph with those of the proteoform graph within
`GraphAlignMng::topmg_align_error_tolerance_` (Da; a constructor argument,
so both drivers must pass it), applied as `getIntTolerance()` in the
integer masses of `convert_ratio_` (274.335215 per Da: 0.1 Da is 28
units) by `graph_align.cpp` (distance matching and cutoffs),
`SpecGraph::compDistCover` and `graph_post_processor.cpp` (PTM masses this
close are merged). It is the argument `topmgAlignErrorTolerance`, default
0.1 Da (the former hard-coded `error_tolerance_`): topmg's hidden
`--align-error-tolerance` (visible in `-h` until 2026-10-08, when it was
moved to the hidden options and topmg_gui's field "Alignment error
tolerance (Da)" on the Advanced Parameters tab was removed, with its
`--align-error-tolerance` entry in `gui/util/command.cpp`, so the GUI leaves
topmg at 0.1 Da; the parameter block still prints it) and
toppic's hidden `--topmg-align-error-tolerance` for `--topmg pipeline` /
`topmg-only`, as `topmgShiftNumber` is `--num-shift` / `--topmg-num-shift`.
Both parsers refuse a value that is not a positive number. It is
independent of `-e` / `-E` (ppm): the precursor tolerance enters the graph
step only through the spectrum graph's zero and precursor PRM peaks
(`spec_graph_reader`) and the precursor refinement after an alignment.
The MCMC E-values (`MCMCMng::error_tolerance_`) keep their own hard-coded
0.1 Da. Facts to keep straight:

- **The post-processing keeps the variable PTM limit.**
  `GraphPostProcessor` explains each variable shift by a PTM combination of
  `mass_ptm_map_`, which holds combinations of up to `max_known_mods_` PTMs,
  a combination being left out when a mass within the tolerance is already
  there. Each shift takes its nearest combination; when these exceed
  `max_known_mods_` (`--var-ptm`) PTMs in all, `choosePtmCombs` chooses the
  combinations together: the smallest sum of mass differences with at most
  `max_known_mods_` PTMs (a dynamic programme over the shifts and the PTM
  number). Below 0.1 Da fewer combinations are left out, and the nearest
  ones gave proteoforms more variable PTMs than `--var-ptm` (6-9 of 5 at
  0.02 and 0.05 Da on the first 200 spectra of 480_F3_02; 4-5 of toppic's
  3), which overran `ptm_mass_vec2d` in `mcmc_dpr_processor.cpp` (indexed
  by the variable PTM number; a segmentation fault). Main 1f88b6e8 kept the
  limit by giving each shift in turn only the PTMs left, which labelled a
  70.06 Da shift of H2B1N as one Phospho (mass split 70.06); the joint
  choice gives Acetyl + 2 Methyl (toppic_test experiment 114). Wherever the
  nearest combinations keep the limit (every PrSM at 0.1 and 0.5 Da) they
  are the choice, so those outputs are unchanged.
- **TopMG's final tables are not repeatable with more than one thread.**
  `comp_pvalue_mcmc.cpp` shuffles with one global `std::mt19937 rng(0)`
  that all MCMC threads share, so the E-values differ between runs (1-2 of
  344 identical between two `-u 16` runs of main 0f938697, 110 and 113
  PrSMs); the intermediate files up to `topmg_graph_post` are repeatable
  and are what to compare (`-k`). toppic's `--topmg` stage with its default
  generating-function E-values is repeatable.

## Protein N-terminal forms in TopMG (`search/graph`, `search/graphalign`)

TopMG's graph alignment (topmg and toppic's `--topmg`, the same code)
reports the protein N-terminal forms of `-n` (`NONE`, `NME`,
`NME_ACETYLATION`, `M_ACETYLATION`) as toppic does: the "Protein N-terminal
form" column, the start residue, the `[Acetyl]-` prefix, the acetylation
not a variable PTM. Facts to keep straight:

- **Graph edges** (`ProteoAnno::anno`, only for the protein's N-terminal
  sequence and only the forms `prot_mod_util::allowMod` accepts for it):
  `M_ACETYLATION` is an acetylated-M edge at residue 0; `NME_ACETYLATION_X`
  is an edge at residue 0 that weighs the acetyl group alone (an empty
  amino acid with the Acetyl PTM: M removed, residue 1 acetylated). Both
  are `PROTEIN_VARIABLE` and take one of the alignment's `--var-ptm`
  layers. Being at residue 0, only an alignment that starts at the protein
  N-terminus can use them. The NME acetylation used to be an acetylated
  residue 1, open to an alignment that kept the M (`M(A)[Acetyl]...`).
  `NME` is the start at vertex 1; `is_nme_` (reset for every protein) says
  whether that start is the protein N-terminus (`getFirstDiagonal`), and
  with `-w` an alignment starts at vertex 1 only then.
- **The form is read back in `GraphPostProcessor`** (`chooseNTermForm`),
  before the variable shifts are named, from the start and the variable
  shift that begins at the first residue (the edge's layer change makes the
  first diagonal end at the start node, so that shift holds the edge's mass,
  possibly with PTMs of the next residues). Start 0: `NONE`, `M_ACETYLATION`
  (shift - acetyl) or `NME_ACETYLATION` (shift - acetyl + M, the
  proteoform then starts at residue 1: residues, fixed and other shifts
  moved by one); start 1: `NME` if allowed, else `NONE` (a truncation), or
  `NME_ACETYLATION` (shift - acetyl); later starts `NONE`. An acetylated
  form is taken when what it leaves (nothing, or the nearest PTM
  combination) explains the shift at least as well as the shift's own
  naming (error in Da within the tolerance, or no larger than the shift's
  own error up to twice the tolerance) with fewer variable PTMs (or as
  many and a smaller error). The form is built as
  `proteoform_factory::geneProtModProteoform` builds it: the `ProtMod`, the
  residue sequence's N-terminal mod (its mass), a `PROTEIN_VARIABLE` shift
  on residue 0 (the annotation). A mass left over (at most twice the
  tolerance) is dropped, so the proteoform mass can differ from the
  adjusted precursor mass by that much. Main used to name the acetylation
  as a variable `Acetyl` (adding Acetyl to the naming whenever an
  acetylated form was allowed; now only the `-i`/`-b` file's PTMs are
  names) and set `NME` for every start at residue 1.
- **Unexpected shifts are not read as forms**: an N-terminal acetylation
  that an alignment explains with an unexpected shift (topmg's default
  `--num-shift 1` allows one) stays one.
- **MCMC E-values** (`comp_pvalue_mcmc`): the residue sequence's N-terminal
  mod is added to the N-terminal fragment masses of the PrSM and of the
  random sequences (`n_term_shift_`); the candidate counts already have
  the forms (`CountTestNum`). TDGF (toppic `--topmg`) handles the forms as
  for toppic's own PrSMs.
- **Truncation N-terminal modifications** (proteoform N-terminal variable
  PTMs, Position `N-term` of the PTM file; toppic's since experiment 100,
  "Truncation N-terminal modifications" in `toppic_search.md`). The graph
  has no edge for
  them: `GraphAlign::getConsistentPairs` also matches each distance of the
  spectrum's zero vertex (the N-terminus) with proteoform distances m
  shorter (m: each N-term mass, `GraphAlignMng::n_term_int_masses_`, filled
  by `GraphAlignProcessor::process`) from the vertices where the mod is
  allowed (`GraphAlign::n_term_starts_`: protein residue 3 or later, an
  amino acid of the line), one variable PTM layer more; the spectrum's
  distance cover marks those distances (`SpecGraph::compDistCover`).
  `chooseNTermForm` reads the mod back from the variable shift that begins
  at the first residue of a proteoform starting at residue 3 or later (the
  same rule as the protein forms) and sets it with
  `Proteoform::setTruncNTermModPtr`. Variable shifts are named only with the
  side-chain PTMs (`readModTxt` list 2), in the post-processing and in the
  MCMC's combinations (an N-term line used to name side-chain shifts, and
  one name on an N-term and a side-chain line made `readPtmTxt` throw); the
  MCMC candidate number of a PrSM with the mod sums over the N-term masses;
  `toppic --topmg topmg-only` sets `nTermVarPtmTypeNum` (the TopMG stage's
  and the cluster correction's generating function multiply by T_N + 1).
  **Ties go to the path without the mod**, so that a protein N-terminal
  form wins a tie against it: `GraphAlign` keeps a flag per DP entry
  (`trunc_paths_`; a mod's consistent pair stores its cell c as -1 - c),
  and wherever two paths score the same (the best predecessor, a variable
  against an unexpected-shift path, the best shift node, the end cell, and
  equal matched fragment numbers over the PTM numbers in `geneResult`) the
  one without the mod is kept. Without N-term lines nothing changes.
  Results (toppic_test experiment 116, report 0116; 480_F3_02, topmg
  `-f C57 -d -t FDR -T FDR`, oxidation M, phospho STY, `Acetyl,*,N-term`):
  topmg reports 12 and 16 truncation N-terminal acetylations at 1% FDR,
  8 and 9 of toppic's 11 (DDX17 from residue 25, JUPI2 44, CLCB 69, TRIR
  28), `--topmg topmg-only` 39 with 9 of them (main: none); PrSMs 1,973
  and 2,022 against 2,005 (MCMC noise), no time cost. The modification
  competes with the protein forms: the alignment keeps the path with the
  most matched fragments whatever its PTM number, so an acetylated
  truncation plus PTMs can displace a protein form (24 PrSMs of the
  ported TopMG start elsewhere than toppic's, 11 where toppic has NME or
  NME_ACETYLATION). With the tie rule (fcba5f51) topmg reports 14 (9 of
  toppic's) and `--topmg topmg-only` 26 (the same 9): where toppic has a
  protein form 6 instead of 11, the other unconfirmed ones 11 instead of
  19 (H2B2E, H2B2F, CLCB now toppic's forms with its E-values). The 6 left
  are not ties: the truncation matches at least as many fragments as the
  best path main finds, or (YBOX3, one spectrum) has the higher DP score in
  the protein form's layer, so the form is dropped before its fragments
  are counted; giving the mod's paths their own DP layers would let both
  reach the fragment count.
- **Results** (toppic_test experiment 115, report 0115; 480_F3_02, topmg
  `-f C57 -d -t FDR -T FDR`, oxidation M, acetyl K, phospho STY): on the
  spectra topmg and toppic both identify with the same protein and toppic
  starts at residue 1-2, the same form and start in 89.9% (main 21.7%;
  `--topmg topmg-only` 81.6%, main 18.5%); the final topmg table has 299
  NME_ACETYLATION and 35 M_ACETYLATION PrSMs and no variable Acetyl where
  an acetylated form fits (main: 188 such PrSMs, 140 Acetyls on ranges
  without K). topmg 2,009 and 1,995 PrSMs against 1,854 (the acetylated
  forms' MCMC E-values, a median 10^0.63 better); `--topmg topmg-only`
  2,976 against 2,979 (TDGF: one PTM fewer, x1/4). No time cost. The
  remaining disagreements with toppic are alignments: histone H2B starts,
  isobaric alternatives (acetyl-S N-terminus against a C-terminal E;
  acetylation + K adduct against Phospho), an acetyl put on a K further in.

## PTMs written in FASTA sequences in TopMG (`search/graph`, `search/graphalign`)

A FASTA sequence can carry PTMs in the extended format, a residue followed
by the PTM's abbreviation of `res/base_data/ptm_base.xml` in brackets
(`M[Oxidation]`, `C[Carbamidomethyl]`; the names of earlier versions, such
as `C[Carbamidomethylation]`, are no longer read, see "ProForma strings"
in `toppic_output.md`).
`FastaSeq::compAcidPtmPairVec` takes a name up to the `]` that closes its
`[`, counting nested pairs, so the 51 Unimod names with brackets can be
written (`E[Cation:Mg[II]]`; it used to stop at the first `]`); an
unpaired `[` stops with "Unpaired '[' after residue <n> in the sequence
of <protein>".
`FastaSeq` reads them into residue/PTM pairs (`getAcidPtmPairVec`), and
`proteoform_factory::geneDbProteoformPtr` builds the modified residues with
an `INPUT` mass shift on each (shown in the Proteoform column as
`(M)[Oxidation]` and in the Database protein sequence; not in the Fixed
PTMs column or the variable PTM counts). Since branch `topmg_fasta_ptm` TopMG's graph alignment (topmg and
toppic `--topmg`) treats them as fixed modifications. Facts to keep
straight:

- **The mass graph is built from the pairs** (`ProteoAnno::anno` takes a
  `StringPairVec`). It used to take the raw string, whose
  `convertStrToResiduePtrVec` dropped the brackets and lower-case letters
  and kept the rest: `M[Oxidation]` became M and O (pyrrolysine), the graph
  had more vertices than the proteoform had residues, and topmg and toppic
  `--topmg` ended with a segmentation fault at the start of the graph
  alignment (toppic_test experiment 122).
- **One edge per modified residue.** A residue with a PTM, of the FASTA or
  of `-f`, is marked `FIXED`, and `getMassGraphPtr` gives a position with a
  `FIXED` residue that edge alone (and `PROTEIN_VARIABLE` edges, below): no
  variable PTM edge. A first residue with a FASTA PTM gets no
  `M_ACETYLATION` edge, which would modify it, and
  `GraphPostProcessor::chooseNTermForm` does not take that form for it. The
  forms that remove it stay: NME (a start, not an edge) and
  `NME_ACETYLATION`, whose edge in place of the M is the position's
  `PROTEIN_VARIABLE` edge; the post-processing drops the FASTA PTM of the
  removed M (an `INPUT` or `FIXED` shift on residues before the form's
  start). A first version took every acetylated form away from such a
  residue: with M1 written oxidized, YBOX3, FUBP2 and H1.5 lost their
  NME_ACETYLATION PrSMs on 480_F3_02 (experiment 122). The truncation
  N-terminal modifications (not edges either) stay open, as for a fixed
  modification's residue. A `-f` modification is not applied to a residue
  with a FASTA PTM (`applyFixedMod` compares with the unmodified residue).
- **The post-processing keeps the `INPUT` shifts** (it rebuilds each
  proteoform from its `FIXED`, `VARIABLE` and `UNEXPECTED` shifts, which
  dropped them from the reported proteoforms).
- **Residues are registered before the threads start.**
  `fasta_util::dbPreprocess` ends with `addPtmResidues`, which adds the
  residues of the database's FASTA PTMs to `ResidueBase` (a bracket is
  looked for first; a database without PTMs is read once more as text and
  nothing else happens). `ResidueBase::getBaseResiduePtr` appends a new
  residue to a shared vector without a lock, which the search threads did
  for a PTM not in the `-i`/`-b` file. An unknown PTM name now throws there.
- **E-values**: `count_util::updateResidueCounts` (the residue frequencies
  of `CountTestNum`: TDGF and the MCMC's candidate numbers) adds a residue
  that is not in the list, so that a residue with a FASTA PTM is counted as
  its own residue, as a fixed modification's residue is; it used to be
  left out (toppic `--topmg` E-values a median 10^0.1 larger with every C
  written as `C[Carbamidomethylation]` than with `-f C57`). The MCMC's
  random residues (`mcmc_mass_table_util`) take the fixed modifications of
  `-f` only.
- **Protein hits** (`prsm/search_fasta_match`, the other proteins of a PrSM
  in the TSV and SQLite tables of toppic and topmg) index the residues, one
  letter each, and a hit must have the same FASTA PTMs. The raw string was
  indexed and cut at the proteoform's residue positions: with brackets
  before the start the cut was a piece of a PTM name and the next residues,
  which matched thousands of proteins (350,183 rows of other proteins in
  topmg's PrSM table of 480_F3_02 with every C written as
  `C[Carbamidomethylation]`, against 269 with `-f C57`). A database
  without PTMs is indexed and searched as before.
- The PTM localization of toppic (`-l`, `stat/local/local_proteoform`)
  still reads the raw string (`getRawSeq`); topmg has no such stage.
- Results (toppic_test experiment 122, report 0122; 480_F3_02, `-d -t FDR
  -T FDR`, oxidation M, phospho STY, `Acetyl,*,N-term`): every C written
  as `C[Carbamidomethylation]` without `-f` gives `-f C57`'s filter,
  alignment and post-processing files (topmg) and all of `--topmg`'s files
  with its E-values and tables (but the Fixed PTMs and Database protein
  sequence columns), with the same index (built in memory; the shared
  index files of toppic_test's `data/`, built on 2026-09-28, give 70 of
  60,135 one-PTM filter candidates differently). Without FASTA PTMs the
  files are main's. With the 249 methionines main reports oxidized written
  as `M[Oxidation]`, 86-90% of the PrSMs that had the oxidation are the
  same proteoform with it from the FASTA (`--topmg` E-values a median
  10^-0.48 smaller), and the unoxidized forms over those sites take an
  unexpected shift or are lost (`--topmg` 2,841 PrSMs against 2,959).

## Precursor refinement of TopMG's graph alignment (`search/graphalign/graph_align`, `search/diag/diagonal_util`)

Branch `topmg_prec_matched_ions` (toppic_test experiment 143, report 0143;
merged into main as 5c9e2e3d on 2026-10-07, after `filter_pairs`; on the first
500 spectra of 480_F3_02 toppic in setting d writes the same files as main
20d5895b, apart from the tables' start and end times). `GraphAlign::geneResult(s, m)`, shared by topmg and toppic `--topmg`,
calls `refinePrecursorByMatchedIons` like `PsAlign` and `PsAlignMultiType`;
`refinePrecursorAndHeaderShift`, its two helpers and
`GraphAlignMng::refine_prec_step_width_` are removed. Facts to keep straight:

- **Where TopMG's precursor error enters.** The spectrum graph's vertices are
  PRM peaks of the observed precursor P (C-terminal ones P - (m - C shift)),
  matched by distance within the alignment error tolerance (0.1 Da); every
  diagonal's C-terminal shift is set from P (`geneHeaders`). The first
  diagonal's shift comes from the zero vertex, the last one's from the
  precursor vertex, an internal one's is the mean of its nodes' shifts.
- **Single diagonals.** An alignment without unexpected shifts and variable
  PTMs is one diagonal, anchored at the zero vertex: its N-terminal ions do not
  move with P. `getDepIonPairs` takes its C-terminal ions as the dependent ions
  and the refinement moves its C-terminal shift (be5fb6ff). With two or more
  diagonals the rule is unchanged (the C-terminal ions of every diagonal but
  the last, the N-terminal ions of the last). toppic's searches never pass one
  diagonal: the one-PTM search asks `PsAlign` for one shift, the multiple-PTM
  search for two or more.
- **Downstream.** The PrSM's adjusted precursor is the refined one;
  `GraphPostProcessor` and the MCMC take it, the generating function of toppic
  `--topmg` works on the observed precursor's PRM spectra (its E-values did not
  change when only the single diagonals' precursor moved).
- **Results on 480_F3_02** (topfd main 7b5112c2, `-f C57 -d -t FDR -T FDR`,
  Oxidation M, Phospho STY, Acetyl N-term; toppic `--topmg topmg-only`):
    - The old method put 37.6% of 9,885 target graph alignments at its +2T
      fallback (decoys 50.9%; 60% of the alignments with fewer than 5
      fragments, 51% of single diagonals) and 44.3% outside [-T, T]; 281 of
      2,981 reported PrSMs at +2T. Now none outside T.
    - Same alignments: 1,024 of 9,238 targets gain matched fragments, 170 lose.
    - Reported PrSMs without unexpected shifts: adjusted precursor within 1 ppm
      of the proteoform mass for 96.1% (old 54.9%), > 10 ppm off 0.2% (6.0%).
    - Every precursor moved by up to +-5 ppm: 82.1% of the reported PrSMs get
      exactly the unmoved run's adjusted precursor (old 0%); strong alignments
      with one diagonal 98.9%, two 90.5%, three or more 35.0%.
    - 1% FDR: 2,968 / 1,748 / 403 PrSMs / proteoforms / proteins against
      2,981 / 1,754 / 398; more targets at every fixed E-value (E <= 0.01:
      2,752 against 2,736) and 4 / 10 more decoys at E <= 10 / 100, so the
      cutoff fell from 55.3 to 24. On the moved precursors 2,957 / 1,744 / 405
      against 2,954 / 1,738 / 402. Same running time.
    - topmg (MCMC E-values, two runs each): 1,979 / 1,965 PrSMs against
      2,021 / 2,002, proteoforms 1,212 / 1,217 against 1,239 / 1,219,
      proteins 390 / 390 against 394 / 386. The refined decoys gain fragments
      as the targets do (the old fallback hit half of the decoy alignments):
      decoys at E <= 0.1 114 / 121 against 93 / 98. Unchanged PrSMs' MCMC
      E-values scatter as between two runs of one binary (sd 0.9 in log10).
      MCMC +1.6 min (15,423 post-processed PrSMs against 14,818).
- **Known effects:**
    - Three or more diagonals recover the precursor for only 35%: an internal
      diagonal's mean shift mixes N-terminal and C-terminal PRM vertices, so it
      carries part of P's error.
    - A variable-PTM diagonal that matched no ion was dropped by
      `refineHeadersBgnEnd` and its mass merged into the unexpected shift; when
      the refined precursor lets one ion match, the diagonal stays and the
      unexpected shift can exceed `max_ptm_mass_` (500 Da), so the PrSM is
      rejected (spectrum 712 of 480, THOC4, E 2.5e-8 with the old method).

## SQLite output (`--no-sql`)

topmg writes its PrSM tables (`prsm_sql_output::write`) into TopFD's
database of the spectrum file (`prsm_post_mass_match::sqlFileName`, the
`_post_ms2` of toppic's post-matched files read as `_ms2`) when it exists;
it has no database of its own. The visible `-z`, `--no-sql` (branch
`no_sql_output`, merged into main as 4cdaa7a1; toppic_test experiment 164;
the same option as toppic's,
see `toppic_output.md`) skips that write, so topmg writes or changes no
SQLite file; the parameter block prints "SQLite database output: No" only
when it is set.

