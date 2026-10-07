# Runtime resources

The tools of the TopPIC Suite load these files at run time. They look for
this directory as `res` next to the executable first, then as `../res` (so
the tools can be run from the repository's `bin/` without installing), then
in the installed shared directory (`share/toppic` on Linux). `make install`
copies the whole directory there.

Two kinds of file are stored with [Git LFS](https://git-lfs.com/) rather than
in the normal Git history: the `*.onnx` models and
`base_data/theo_patt.txt`. Without Git LFS a clone leaves small pointer
files in their place and the tools fail to load them; run `git lfs pull` to
fetch the real files.

## `base_data/`

XML tables of the chemical and analytical building blocks, each read once at
start-up by the matching `*_base` class in `src/common/base`:

| File | Contents |
|---|---|
| `amino_acid_base.xml` | The amino acids as residues (the amino acid minus H2O): one-letter and three-letter codes, composition, monoisotopic and average masses. The 20 standard ones, U (selenocysteine, monoisotopic mass with ⁸⁰Se, as in Unimod) and O (pyrrolysine, the IUPAC-IUBMB and UniProt letter). Monoisotopic masses are sums of the AME2020 isotope masses (¹H, ¹²C, ¹⁴N, ¹⁶O, ³²S, ⁸⁰Se), written exactly with all their decimals. Average masses are sums of IUPAC's 2021 abridged standard atomic weights (H 1.0080, C 12.011, N 14.007, O 15.999, S 32.06, Se 78.971), written exactly with all their decimals. The frequencies are those of the UniProt human reference proteome UP000005640, release 2026_03 (20,652 proteins, 11,456,352 residues; toppic_test experiment 108). No code reads the average masses or the frequencies: TopPIC's E-values and candidate counts use the residue frequencies of the search database (`stat/count/count_test_num.cpp`). |
| `ptm_base.xml` | Post-translational modifications: name, abbreviation, monoisotopic mass and Unimod id. It holds every Unimod modification, written by toppic_test's `tests/202610/scripts/proforma_output/gen_ptm_base.py` (experiment 128) from `unimod.obo` (release of 17 February 2026: 1,560 modifications) in Unimod id order: the abbreviation, the name TopPIC writes and reads in modification files, is the Unimod (PSI-MS) name; the name is the Unimod full name; the mass is the exact sum of the AME2020 nuclide masses (CIAAW table; H 1H, C 12C, Se 80Se, ...) over the Unimod composition, the bricks (Hex, HexNAc, ...) expanded with Unimod's element lists, written with all its decimals. Before them `No PTM` (the empty PTM, found by that name), the file's only entry that is not Unimod's. Experiment 128 renamed the earlier entries to their Unimod names and dropped the old names (Carbamidomethylation, now Carbamidomethyl, and seven more); they are not kept as aliases. The ProForma column writes a PTM by its Unimod name and id ([Cation:Na|UNIMOD:30]; `PtmBase::getUnimodName` looks the name up by the id). To follow a new Unimod release, rerun the script on the new `unimod.obo`; do not edit the file by hand. |
| `residue_base.xml` | Residues, i.e. amino acids with an attached PTM, used to build proteoform sequences. |
| `mod_base.xml` | Modifications as residue-to-residue changes (the form used by the fixed / variable modification options of TopPIC and TopMG). |
| `prot_mod_base.xml` | Protein N-terminal forms: methionine excision, N-terminal acetylation and their combinations. |
| `trunc_base.xml` | N- and C-terminal truncations allowed for proteoforms. |
| `activation_base.xml` | Fragmentation methods (CID, HCD, ETD, ETCID, ETHCD, UVPD, MPD) with the ion types each produces; ETCID (ETciD) and ETHCD (EThcD) have the ion types of ETD. |
| `ion_type_base.xml` | Fragment ion types (a, b, c, x, y, z, ...) and their mass shifts. The shifts are sums of the AME2020 isotope masses (¹H, ¹²C, ¹⁴N, ¹⁶O), written exactly with all their decimals. |
| `neutral_loss_base.xml` | Neutral losses considered for fragment ions. Their masses are sums of the AME2020 isotope masses (water 2 ¹H + ¹⁶O, ammonia ¹⁴N + 3 ¹H), written exactly with all their decimals, as in `mass_constant.hpp`. |
| `support_peak_type_base.xml` | Peak types that count as support for a fragment in scoring. |

Data tables:

| File | Contents |
|---|---|
| `theo_patt.txt` | Theoretical isotopic distributions, tabulated by mass, used by TopFD's envelope detection (`EnvBase`): 11,000 entries, one per 10 Da of monoisotopic mass (entry i for 10 i Da; entry 0 repeats entry 1). Each is the human averagine of the UniProt human reference proteome (UP000005640, release 2026_03: C 4.9202, H 7.7647, N 1.3745, O 1.4842, S 0.0443 per residue of the 20 standard amino acids) scaled to the mass with each of the five atom numbers rounded (no hydrogen filling the mass, so an entry's exact mass, in its formula line, lies up to about 35 Da from 10 i Da), computed by the emass method with AME2020 isotope masses and IUPAC representative abundances (toppic_test experiments 107, 108 and 111). Every peak from the monoisotopic one is listed (`Env::Env` takes the first listed peak as the monoisotopic one), up to the last of at least 5e-9 of the apex; intensities in % of the apex. The previous table was emass output with emass's old abundances, another composition and no monoisotopic peak from 36 kDa. `common/emass` reproduces the table byte for byte: `make theo_patt_gen && bin/theo_patt_gen -o res/base_data/theo_patt.txt` (toppic_test experiment 112). Stored with Git LFS (about 21 MB). |
| `mass_table.txt` | Table of residue combinations by mass, used by TopMG's TopMCMC p-value estimation. |
| `env_rescore_para.txt` | Parameters of a former envelope re-scoring model; not read by the tools. |

## `envcnn_models/`

`envcnn_two_block.onnx`: the EnvCNN neural network that scores isotopic
envelopes in TopFD and TopDIA (and TopPIC's post mass matching), run through
the ONNX Runtime. Stored with Git LFS.

## `ecscore_models/`

`ecscore_seven_attr.onnx`: the ECScore neural network that scores proteoform
features detected in the LC-MS map by TopFD and TopDIA. Stored with Git LFS.
