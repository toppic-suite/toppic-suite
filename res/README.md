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
| `ptm_base.xml` | Post-translational modifications: name, abbreviation, monoisotopic mass and Unimod id. Monoisotopic masses are sums of the AME2020 isotope masses (¹H, ¹²C, ¹⁴N, ¹⁶O, ³¹P, ³²S, ¹²⁷I) over the composition of the Unimod record, written exactly with all their decimals; the entries without a Unimod id have Carbamidomethylation H(3)C(2)NO, Decarboxy C(-1)O(-2), Citrullination H(-1)N(-1)O and TPQ H(-2)O(2). |
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
the ONNX Runtime; the file the tools load. It is the dynamically quantized
form of `envcnn_two_block_full.onnx` (22 MB): the three fully connected
layers have int8 weights with one scale per output channel, and their inputs
are quantized to uint8 at run time with one scale per envelope (the row
maximum / 255; the inputs are ReLU outputs, >= 0), so that an envelope's
score does not depend on the other envelopes of its batch. ONNX Runtime's
`quantize_dynamic` scales a whole batch by its largest value, which moved the
deconvolution's scores by 0.06 on average; this model moves them by 0.002.
The convolution layers stay float. Made by `quantize.py` of toppic_test
experiment 131 (variant `fc_row`); with it TopFD takes 0.70 of the float
model's time on 480_F3_02 (toppic_test report 0131). Stored with Git LFS.

`envcnn_two_block_full.onnx`: the float EnvCNN model (PyTorch export, opset
11, 87 MB), the source of the quantized model; not loaded by the tools.
Stored with Git LFS.

## `ecscore_models/`

`ecscore_seven_attr.onnx`: the ECScore neural network that scores proteoform
features detected in the LC-MS map by TopFD and TopDIA; the file the tools
load. It is the dynamically quantized form of `ecscore_seven_attr_full.onnx`
(ONNX Runtime's `quantize_dynamic` after each Gemm was rewritten as MatMul +
Add: int8 weights with one scale per tensor, the inputs quantized to uint8
at run time; 122 kB). Made by `quantize.py` of toppic_test experiment 131
(variant `fc`). Its scores differ from the float model's by 0.027 on average
(TopFD reported 6% fewer features on 480_F3_02, toppic_test report 0131).
Stored with Git LFS.

`ecscore_seven_attr_full.onnx`: the float ECScore model (PyTorch export,
opset 11, 454 kB), the source of the quantized model; not loaded by the
tools. Stored with Git LFS.
