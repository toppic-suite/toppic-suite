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
| `amino_acid_base.xml` | The amino acids: one-letter and three-letter codes, composition, monoisotopic and average masses. |
| `ptm_base.xml` | Post-translational modifications: name, abbreviation, monoisotopic mass and Unimod id. |
| `residue_base.xml` | Residues, i.e. amino acids with an attached PTM, used to build proteoform sequences. |
| `mod_base.xml` | Modifications as residue-to-residue changes (the form used by the fixed / variable modification options of TopPIC and TopMG). |
| `prot_mod_base.xml` | Protein N-terminal forms: methionine excision, N-terminal acetylation and their combinations. |
| `trunc_base.xml` | N- and C-terminal truncations allowed for proteoforms. |
| `activation_base.xml` | Fragmentation methods (CID, HCD, ETD, UVPD, MPD) with the ion types each produces. |
| `ion_type_base.xml` | Fragment ion types (a, b, c, x, y, z, ...) and their mass shifts. |
| `neutral_loss_base.xml` | Neutral losses considered for fragment ions. |
| `support_peak_type_base.xml` | Peak types that count as support for a fragment in scoring. |

Data tables:

| File | Contents |
|---|---|
| `theo_patt.txt` | Theoretical isotopic distributions, tabulated by mass, used by TopFD's envelope detection (`EnvBase`). Stored with Git LFS (about 18 MB). |
| `mass_table.txt` | Table of residue combinations by mass, used by TopMG's TopMCMC p-value estimation. |
| `unimod_ptm.xml` | The Unimod PTM list in the `ptm_base.xml` format, kept for reference; not read by the tools. |
| `env_rescore_para.txt` | Parameters of a former envelope re-scoring model; not read by the tools. |

## `envcnn_models/`

`envcnn_two_block.onnx`: the EnvCNN neural network that scores isotopic
envelopes in TopFD and TopDIA (and TopPIC's post mass matching), run through
the ONNX Runtime. Stored with Git LFS.

## `ecscore_models/`

`ecscore_seven_attr.onnx`: the ECScore neural network that scores proteoform
features detected in the LC-MS map by TopFD and TopDIA. Stored with Git LFS.
