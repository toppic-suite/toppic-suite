# Proteoform annotation in TopPIC and TopMG results

TopPIC and TopMG report every identified proteoform as an annotated
sequence string, for example

```text
MKVRASVK(KL(C)[Carbamidomethylation]RN(C)[Carbamidomethylation]K)[-116.0905]IVKRDGVIRVI(C)[Carbamidomethylation]SAEPKHKQRQG
```

This page explains how to read the string: what the parentheses, brackets,
nesting and leading `[...]-` mean, which modifications are shown by name
and which as a mass, and where the names and masses come from.

## 1. Where the annotation appears

- The `Proteoform` column of the PrSM, proteoform and protein tables
  (`*_prsm.tsv`, `*_proteoform.tsv`, `*_protein.tsv` and their `_single`
  versions). The residues flanking the proteoform in the protein are in the
  separate `Previous amino acid` and `Next amino acid` columns; each is `-`
  when the proteoform reaches that end of the protein.
- The `proteo_match_seq` element of the result XML files and the
  `proteoform` column of the identification tables in the SQLite database,
  which hold the same string.
- Some other outputs (visualization tools and older TopPIC releases) join
  the three columns into one string, `previous.proteoform.next`. A dot with
  nothing beside it then means the proteoform reaches the protein terminus:
  `.MKVR...QG.` is a complete protein.

## 2. Reading the string

The plain letters are the amino acids of the proteoform, from its first
residue (`First residue` column) to its last (`Last residue`), as in the
database sequence. Modifications are annotated on top of them:

| Notation | Meaning |
|---|---|
| `(X)[Name]` | Residue `X` carries the modification `Name`. |
| `(XYZ...)[label]` | The mass shift or modification `label` is located somewhere on the residues in parentheses; the fragment ions do not determine which one. |
| `(...)[+123.4567]`, `(...)[-116.0905]` | An unexpected mass shift, given as a mass in Da (see section 4). |
| `[Acetyl]-M...` | Acetylation of the protein N-terminus, written as a prefix before the first residue. |

Every annotation is written as a pair of parentheses around the residues it
applies to, immediately followed by its label in brackets. The label is
either a modification name or a signed mass.

**Nested parentheses.** Annotations are added in layers: first the
modifications present in the input, then fixed modifications, then protein
N-terminal modifications, then variable modifications, and finally the
unexpected mass shifts. When an unexpected shift spans residues that carry a
named modification, the named modification keeps its own `(X)[Name]` and the
shift's parentheses enclose it, giving nested parentheses. In the example
above,

```text
(KL(C)[Carbamidomethylation]RN(C)[Carbamidomethylation]K)[-116.0905]
```

means: both cysteines are carbamidomethylated (a fixed modification), and an
additional unexpected mass shift of -116.0905 Da is located somewhere on the
seven residues `KLCRNCK`. The inner annotations are per residue; the outer
one is a range.

**Position ranges.** The residues inside an unexpected shift's parentheses
are those between the nearest sequence positions on either side that
matched fragment ions of the spectrum. A single residue in parentheses,
`(A)[-0.0639]`, means the shift is pinned to that residue; a long range
means the spectrum has no fragment ions inside it.

## 3. Which modifications are shown by name

A modification is written by name when the search knew it in advance, and
as a mass when it did not:

| Kind | Given by | In the string | Columns |
|---|---|---|---|
| Fixed modification | TopPIC/TopMG `-f` (`C57`, `C58` or a modification file) | `(C)[Carbamidomethylation]` on every modified residue | `Fixed PTMs` |
| Protein N-terminal form | `-n` (`NME`, `NME_ACETYLATION`, `M_ACETYLATION`) | Acetylation as the prefix `[Acetyl]-`; methionine excision is not annotated, the proteoform simply starts at residue 2 | `Protein N-terminal form` |
| Variable modification | TopPIC `-b`, TopMG `-i` (a modification file) | `(S)[Phospho]` on the modified residue | `# Variable PTMs`, `Variable PTMs` |
| Unexpected mass shift | TopPIC `-s`, `-m`, `-M`; TopMG `-s` | `(...)[+79.9663]`, a signed mass over a residue range | `# Unexpected modifications`, `Unexpected modifications` |
| Characterized unexpected shift | TopPIC `-B` (a common modification file) with PTM characterization | The shift is replaced by the matching modification name, `(S)[Phospho]`, over the residues the characterization allows | `Variable PTMs`, `MIScore` |

The `Fixed PTMs`, `Unexpected modifications` and `Variable PTMs` columns
list the same annotations as `name:[position]` or `name:[first-last]`
(1-based positions in the proteoform), separated by `;`; for unexpected
shifts the "name" is the mass, for example `-0.0467:[7-33]`.

**Where the names and masses come from.** The name in brackets is the
modification's abbreviation, and its mass is what TopPIC/TopMG add to the
residue mass:

- `C57` and `C58` are the built-in modifications `Carbamidomethylation`
  (57.021464 Da) and `Carboxymethyl` (58.005479 Da) on cysteine, and
  N-terminal acetylation is the built-in `Acetyl` (42.010565 Da). Their
  definitions are in the resource files `res/base_data/ptm_base.xml`,
  `mod_base.xml` and `prot_mod_base.xml`.
- A modification file (`-f <file>`, `-b`, `-B`, `-i`) supplies a name and a
  mass per line (`Name,Mass,Residues,Position,UnimodID`, see the TopPIC
  manual, section 5). The file is **not** the only source of the values,
  though: when the name is one of the abbreviations in
  `res/base_data/ptm_base.xml` (`Phospho`, `Acetyl`, `Oxidation`, `Methyl`,
  `Dimethyl`, `Trimethyl`, `Carbamidomethylation`, ...), the built-in entry
  is used and the mass in the file is ignored. Only a name that is not in
  that list is created from the file with the file's mass.
  So a line `Phospho,79.97,STY,any,21` in `common_mods.txt` still uses
  79.966331 Da, and a differently spelled name, say `Phosphorylation`,
  would use the mass from the file. To give a modification a different
  mass, use a name that is not in `ptm_base.xml`.

## 4. The mass in brackets

An unexpected shift is written as a mass in Da with four decimals and an
explicit sign (`+` for positive shifts). It is the part of the observed
proteoform mass that the sequence with all its **named** modifications does
not explain: fixed, N-terminal and variable modifications are already
included in the residue masses, so the number never contains them. In the
example, the -116.0905 Da is in addition to the two carbamidomethylations
inside the parentheses.

Shifts close to zero (`-0.0467`, `+0.9477`) are usually not modifications
but the precursor mass error or a one-Dalton isotope error of the
deconvolution. Shifts close to a known modification mass (`+79.9663`,
`+42.0106`, `+15.9949`) can be turned into names by running TopPIC with a
common modification file (`-B`); the `MIScore` column then gives the
confidence of the assigned position.

The `Proteoform mass` column is the mass of the annotated proteoform,
including all named modifications and the mass shifts, and
`Adjusted precursor mass` is the precursor mass it was matched to.

## 5. Examples

The first two strings are from TopPIC result tables; the others are
constructed to show N-terminal acetylation and fixed and variable
modifications (`Previous amino acid` | `Proteoform` | `Next amino acid`):

```text
- | MKRTFQ(PSVLKRNRSHGFRARMATKNGRQVLAR)[-0.0467]RRAKGRARLTVSK | -
```
A complete protein (flanked by `-` on both sides) with a -0.0467 Da shift
located somewhere on residues 7 to 33 (`Unexpected modifications` =
`-0.0467:[7-33]`): a precursor mass error, not a modification.

```text
M | SLSTEATAKIVSEFGRDANDTGSTDVQVA(L)[+97.9219]LTAQINHLQGHF... | -
```
The initiator methionine is removed (`Protein N-terminal form` = `NME`,
`First residue` = 2, previous amino acid `M`), and a +97.9219 Da shift is
pinned to the leucine at position 30 of the proteoform.

```text
M | [Acetyl]-SLSTEATAKIVSEFGRDANDTGSTDVQVALLTAQINHLQGHF... | -
```
Methionine excision followed by N-terminal acetylation
(`Protein N-terminal form` = `NME_ACETYLATION`).

```text
K | (C)[Carbamidomethylation]SAEPKHKQ(S)[Phospho]RQG | -
```
A C-terminal fragment of a protein searched with `-f C57` and a variable
modification file containing `Phospho`: the cysteine is carbamidomethylated
(`Fixed PTMs` = `Carbamidomethylation:[1]`) and the serine phosphorylated
(`Variable PTMs` = `Phospho:[10]`).

```text
MKVRASVK(KL(C)[Carbamidomethylation]RN(C)[Carbamidomethylation]K)[-116.0905]IVKRDGVIRVI(C)[Carbamidomethylation]SAEPKHKQRQG
```
Three carbamidomethylated cysteines (fixed, `-f C57`) and one unexpected
shift of -116.0905 Da somewhere on `KLCRNCK`; the shift's parentheses
enclose two of the fixed modifications.
