# Changelog

All notable changes to the Martini 3 Metabolome repository are documented in this file.

## [Unreleased] - Resolving Issue #4: Molecular Data Integrity

Reference: [GitHub Issue #4](https://github.com/Martini-Force-Field-Initiative/M3-Metabolome/issues/4)  
Dataset: [Zenodo 18622430 (v1)](https://zenodo.org/records/18622430)

### 1. Corrected Swapped L-Glutamine (LGLN) and L-Glutamate (LGLU) Entries
- **Problem**: In `misc/database.csv`, the entries for `LGLN` (L-Glutamine) and `LGLU` (L-Glutamate) had swapped `Formula`, `Canonical SMILES`, `MW`, and computational logP values (`iLOGP`, `XLOGP3`, `WLOGP`, `MLOGP`, `Silicos-IT Log P`, `Consensus Log P`).
- **Correction**: Swapped these columns back to their appropriate rows:
  - `LGLN`: Formula `C5H10N2O3`, Canonical SMILES `NC(=O)CC[C@H]([NH3+])C(=O)[O-]`, MW `146.14`, WLOGP `-3.39`, MLOGP `-5.05`.
  - `LGLU`: Formula `C5H8NO4-`, Canonical SMILES `[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]`, MW `146.12`, WLOGP `-4.12`, MLOGP `-4.64`.
- **Validation**: Verified molecular identities, atom connectivity, formal charges, and stereocenters using RDKit and SwissADME.

### 2. Corrected Choline (CHT) Molecular Identity
- **Problem**: In `misc/database.csv`, row `CHT` accidentally contained carnosine metadata (Formula `C9H14N4O3`, SMILES `[NH3+]CCC(=O)N[C@H](C(=O)[O-])Cc1nc[nH]c1`, MW `226.23`).
- **Correction**:
  - Formula: `C5H14NO+`
  - Canonical SMILES: `C[N+](C)(C)CCO`
  - MW: `104.17`
  - Computed logP values updated via SwissADME: iLOGP `-3.43`, XLOGP3 `-0.40`, WLOGP `-0.32`, MLOGP `-3.46`, Silicos-IT `-3.00`, Consensus `-2.12`.
- **Validation**: Consistent with `itps/molecules/CHT.itp` (ALC1 mass 36, NC3 mass 72, +1.0 charge), `mappings/map_files/CHT.map`, and `gros/CHT.gro`.

### 3. Specified Explicit L-Stereochemistry for Amino Acids & L-Lactate
- **Problem**: In `misc/database.csv`, SMILES entries for `LALA`, `LASP`, `LCYS`, `LGLN`, `LGLU`, `LILE`, `LSER`, `LVAL`, and `LAC` lacked stereochemical specification. `LHIS` had inverted (R) stereochemistry in its SMILES.
- **Correction**: Updated Canonical SMILES to explicitly describe verified L-stereochemistry:
  - `LALA`: `C[C@H]([NH3+])C(=O)[O-]` -> (S)
  - `LASP`: `[NH3+][C@@H](CC(=O)[O-])C(=O)[O-]` -> (S)
  - `LCYS`: `[NH3+][C@@H](CS)C(=O)[O-]` -> (R)
  - `LGLN`: `NC(=O)CC[C@H]([NH3+])C(=O)[O-]` -> (S)
  - `LGLU`: `[NH3+][C@@H](CCC(=O)[O-])C(=O)[O-]` -> (S)
  - `LHIS`: `[NH3+][C@@H](Cc1c[nH]c[nH+]1)C(=O)[O-]` -> (S)
  - `LILE`: `CC[C@H](C)[C@H]([NH3+])C(=O)[O-]` -> (2S, 3S)
  - `LSER`: `[NH3+][C@@H](CO)C(=O)[O-]` -> (S)
  - `LVAL`: `CC(C)[C@H]([NH3+])C(=O)[O-]` -> (S)
  - `LAC`: `C[C@H](O)C(=O)[O-]` -> (S)
- **Validation**: Evaluated chiral centers using RDKit `FindMolChiralCenters` to confirm full assignment and CIP configurations. Coarse-grained models in `gros/` and `itps/` are 2-bead or 3-bead achiral models invariant to enantiomeric reflection.

### 4. Fixed Guanine (GUA) Mapping
- **Problem**: `mappings/map_files/GUA.map` had a duplicate `GUA4` entry in the `[ martini ]` header line (`GUA1 GUA2 GUA3 GUA4 GUA4 GUA5 GUA6`), causing map parsers to fail.
- **Correction**: Removed the duplicate `GUA4`, leaving the 6 beads matching `GUA.itp` and `GUA.gro`:
  `GUA1 GUA2 GUA3 GUA4 GUA5 GUA6`.
- **Validation**: Verified against `itps/molecules/GUA.itp` and atom definitions in `GUA.map`.

### 5. Fixed UDP-GlcNAc (GUDP) Mapping & Topology
- **Problem**:
  - `mappings/map_files/GUDP.map` had duplicate `GYN2` and missing `GYN3` in the `[ martini ]` header line.
  - `itps/molecules/GUDP.itp` and `all_metabolites.itp` defined atom 16 with duplicate name `GYN4` and cgnr `15`.
  - `gros/GUDP.gro` labeled atom 16 as `GYN4`.
- **Correction**:
  - `GUDP.map`: Corrected `[ martini ]` line to `URA1 URA2 URA3 URA4 URA5 RIB1 RIB2 RIB3 RIB4 P1 P2 GYN1 GYN2 GYN3 GYN4 GYN5`.
  - `GUDP.itp` & `all_metabolites.itp`: Renamed atom 16 to `GYN5` with cgnr `16`.
  - `gros/GUDP.gro`: Renamed bead 16 to `GYN5`.
- **Validation**: Confirmed GlcNAc ring mapping and virtual site 16 construction on beads 12, 13, 14.

### 6. Fixed Urate (URT) Mapping and Bead Count
- **Problem**: `mappings/map_files/URT.map` had 7 beads in `[ martini ]` (`URT1`..`URT7`), with atom 3 (N1) mapped to `URT7`, while `URT.itp`, `all_metabolites.itp`, and `gros/URT.gro` define only 6 beads.
- **Correction**:
  - Mapped atom 3 (`N1`) to `URT6` (consistent with atom 13 `H1` being excluded from `URT6` via `!URT6`).
  - Removed extraneous `URT7` from `[ martini ]` line (`URT1 URT2 URT3 URT4 URT5 URT6`).
- **Validation**: Confirmed agreement with purine 6-bead architecture (identical to ADE and GUA) and topology virtual site definitions.

### 7. Added Automated Validation Suite
- Added `tests/validate_metabolome.py` and `tests/test_metabolome.py` to automate 10 scientific integrity checks across topologies, mappings, coordinates, and database records.

---

### Unresolved Scientific Ambiguities Documented for Maintainer Decision
1. **LHIS Charge State**: `itps/molecules/LHIS.itp` is parameterized with neutral beads `TN6d` and `TN5a` (net charge 0.0), matching typical physiological histidine (pKa ~ 6.0). In contrast, `database.csv`, Zenodo `AA/LHIS/LHIS.gro`, and `LHIS.map` represent the protonated HIP (+1) cation. Changing topology charges arbitrarily without force field reparameterization would violate Martini 3 bead definitions; maintainer confirmation is requested on whether to align database metadata to neutral or maintain dual documentation.
2. **PO4 Charge State**: `itps/molecules/PO4.itp` carries charge `0.0` for bead type `D` (annotated `; Note: other study`). In all other monophosphates in the metabolome, terminal `D` beads carry charge `-2.0` (matching divalent character). Because inorganic phosphate exists in physiological equilibrium between HPO4^2- and H2PO4^-, maintainer confirmation is needed to determine the intended protonation state before altering simulation charge parameters.
3. **Zenodo v1 All-Atom Starting Coordinates**: Initial gro files for LALA, LASP, LCYS, LGLU, LHIS, LILE, LSER, LVAL, and LAC in Zenodo v1 `AA.zip` were D-stereoisomers. Trajectories cannot be safely regenerated without original simulation setup files and workflows, so this is documented as an upstream data limitation.

