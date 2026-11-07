# Automated Validation Suite for Martini 3 Metabolome

This directory contains automated scientific validation tools for checking molecular topologies (`.itp`), mappings (`.map`), coarse-grained coordinates (`.gro`), and database entries (`misc/database.csv` and `database.md`).

## Requirements

The validation utility requires Python 3.9+ and the following packages:
- `rdkit`
- `pandas`
- `tabulate`
- `pytest` (optional, for pytest runner)

Install dependencies:
```bash
pip install rdkit pandas tabulate pytest
```
or via `uv`:
```bash
uv pip install rdkit pandas tabulate pytest
```

## Running the Validation

### Standalone CLI Utility
Run all 10 validation checks on the resolved Issue #4 metabolites:
```bash
python tests/validate_metabolome.py --targets LASP LGLU LGLN LALA LSER LVAL LCYS LHIS LILE LAC CHT PO4 GUA GUDP URT
```

Run validation across the entire repository:
```bash
python tests/validate_metabolome.py
```

Options:
- `--repo-root <path>`: Specify the root directory of the repository (default: `.`)
- `--targets <resname ...>`: Limit validation to specific metabolite resnames.

### Pytest Runner
Execute the test suite using `pytest`:
```bash
pytest tests/
```

## Validation Checks Performed

1. **Duplicate atom names in ITP files**: Ensures all bead/atom names within each `[ atoms ]` block are unique.
2. **Duplicate mapping entries in MAP files**: Ensures all beads declared under `[ martini ]` appear exactly once.
3. **Missing mapping references**: Ensures every atom in `[ atoms ]` references a bead defined under `[ martini ]`.
4. **ITP / MAP / GRO bead consistency**: Ensures bead counts and names match across `.itp`, `.map`, and `.gro` representations.
5. **Molecular formula validity**: Validates chemical formula syntax in `database.csv`.
6. **SMILES parsing via RDKit**: Ensures all canonical SMILES parse cleanly.
7. **Formula consistency with SMILES**: Confirms database formula matches the calculated molecular formula from SMILES.
8. **Stereochemical specification**: Validates that chiral centers are fully specified (no unassigned `?` centers) and verify absolute CIP configurations ((S) for L-amino acids / L-lactate, (R) for L-cysteine, (2S, 3S) for L-isoleucine).
9. **Topology charge consistency**: Verifies that total charge summed from ITP beads matches the expected formal charge from chemical definitions, flagging known scientific ambiguities.
10. **Cross-file consistency**: Ensures `database.md` is synchronized with `misc/database.csv`, and all individual topologies are present in `all_metabolites.itp`.

