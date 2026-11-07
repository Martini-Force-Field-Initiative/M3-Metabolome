#!/usr/bin/env python3
"""
M3-Metabolome Validation Utility
=================================
Automated validation suite for the Martini 3 Metabolome dataset.

Validates:
 1. Duplicate atom names in ITP files.
 2. Duplicate mapping entries in MAP files.
 3. Missing mapping references in MAP files.
 4. ITP / MAP / GRO bead-count and name consistency.
 5. Molecular formula validity in database.csv.
 6. SMILES parsing using RDKit.
 7. Molecular formula consistency with SMILES.
 8. Stereochemical specification for chiral metabolites.
 9. Total topology charge consistency against molecular definitions.
 10. Cross-file consistency (database.csv, database.md, itp, map, gro, all_metabolites.itp).

Usage:
    python tests/validate_metabolome.py [--repo-root <path>] [--target <resname> ...]
    pytest tests/test_metabolome.py
"""

import argparse
import glob
import os
import re
import sys
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


# Authoritative L-amino acid and L-hydroxy acid stereochemical definitions:
# CIP stereocenters for standard L-enantiomers
EXPECTED_STEREOCHEMISTRY = {
    "LALA": [(1, "S")],             # L-Alanine
    "LASP": [(4, "S")],             # L-Aspartate
    "LCYS": [(2, "R")],             # L-Cysteine (priority of -CH2SH > -COOH reverses CIP descriptor to R)
    "LGLN": [(5, "S")],             # L-Glutamine
    "LGLU": [(5, "S")],             # L-Glutamate
    "LSER": [(2, "S")],             # L-Serine
    "LVAL": [(3, "S")],             # L-Valine
    "LILE": [(2, "S"), (4, "S")],   # L-Isoleucine (2S, 3S)
    "LHIS": [(3, "S")],             # L-Histidine
    "LAC":  [(1, "S")],             # L-Lactate ((S)-2-hydroxypropanoate)
}

# Known scientific ambiguities flagged for maintainer decision:
KNOWN_SCIENTIFIC_AMBIGUITIES = {
    "LHIS": "LHIS.itp is parameterized as neutral histidine (zwitterion, net charge 0.0 with TN6d/TN5a beads), whereas database metadata describes protonated histidine (C6H10N3O2+, +1).",
    "PO4":  "PO4.itp has charge 0.0 (marked 'Note: other study'), whereas physiological phosphate exists as an anion (predominantly HPO4^2- / H2PO4-).",
}


def normalize_formula(formula: str) -> str:
    """Normalize chemical formula string to allow comparison between
    repeated signs (e.g. '--') and numeric charges (e.g. '-2')."""
    if not isinstance(formula, str):
        return ""
    f = formula.strip()
    # Replace trailing repeated + or - with numeric charge:
    # e.g. '----' -> '-4', '---' -> '-3', '--' -> '-2', '++++' -> '+4'
    match_minus = re.search(r"(-+)$", f)
    if match_minus and len(match_minus.group(1)) > 1:
        count = len(match_minus.group(1))
        f = f[:match_minus.start()] + f"-{count}"
    match_plus = re.search(r"(\++)$", f)
    if match_plus and len(match_plus.group(1)) > 1:
        count = len(match_plus.group(1))
        f = f[:match_plus.start()] + f"+{count}"
    return f


class MetabolomeValidator:
    def __init__(self, repo_root: str):
        self.repo_root = os.path.abspath(repo_root)
        self.db_csv_path = os.path.join(self.repo_root, "misc", "database.csv")
        self.db_md_path = os.path.join(self.repo_root, "database.md")
        self.itp_dir = os.path.join(self.repo_root, "itps", "molecules")
        self.map_dir = os.path.join(self.repo_root, "mappings", "map_files")
        self.gro_dir = os.path.join(self.repo_root, "gros")
        self.all_itp_path = os.path.join(self.itp_dir, "all_metabolites.itp")

        self.errors: List[str] = []
        self.warnings: List[str] = []

    def log_error(self, check_name: str, file_path: str, message: str, line_no: Optional[int] = None):
        rel_path = os.path.relpath(file_path, self.repo_root)
        loc = f"{rel_path}:{line_no}" if line_no else rel_path
        self.errors.append(f"[{check_name}] {loc} - {message}")

    def log_warning(self, check_name: str, file_path: str, message: str, line_no: Optional[int] = None):
        rel_path = os.path.relpath(file_path, self.repo_root)
        loc = f"{rel_path}:{line_no}" if line_no else rel_path
        self.warnings.append(f"[{check_name}] {loc} - {message}")

    # 1. Duplicate atom names in ITP files
    def check_duplicate_atoms_in_itp(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify no ITP file has duplicate bead/atom names within a [ moleculetype ]."""
        passed = True
        itp_files = sorted(glob.glob(os.path.join(self.itp_dir, "*.itp")))
        for itp_file in itp_files:
            if os.path.basename(itp_file) == "all_metabolites.itp":
                continue
            mol_name = os.path.basename(itp_file)[:-4]
            if target_molecules and mol_name not in target_molecules:
                continue

            with open(itp_file, "r", encoding="utf-8") as f:
                in_atoms = False
                seen_atoms: Dict[str, int] = {}
                for line_no, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith(";"):
                        continue
                    if line.startswith("["):
                        in_atoms = bool(re.search(r"\[\s*atoms\s*\]", line))
                        continue
                    if in_atoms:
                        parts = line.split()
                        if len(parts) >= 5:
                            atom_name = parts[4]
                            if atom_name in seen_atoms:
                                passed = False
                                self.log_error(
                                    "ITP_DUPLICATE_ATOM",
                                    itp_file,
                                    f"Offending atom '{atom_name}': duplicate atom name (first defined on line {seen_atoms[atom_name]})",
                                    line_no=line_no,
                                )
                            else:
                                seen_atoms[atom_name] = line_no
        return passed

    # 2. Duplicate mapping entries in MAP files
    def check_duplicate_beads_in_map(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify no MAP file lists duplicate beads in its [ martini ] header section."""
        passed = True
        map_files = sorted(glob.glob(os.path.join(self.map_dir, "*.map")))
        for map_file in map_files:
            mol_name = os.path.basename(map_file)[:-4]
            if target_molecules and mol_name not in target_molecules:
                continue

            with open(map_file, "r", encoding="utf-8") as f:
                in_martini = False
                seen_beads: Dict[str, int] = {}
                for line_no, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith(";"):
                        continue
                    if line.startswith("["):
                        in_martini = bool(re.search(r"\[\s*martini\s*\]", line))
                        continue
                    if in_martini:
                        beads = line.split()
                        for bead in beads:
                            if bead in seen_beads:
                                passed = False
                                self.log_error(
                                    "MAP_DUPLICATE_BEAD",
                                    map_file,
                                    f"Offending bead '{bead}': duplicate bead in [ martini ] header (first seen on line {seen_beads[bead]})",
                                    line_no=line_no,
                                )
                            else:
                                seen_beads[bead] = line_no
        return passed

    # 3. Missing mapping references in MAP files
    def check_missing_mapping_references(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify every bead referenced in the [ atoms ] section of a MAP file exists in [ martini ]."""
        passed = True
        map_files = sorted(glob.glob(os.path.join(self.map_dir, "*.map")))
        for map_file in map_files:
            mol_name = os.path.basename(map_file)[:-4]
            if target_molecules and mol_name not in target_molecules:
                continue

            with open(map_file, "r", encoding="utf-8") as f:
                in_martini = False
                in_atoms = False
                martini_beads: Set[str] = set()
                atom_references: List[Tuple[int, str, str]] = []  # line_no, atom_name, bead_ref

                for line_no, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith(";"):
                        continue
                    if line.startswith("["):
                        in_martini = bool(re.search(r"\[\s*martini\s*\]", line))
                        in_atoms = bool(re.search(r"\[\s*atoms\s*\]", line))
                        continue
                    if in_martini:
                        martini_beads.update(line.split())
                    elif in_atoms:
                        parts = line.split()
                        if len(parts) >= 3:
                            atom_name = parts[1]
                            for bead_ref in parts[2:]:
                                clean_ref = bead_ref.lstrip("!")
                                atom_references.append((line_no, atom_name, clean_ref))

            for line_no, atom_name, bead_ref in atom_references:
                if bead_ref not in martini_beads:
                    passed = False
                    self.log_error(
                        "MAP_MISSING_REFERENCE",
                        map_file,
                        f"Offending atom '{atom_name}': references bead '{bead_ref}' which is missing from [ martini ]",
                        line_no=line_no,
                    )
        return passed

    # 4. ITP / MAP / GRO bead-count and name consistency
    def check_bead_count_consistency(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify that bead counts and names match across .itp, .map, and .gro files."""
        passed = True
        itp_files = sorted(glob.glob(os.path.join(self.itp_dir, "*.itp")))
        for itp_file in itp_files:
            if os.path.basename(itp_file) == "all_metabolites.itp":
                continue
            mol_name = os.path.basename(itp_file)[:-4]
            if target_molecules and mol_name not in target_molecules:
                continue

            map_file = os.path.join(self.map_dir, f"{mol_name}.map")
            gro_file = os.path.join(self.gro_dir, f"{mol_name}.gro")

            # Parse ITP beads
            itp_beads: List[str] = []
            with open(itp_file, "r", encoding="utf-8") as f:
                in_atoms = False
                for line in f:
                    line = line.strip()
                    if not line or line.startswith(";"):
                        continue
                    if line.startswith("["):
                        in_atoms = bool(re.search(r"\[\s*atoms\s*\]", line))
                        continue
                    if in_atoms:
                        parts = line.split()
                        if len(parts) >= 5:
                            itp_beads.append(parts[4])

            # Parse MAP beads
            if os.path.exists(map_file):
                map_beads: List[str] = []
                with open(map_file, "r", encoding="utf-8") as f:
                    in_martini = False
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith(";"):
                            continue
                        if line.startswith("["):
                            in_martini = bool(re.search(r"\[\s*martini\s*\]", line))
                            continue
                        if in_martini:
                            map_beads.extend(line.split())

                if len(itp_beads) != len(map_beads):
                    passed = False
                    self.log_error(
                        "BEAD_COUNT_MISMATCH",
                        map_file,
                        f"Bead count mismatch: ITP has {len(itp_beads)} beads {itp_beads}, but MAP has {len(map_beads)} beads {map_beads}",
                    )

            # Parse GRO beads
            if os.path.exists(gro_file):
                gro_beads: List[str] = []
                with open(gro_file, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                    if len(lines) >= 3:
                        try:
                            num_atoms = int(lines[1].strip())
                            for line in lines[2 : 2 + num_atoms]:
                                atom_name = line[10:15].strip()
                                gro_beads.append(atom_name)
                        except ValueError:
                            pass
                if len(itp_beads) != len(gro_beads):
                    passed = False
                    self.log_error(
                        "GRO_BEAD_COUNT_MISMATCH",
                        gro_file,
                        f"Bead count mismatch: ITP has {len(itp_beads)} beads {itp_beads}, but GRO has {len(gro_beads)} beads {gro_beads}",
                    )
        return passed

    # 5. Invalid molecular formulas in database.csv
    def check_formula_syntax(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Validate syntax of molecular formulas in database.csv."""
        passed = True
        df = pd.read_csv(self.db_csv_path)
        formula_regex = re.compile(r"^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*([+-]+|\d*[+-])?$")
        for idx, row in df.iterrows():
            res = str(row["resname"]).strip()
            if target_molecules and res not in target_molecules:
                continue
            formula = row["Formula"]
            if pd.isna(formula) or not str(formula).strip():
                continue
            f_str = str(formula).strip()
            if not formula_regex.match(f_str):
                passed = False
                self.log_error(
                    "INVALID_FORMULA_SYNTAX",
                    self.db_csv_path,
                    f"Offending metabolite '{res}': malformed formula string '{f_str}'",
                    line_no=idx + 2,
                )
        return passed

    # 6. SMILES parsing using RDKit
    def check_smiles_parsing(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify that every Canonical SMILES in database.csv parses successfully with RDKit."""
        passed = True
        df = pd.read_csv(self.db_csv_path)
        for idx, row in df.iterrows():
            res = str(row["resname"]).strip()
            if target_molecules and res not in target_molecules:
                continue
            smiles = row["Canonical SMILES"]
            if pd.isna(smiles) or not str(smiles).strip():
                continue
            s_str = str(smiles).strip()
            mol = Chem.MolFromSmiles(s_str)
            if mol is None:
                passed = False
                self.log_error(
                    "SMILES_PARSE_FAILURE",
                    self.db_csv_path,
                    f"Offending metabolite '{res}': RDKit failed to parse SMILES '{s_str}'",
                    line_no=idx + 2,
                )
        return passed

    # 7. Molecular formula consistency with SMILES
    def check_formula_smiles_consistency(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify that database formula matches the formula computed from SMILES."""
        passed = True
        df = pd.read_csv(self.db_csv_path)
        for idx, row in df.iterrows():
            res = str(row["resname"]).strip()
            if target_molecules and res not in target_molecules:
                continue
            smiles = row["Canonical SMILES"]
            formula = row["Formula"]
            if pd.isna(smiles) or pd.isna(formula):
                continue
            s_str = str(smiles).strip()
            f_str = str(formula).strip()
            mol = Chem.MolFromSmiles(s_str)
            if mol is None:
                continue
            calc_formula = CalcMolFormula(mol)
            if normalize_formula(f_str) != normalize_formula(calc_formula):
                passed = False
                self.log_error(
                    "FORMULA_SMILES_MISMATCH",
                    self.db_csv_path,
                    f"Offending metabolite '{res}': database formula '{f_str}' != calculated SMILES formula '{calc_formula}'",
                    line_no=idx + 2,
                )
        return passed

    # 8. Stereochemical specification for relevant chiral molecules
    def check_stereochemistry(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Validate chiral center specification and absolute configuration for relevant chiral molecules."""
        passed = True
        df = pd.read_csv(self.db_csv_path)
        for res, expected in EXPECTED_STEREOCHEMISTRY.items():
            if target_molecules and res not in target_molecules:
                continue
            sub = df[df["resname"] == res]
            if sub.empty:
                continue
            idx = sub.index[0]
            smiles = sub.iloc[0]["Canonical SMILES"]
            if pd.isna(smiles) or not str(smiles).strip():
                passed = False
                self.log_error("CHIRALITY_MISSING", self.db_csv_path, f"Chiral molecule '{res}' has no SMILES", line_no=idx + 2)
                continue
            mol = Chem.MolFromSmiles(str(smiles).strip())
            if mol is None:
                continue

            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

            # Check unassigned centers
            unassigned = [c for c in centers if c[1] not in ("R", "S")]
            if unassigned:
                passed = False
                self.log_error(
                    "CHIRALITY_UNASSIGNED",
                    self.db_csv_path,
                    f"Offending metabolite '{res}': SMILES lacks specified stereocenter(s) {unassigned}",
                    line_no=idx + 2,
                )
                continue

            # Compare against expected CIP configuration
            actual_configs = [c[1] for c in centers]
            expected_configs = [c[1] for c in expected]
            if actual_configs != expected_configs:
                passed = False
                self.log_error(
                    "CHIRALITY_INVERSION",
                    self.db_csv_path,
                    f"Offending metabolite '{res}': stereochemistry {actual_configs} does not match expected {expected_configs}",
                    line_no=idx + 2,
                )
        return passed

    # 9. Total topology charge consistency
    def check_topology_charges(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Validate that total charge calculated from ITP matches expected formal charge."""
        passed = True
        df = pd.read_csv(self.db_csv_path)
        for idx, row in df.iterrows():
            res = str(row["resname"]).strip()
            if target_molecules and res not in target_molecules:
                continue
            itp_file = os.path.join(self.itp_dir, f"{res}.itp")
            if not os.path.exists(itp_file):
                continue

            # Calculate total charge from ITP
            itp_charge = 0.0
            with open(itp_file, "r", encoding="utf-8") as f:
                in_atoms = False
                for line in f:
                    line = line.strip()
                    if not line or line.startswith(";"):
                        continue
                    if line.startswith("["):
                        in_atoms = bool(re.search(r"\[\s*atoms\s*\]", line))
                        continue
                    if in_atoms:
                        p = line.split()
                        if len(p) >= 7:
                            itp_charge += float(p[6])

            # Calculate expected charge from SMILES if available
            smiles = row["Canonical SMILES"]
            if pd.notna(smiles) and str(smiles).strip():
                mol = Chem.MolFromSmiles(str(smiles).strip())
                if mol is not None:
                    expected_charge = float(Chem.GetFormalCharge(mol))
                    if abs(itp_charge - expected_charge) > 1e-3:
                        if res in KNOWN_SCIENTIFIC_AMBIGUITIES:
                            self.log_warning(
                                "KNOWN_CHARGE_AMBIGUITY",
                                itp_file,
                                f"Known scientific ambiguity for '{res}': topology charge ({itp_charge:.1f}) != SMILES charge ({expected_charge:.1f}). Context: {KNOWN_SCIENTIFIC_AMBIGUITIES[res]}",
                            )
                        else:
                            passed = False
                            self.log_error(
                                "TOPOLOGY_CHARGE_MISMATCH",
                                itp_file,
                                f"Offending metabolite '{res}': topology charge ({itp_charge:.1f}) != database SMILES formal charge ({expected_charge:.1f})",
                            )
            elif res in KNOWN_SCIENTIFIC_AMBIGUITIES:
                self.log_warning(
                    "KNOWN_CHARGE_AMBIGUITY",
                    itp_file,
                    f"Known scientific ambiguity for '{res}': topology charge is {itp_charge:.1f}. Context: {KNOWN_SCIENTIFIC_AMBIGUITIES[res]}",
                )
        return passed

    # 10. Cross-file molecular identity consistency
    def check_cross_file_consistency(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Verify cross-file consistency between database.csv, database.md, and all_metabolites.itp."""
        passed = True

        # Check all_metabolites.itp synchronization with individual ITPs
        if os.path.exists(self.all_itp_path):
            with open(self.all_itp_path, "r", encoding="utf-8") as f:
                all_itp_content = f.read()

            itp_files = sorted(glob.glob(os.path.join(self.itp_dir, "*.itp")))
            for itp_file in itp_files:
                if os.path.basename(itp_file) == "all_metabolites.itp":
                    continue
                res = os.path.basename(itp_file)[:-4]
                if target_molecules and res not in target_molecules:
                    continue

                # Check that molecule is present in all_metabolites.itp
                if f"[ moleculetype ]\n{res} " not in all_itp_content and f"[ moleculetype ]\n{res}\t" not in all_itp_content:
                    passed = False
                    self.log_error("CROSS_FILE_MISSING", self.all_itp_path, f"Molecule '{res}' missing from all_metabolites.itp")

        # Check database.md synchronization with database.csv
        if os.path.exists(self.db_csv_path) and os.path.exists(self.db_md_path):
            df_csv = pd.read_csv(self.db_csv_path)
            for col in ["Metabolite name", "class"]:
                df_csv[col] = df_csv[col].str.replace("_", " ", regex=False)
            expected_md = df_csv.astype(object).fillna("").to_markdown(index=False, tablefmt="github")
            with open(self.db_md_path, "r", encoding="utf-8") as f:
                current_md = f.read()
            if current_md.strip() != expected_md.strip():
                passed = False
                self.log_error(
                    "DATABASE_MD_OUT_OF_SYNC",
                    self.db_md_path,
                    "database.md is out of sync with misc/database.csv. Run 'python misc/database_md.py' to regenerate.",
                )

        return passed

    def run_all(self, target_molecules: Optional[Set[str]] = None) -> bool:
        """Run all 10 validation checks and return overall pass/fail boolean."""
        self.errors.clear()
        self.warnings.clear()

        checks = [
            ("1. Duplicate atom names in ITP", self.check_duplicate_atoms_in_itp),
            ("2. Duplicate mapping entries in MAP", self.check_duplicate_beads_in_map),
            ("3. Missing mapping references in MAP", self.check_missing_mapping_references),
            ("4. ITP/MAP/GRO bead count consistency", self.check_bead_count_consistency),
            ("5. Molecular formula validity", self.check_formula_syntax),
            ("6. SMILES parsing via RDKit", self.check_smiles_parsing),
            ("7. Formula consistency with SMILES", self.check_formula_smiles_consistency),
            ("8. Stereochemical specification", self.check_stereochemistry),
            ("9. Topology charge consistency", self.check_topology_charges),
            ("10. Cross-file consistency", self.check_cross_file_consistency),
        ]

        print("=" * 75)
        print("M3-METABOLOME VALIDATION SUITE")
        print("=" * 75)

        all_passed = True
        for name, func in checks:
            check_passed = func(target_molecules)
            status = "PASS" if check_passed else "FAIL"
            print(f"[{status}] {name}")
            if not check_passed:
                all_passed = False

        print("\n" + "=" * 75)
        print("VALIDATION SUMMARY")
        print("=" * 75)

        if self.warnings:
            print(f"\nWARNINGS ({len(self.warnings)}):")
            for w in self.warnings:
                print(f"  {w}")

        if self.errors:
            print(f"\nERRORS ({len(self.errors)}):")
            for e in self.errors:
                print(f"  {e}")
            print(f"\nStatus: FAILED with {len(self.errors)} error(s).")
            return False
        else:
            print(f"\nStatus: PASSED (0 errors, {len(self.warnings)} documented warning(s)).")
            return True


def main():
    parser = argparse.ArgumentParser(description="Automated validation utility for M3-Metabolome")
    parser.add_argument("--repo-root", default=".", help="Root path of the M3-Metabolome repository")
    parser.add_argument("--targets", nargs="*", default=None, help="Specific metabolite resnames to validate (e.g. LASP CHT GUDP)")
    args = parser.parse_args()

    target_set = set(args.targets) if args.targets else None
    validator = MetabolomeValidator(args.repo_root)
    success = validator.run_all(target_molecules=target_set)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
