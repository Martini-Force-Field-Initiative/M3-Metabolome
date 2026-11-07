"""
Pytest integration for M3-Metabolome Validation Suite
======================================================
Tests all integrity constraints for the resolved GitHub Issue #4 metabolites,
and provides parameterized tests across the full metabolome database.
"""

import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

import pytest
from tests.validate_metabolome import (
    EXPECTED_STEREOCHEMISTRY,
    MetabolomeValidator,
)
TARGETS_ISSUE_4 = {"LASP", "LGLU", "LGLN", "LALA", "LSER", "LVAL", "LCYS", "LHIS", "LILE", "LAC", "CHT", "PO4", "GUA", "GUDP", "URT"}


@pytest.fixture(scope="module")
def validator():
    return MetabolomeValidator(REPO_ROOT)


def test_issue_4_duplicate_atoms_in_itp(validator):
    """Verify no duplicate atom names exist in ITP files for Issue #4 metabolites."""
    assert validator.check_duplicate_atoms_in_itp(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_duplicate_beads_in_map(validator):
    """Verify no duplicate bead names exist in MAP files for Issue #4 metabolites."""
    assert validator.check_duplicate_beads_in_map(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_missing_mapping_references(validator):
    """Verify all referenced atoms in MAP files map to valid beads for Issue #4 metabolites."""
    assert validator.check_missing_mapping_references(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_bead_count_consistency(validator):
    """Verify bead counts and names match across ITP, MAP, and GRO for Issue #4 metabolites."""
    assert validator.check_bead_count_consistency(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_formula_syntax(validator):
    """Verify molecular formulas have valid chemical syntax for Issue #4 metabolites."""
    assert validator.check_formula_syntax(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_smiles_parsing(validator):
    """Verify SMILES strings parse cleanly with RDKit for Issue #4 metabolites."""
    assert validator.check_smiles_parsing(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_formula_smiles_consistency(validator):
    """Verify database formula matches formula calculated from SMILES for Issue #4 metabolites."""
    assert validator.check_formula_smiles_consistency(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_stereochemistry(validator):
    """Verify stereochemical specification and CIP configuration for chiral Issue #4 metabolites."""
    assert validator.check_stereochemistry(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_topology_charges(validator):
    """Verify topology charge consistency against molecular definitions for Issue #4 metabolites."""
    assert validator.check_topology_charges(target_molecules=TARGETS_ISSUE_4)


def test_issue_4_cross_file_consistency(validator):
    """Verify cross-file consistency between database.csv, database.md, and all_metabolites.itp."""
    assert validator.check_cross_file_consistency(target_molecules=TARGETS_ISSUE_4)


@pytest.mark.parametrize("resname", sorted(EXPECTED_STEREOCHEMISTRY.keys()))
def test_individual_chiral_center(validator, resname):
    """Individually test each chiral metabolite to confirm correct stereochemical assignment."""
    assert validator.check_stereochemistry(target_molecules={resname})
