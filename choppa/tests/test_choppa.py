"""
Unit and regression test for the choppa package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest


def test_choppa_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "choppa" in sys.modules


from choppa.data.toy_data.resources import (
    TOY_COMPLEX,
    TOY_FITNESS_DATA_COMPLETE,
    TOY_FITNESS_DATA_TRUNCATED,
)
from choppa.data.toy_data.resources import (
    TOY_FITNESS_DATA_SECTIONED,
    TOY_FITNESS_DATA_COMPLETE_NOCONF,
)

from choppa.IO.input import FitnessFactory, ComplexFactory
from choppa.align.align import AlignFactory
from choppa.render import render


def test_choppa_render_toy_mac1_sectioned():
    """Tests that `choppa` is able to render views on a fitness dataset that
    has poor (fractioned) overlap with a toy PDB file."""

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_SECTIONED, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()


def test_choppa_correct_custom_columns():
    """Tests that `choppa` is able to read in correct columns in input fitness data."""

    FitnessFactory(TOY_FITNESS_DATA_SECTIONED).get_fitness_basedict()


def test_choppa_incorrect_custom_columns():
    """Tests that `choppa` throws an error for incorrectly specified columns in input fitness data."""
    with pytest.raises(KeyError, match="not found in"):
        FitnessFactory(
            TOY_FITNESS_DATA_SECTIONED,
            wildtype_colname="foo",
            mutant_colname="foo",
            fitness_colname="foo",
            resindex_colname="foo",
        ).get_fitness_basedict()


def test_choppa_render_toy_mac1_sectioned_noconf():
    """Tests that `choppa` is able to render views on a fitness dataset that
    has poor (fractioned) overlap with a toy PDB file, while not adding confidence."""

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()


def test_choppa_render_toy_mac1_full():
    """Tests that `choppa` is able to render views on a fitness dataset that
    has complete overlap with a toy PDB file."""

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_COMPLETE, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()


def test_choppa_render_toy_mac1_truncated():
    """Tests that `choppa` is able to render views on a fitness dataset that
    has decenr (trunacted at either end) overlap with a toy PDB file."""

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_TRUNCATED, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7
    ).render()
