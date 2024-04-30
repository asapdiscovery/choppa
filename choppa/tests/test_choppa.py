"""
Unit and regression test for the choppa package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import choppa


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


def test_choppa_render():
    """Sample test, will always pass so long as import statement worked."""

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
