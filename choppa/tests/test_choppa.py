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
from choppa.align import AlignFactory
from choppa import render


def test_choppa_render_toy_mac1_sectioned(tmp_path):
    """Tests that `choppa` is able to render views on a fitness dataset that
    has poor (fractioned) overlap with a toy PDB file."""

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_SECTIONED, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    # also test that we're actually reproducing the PDB sequence in the alignment.
    # prevents us from accidentally writing out the wrong logoplots.
    alignment = [
        resdict["wildtype"] for _, resdict in filled_aligned_fitness_dict.items()
    ]

    assert "".join([val.pop() for val in alignment[:4]]) == "SFSG"
    # assert "".join([val.pop() for val in alignment[6:25]]) == "KLTDNVYIKNADIVEEAKK"
    # TODO: fix this test
    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.pse"
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.html"
    ).render()

    # we know the intended output, test on this.
    assert len(filled_aligned_fitness_dict) == 164
    assert 4 not in filled_aligned_fitness_dict
    assert 5 in filled_aligned_fitness_dict

    # check files exist
    assert (tmp_path / "test.pse").exists()
    assert (tmp_path / "test.html").exists()


def test_choppa_correct_custom_columns():
    """Tests that `choppa` is able to read in correct columns in input fitness data."""

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED).get_fitness_basedict()
    # check that some of the correct keys are in the dict
    assert "fitness_csv_index" in fitness_dict[5]
    assert "wildtype" in fitness_dict[5]
    assert "mutants" in fitness_dict[5]

    # test that one of the residues has the intended number of mutants
    assert len(fitness_dict[5]["wildtype"]) == 3
    assert len(fitness_dict[5]["mutants"]) == 20


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


def test_choppa_render_toy_mac1_sectioned_noconf(tmp_path):
    """Tests that `choppa` is able to render views on a fitness dataset that
    has poor (fractioned) overlap with a toy PDB file, while not adding confidence."""

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.pse"
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.html"
    ).render()

    assert len(filled_aligned_fitness_dict) == 164

    # check files exist
    assert (tmp_path / "test.pse").exists()
    assert (tmp_path / "test.html").exists()


def test_choppa_render_toy_mac1_full(tmp_path):
    """Tests that `choppa` is able to render views on a fitness dataset that
    has complete overlap with a toy PDB file."""

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_COMPLETE, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.pse"
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.html"
    ).render()

    assert len(filled_aligned_fitness_dict) == 164

    # check files exist
    assert (tmp_path / "test.pse").exists()
    assert (tmp_path / "test.html").exists()

def test_choppa_render_toy_mac1_truncated(tmp_path):
    """Tests that `choppa` is able to render views on a fitness dataset that
    has decenr (trunacted at either end) overlap with a toy PDB file."""

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_TRUNCATED, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
    complex_rdkit = ComplexFactory(TOY_COMPLEX).load_pdb_rdkit()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    render.PublicationView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.pse"
    ).render()

    render.InteractiveView(
        filled_aligned_fitness_dict, complex, complex_rdkit, fitness_threshold=0.7, output_session_file=tmp_path/"test.html"
    ).render()

    assert len(filled_aligned_fitness_dict) == 164

    # check files exist
    assert (tmp_path / "test.pse").exists()
    assert (tmp_path / "test.html").exists()
