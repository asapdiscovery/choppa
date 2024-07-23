import pandas as pd
from typing import Optional
from pathlib import Path
import numpy as np
import logging, sys
from collections import OrderedDict

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()


class FitnessFactory:
    """
    Base class for handling Fitness data within `choppa`.
    """

    def __init__(
        self,
        input_fitness_csv: Path,
        resindex_colname: Optional[str] = "residue_index",
        wildtype_colname: Optional[str] = "wildtype",
        mutant_colname: Optional[str] = "mutant",
        fitness_colname: Optional[str] = "fitness",
        confidence_colname: Optional[str] = None,
    ):
        self.input_fitness_csv = input_fitness_csv
        self.resindex_colname = resindex_colname
        self.wildtype_colname = wildtype_colname
        self.mutant_colname = mutant_colname
        self.fitness_colname = fitness_colname
        self.confidence_colname = confidence_colname
        if self.confidence_colname is not None:
            # makes it easier to track whether confidence values are provided or set to NaN by us
            self.confidence_set = True
        else:
            self.confidence_set = False
        self.fitness_df = None

    def check_validity(self, fitness_df):
        """
        Does some quick checks to make sure the imported CSV file is valid.
        """
        if self.confidence_colname not in fitness_df.columns:
            raise KeyError(
                f"Column {self.confidence_colname} not found in {self.input_fitness_csv}"
            )

        # quick check for columns
        non_confidence_columns = [
            self.resindex_colname,
            self.wildtype_colname,
            self.mutant_colname,
            self.fitness_colname,
        ]
        missing_columns = []
        for colname in non_confidence_columns:
            if not colname in fitness_df.columns:
                missing_columns.append(colname)
        if missing_columns:
            raise KeyError(
                f"Column(s) {missing_columns} not found in {self.input_fitness_csv}"
            )

        # keep only the requested columns
        fitness_df = fitness_df[non_confidence_columns + [self.confidence_colname]]

        # if residue_index is a float (ran into this use-case), convert to int.
        fitness_df[self.resindex_colname] = fitness_df[self.resindex_colname].apply(
            np.int64
        )

        # check that there aren't any NaNs and that fitness (and confidence) data is scalar
        if fitness_df[non_confidence_columns].isnull().values.any():
            raise ValueError(
                f"Found missing values in input CSV: {fitness_df[fitness_df[non_confidence_columns].isnull().any(axis=1)]}"
            )
        if (
            len(
                fitness_df[
                    pd.to_numeric(
                        fitness_df[self.fitness_colname], errors="coerce"
                    ).isnull()
                ]
            )
            > 0
        ):
            raise ValueError(
                f"Found non-numeric fitness values in input CSV: {fitness_df[pd.to_numeric(fitness_df[self.fitness_colname], errors='coerce').isnull()]}"
            )
        if self.confidence_colname is not None:
            if self.confidence_set:
                if (
                    len(
                        fitness_df[
                            pd.to_numeric(
                                fitness_df[self.confidence_colname], errors="coerce"
                            ).isnull()
                        ]
                    )
                    > 0
                ):
                    raise ValueError(
                        f"Found non-numeric confidence values in input CSV: {fitness_df[pd.to_numeric(fitness_df[self.confidence_colname], errors='coerce').isnull()]}"
                    )

        # if checks reach this point then the input data should be correctly formatted. Rename and adopt.
        fitness_df = fitness_df.rename(
            columns={
                self.resindex_colname: "residue_index",
                self.wildtype_colname: "wildtype",
                self.mutant_colname: "mutant",
                self.fitness_colname: "fitness",
            }
        )

        for resi, res_data in fitness_df.groupby(by="residue_index"):
            if len(res_data) > 25:
                raise ValueError(
                    f"Found residue indices in input fitness CSV (at index {resi}) with more mutants than expected ({len(res_data)})! Does your fitness data have "
                    f"multiple chains in it with overlapping residue indices? Please resolve the input data so that there is no overlap in "
                    f"residue indices between chains."
                )

        self.fitness_df = fitness_df
        return True

    def read_fitness_csv(self):
        """
        Reads in a fitness CSV file and checks that all requested columns are present, complete and numeric.
        """
        # read in the complete fitness data. This may have many columns
        logger.info(f"Reading in fitness data from {self.input_fitness_csv}")
        fitness_df = pd.read_csv(self.input_fitness_csv)

        # if no confidence is provided, just set to NaN
        if self.confidence_colname is None:
            fitness_df["confidence"] = np.nan
            self.confidence_colname = "confidence"

        # check whether the CSV file is correct
        if self.check_validity(fitness_df):
            logger.info(f"Successfully read fitness data:\n{self.fitness_df}")
            return self.fitness_df

    def get_fitness_basedict(self):
        """
        Converts a `pandas` fitness dataframe (read by `FitnessFactory.read_fitness_csv`) into
        a `fitness basedict` which is essentially just an `OrderedDict`.

        We want the dict to have the form:
        {
        residue_index: {
                fitness_csv_index, # this is the original index in the fitness CSV for provenance
                wildtype: {AA, fitness, confidence},
                mutants: [{AA, fitness, confidence}, {AA, fitness, confidence}, etc]
            }
        }

        """
        fitness_basedict = OrderedDict()

        for residx, res_df in self.read_fitness_csv().groupby(by="residue_index"):

            # first construct the dict entry for the wildtype
            wildtype = res_df["wildtype"].values[0]

            # if there is no wildtype mutation available, the experimentalists have omitted it and set
            # the value implicitly as 0.
            if len(res_df[res_df.mutant == wildtype]) == 0:
                res_df.loc[-1] = [
                    res_df["residue_index"].values[0],
                    res_df["wildtype"].values[0],
                    res_df["wildtype"].values[0],  # this is where we add the wildtype
                    0.0,
                    res_df["confidence"].values[0],
                ]

            wt_dict = {
                "aa": wildtype,
                "fitness": res_df[res_df.mutant == wildtype]["fitness"].values[0],
                "confidence": res_df[res_df.mutant == wildtype]["confidence"].values[0],
            }
            # now construct the mutant list of dict entries for this residue index while excluding wildtype
            mutant_list = []
            for mut, fitness, conf in res_df[
                [
                    "mutant",
                    "fitness",
                    "confidence",
                ]
            ].values:
                if not mut == wildtype:  # ignore wildtype entry
                    mutant_list.append(
                        {"aa": mut, "fitness": fitness, "confidence": conf}
                    )

            fitness_basedict[residx] = {
                "fitness_csv_index": residx,  # double now, but will be helpful for provenance after alignment
                "wildtype": wt_dict,
                "mutants": mutant_list,
            }

        logger.info(
            f"Created fitness dictionary as `FitnessFactory` of length {len(fitness_basedict)}"
        )
        return fitness_basedict


from Bio.PDB import PDBParser
from rdkit import Chem


class ComplexFactory:
    """ """

    def __init__(
        self,
        path_to_pdb_file: Path,
    ):
        self.path_to_pdb_file = path_to_pdb_file

    def remove_waters(system):
        """[Placeholder] Returns a system with water entries removed"""
        return system

    def extract_ligands(system):
        """[Placeholder] Returns a system's ligands"""
        return system

    def reset_complex_sequence(self, complex):
        """Adjusts protein sequence indexing to run from 1 to n, rather than whatever wonky indexing
        the crystallographer may have come up with. BioPython can do this but there are some
        protections built in against hard re-indexing multi-chain proteins. By first setting
        the indexing super high (starting at 100,000) and then re-indexing starting from 1 we can
        circumvent these protections. For more details see https://github.com/biopython/biopython/pull/4623
        """
        # first set indexing to an unphysically high number
        original_index = []
        residue_N = 100000
        for residue in complex.get_residues():
            original_index.append(residue.id[1])
            residue.id = (residue.id[0], residue_N, residue.id[2])
            residue_N += 1

        # now renumber residue in complex starting from 1
        residue_N = 1
        for residue in complex.get_residues():
            residue.id = (residue.id[0], residue_N, residue.id[2])
            residue_N += 1
        return complex, original_index

    def load_pdb(self):
        """
        Loads an input PDB file
        """
        complex = PDBParser(QUIET=False).get_structure("COMPLEX", self.path_to_pdb_file)

        self.reset_complex_sequence(complex)
        return complex

    def load_pdb_rdkit(self):
        """
        Loads an input PDB file to an RDKit object for easier string retrieval
        """
        return Chem.MolFromPDBFile(self.path_to_pdb_file, sanitize=False)


if __name__ == "__main__":
    from choppa.data.toy_data.resources import (
        TOY_COMPLEX,
        TOY_FITNESS_DATA_COMPLETE,
        TOY_FITNESS_DATA_COMPLETE_NOCONF,
    )

    fitness_df = FitnessFactory(
        TOY_FITNESS_DATA_COMPLETE, confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()
