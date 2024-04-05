
import pandas as pd
from typing import Optional
from pathlib import Path
import numpy as np 
import logging, sys
from pydantic import BaseModel

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger()

class FitnessFactory(BaseModel):
    """
    Base class for handling Fitness data within `choppa`.
    """

    def read_fitness_csv(
            input_fitness_csv : Path,
            resindex_colname: Optional[str] = "residue_index",
            wildtype_colname: Optional[str] = "wildtype",
            mutant_colname: Optional[str] = "mutant",
            fitness_colname: Optional[str] = "fitness",
            confidence_colname: Optional[str] = None,
                        ):
        """
        Reads in a fitness CSV file and checks that all requested columns are present.
        """
        # read in the complete fitness data. This may have many columns
        logger.info(f"Reading in fitness data from {input_fitness_csv}")
        fitness_df = pd.read_csv(input_fitness_csv)

        # if no confidence is provided, just set to NaN
        if confidence_colname is None:
            fitness_df["confidence"] = np.nan
            confidence_colname = "confidence"
        elif confidence_colname not in fitness_df.columns:
            raise KeyError(f"Column {confidence_colname} not found in {input_fitness_csv}")

        # quick check for columns
        missing_columns = []
        for colname in [
            resindex_colname,
            wildtype_colname,
            mutant_colname,
            fitness_colname,
            ]:
            if not colname in fitness_df.columns:
                missing_columns.append(colname)
        if missing_columns:
            raise KeyError(f"Column(s) {missing_columns} not found in {input_fitness_csv}")
            
        # keep only the requested columns
        fitness_df = fitness_df[[
            resindex_colname,
            wildtype_colname,
            mutant_colname,
            fitness_colname,
            confidence_colname,
        ]]
        logger.info(f"Successfully read fitness data:\n{fitness_df}")

        return fitness_df
    
    def df_to_basedict(fitness_df):
        """
        Converts a `pandas` fitness dataframe (read by `FitnessFactory.read_fitness_csv`) into
        a fitness basedict.
        """




class ComplexFactory(BaseModel):
    """
    """
    def check_validity():
        """
        Does some quick checks to make sure the imported PDB structure is valid.
        """


if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE
    FitnessFactory.read_fitness_csv(TOY_FITNESS_DATA_COMPLETE, confidence_colname="confidence")