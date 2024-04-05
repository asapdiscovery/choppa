
import pandas as pd
from typing import Optional
from pathlib import Path
import numpy as np 
import logging, sys
from collections import OrderedDict

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

class FitnessFactory():
    """
    Base class for handling Fitness data within `choppa`.
    """
    def __init__(self,
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
   
    def check_validity(self, fitness_df):
        """
        Does some quick checks to make sure the imported CSV file is valid.
        """
        if self.confidence_colname not in fitness_df.columns:
            raise KeyError(f"Column {self.confidence_colname} not found in {self.input_fitness_csv}")

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
            raise KeyError(f"Column(s) {missing_columns} not found in {self.input_fitness_csv}")
            
        # keep only the requested columns
        fitness_df = fitness_df[non_confidence_columns + [self.confidence_colname]]

        # check that there aren't any NaNs and that fitness (and confidence) data is scalar
        if fitness_df[non_confidence_columns].isnull().values.any():
            raise ValueError(f"Found missing values in input CSV: {fitness_df[fitness_df[non_confidence_columns].isnull().any(axis=1)]}")
        if len(fitness_df[pd.to_numeric(fitness_df[self.fitness_colname], errors='coerce').isnull()]) > 0:
            raise ValueError(f"Found non-numeric fitness values in input CSV: {fitness_df[pd.to_numeric(fitness_df[self.fitness_colname], errors='coerce').isnull()]}")
        if self.confidence_colname is not None:
            if self.confidence_set:
                if len(fitness_df[pd.to_numeric(fitness_df[self.confidence_colname], errors='coerce').isnull()]) > 0:
                    raise ValueError(f"Found non-numeric confidence values in input CSV: {fitness_df[pd.to_numeric(fitness_df[self.confidence_colname], errors='coerce').isnull()]}")
            
        # if checks reach this point then the input data should be correctly formatted
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
            logger.info(f"Successfully read fitness data:\n{fitness_df}")
            return fitness_df
        
    
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
        for residx, res_df in self.read_fitness_csv().groupby(by=self.resindex_colname):
            # first construct the dict entry for the wildtype
            wildtype = res_df[self.wildtype_colname].values[0]
            wt_dict = {
                        "aa" : wildtype, 
                        "fitness" : res_df[res_df.mutant == wildtype]["fitness"].values[0],
                        "confidence" : res_df[res_df.mutant == wildtype]["confidence"].values[0]
                             }
            # now construct the mutant list of dict entries for this residue index while excluding wildtype
            mutant_list = []
            for mut, fitness, conf in res_df[[
                    self.mutant_colname,
                    self.fitness_colname,
                    self.confidence_colname,
                    ]].values:
                if not mut == wildtype: # ignore wildtype entry
                    mutant_list.append({
                        "aa" : mut, 
                        "fitness" : fitness,
                        "confidence" : conf})
                    
            fitness_basedict[residx] = {
                "fitness_csv_index" : residx, # double now, but will be helpful for provenance after alignment
                "wildtype": wt_dict,
                "mutants": mutant_list
            }
        
        logger.info(f"Created fitness dictionary as `FitnessFactory` of length {len(fitness_basedict)}")

        return fitness_basedict


class ComplexFactory():
    """
    """
    def check_validity():
        """
        Does some quick checks to make sure the imported PDB structure is valid.
        """


if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_COMPLETE_NOCONF
    fitness_df = FitnessFactory(TOY_FITNESS_DATA_COMPLETE_NOCONF, 
                                    # confidence_colname="confidence"
                                    ).get_fitness_basedict()
    