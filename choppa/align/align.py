
import pandas as pd
from typing import Optional
from pathlib import Path
import numpy as np 
import logging, sys
from collections import OrderedDict
from Bio import PDB
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

class AlignFactory():
    """
    Base class for aligning Fitness data with PDB complexes within `choppa`.
    """
    def __init__(self,
        fitness_dict: OrderedDict,
        complex: PDB.Structure.Structure,
        ):
        self.fitness_input = fitness_dict
        self.complex = complex
   
    def validate_alignment(self):
        """
        [Placeholder] validates the alignment of a fitness OrderedDict to a complex object
        """
    def align_fitness(self):
        """
        [Placeholder] align a fitness OrderedDict to a complex object
        """


if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_COMPLETE_NOCONF
    from choppa.IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_COMPLETE, 
                                    confidence_colname="confidence"
                                    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()

    AlignFactory(fitness_dict, complex)
    # print(type(fitness_dict), type(complex))
    # AlignFactory(fitness_dict, complex)

    
    