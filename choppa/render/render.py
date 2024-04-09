import logging, sys
import pymol2
from choppa.render.utils import show_contacts

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()

class PYMOL():
    """
    Uses the PyMOL API to create a session file for publication-ready views of the fitness data on 
    top of the complex PDB
    """

class HTML():
    """
    Uses 3DMol and Jinja to create a single HTML file that can be hosted anywhere to enable shareable 
    interactive views of the fitness data on top of the complex PDB
    """


if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_TRUNCATED
    from choppa.data.toy_data.resources import TOY_FITNESS_DATA_SECTIONED, TOY_FITNESS_DATA_COMPLETE_NOCONF

    from choppa.IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED, 
                                    # confidence_colname="confidence"
                                    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()

    from choppa.align.align import AlignFactory
    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()

    print(filled_aligned_fitness_dict)