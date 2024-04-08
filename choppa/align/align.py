
import logging, sys
from collections import OrderedDict
from Bio import PDB, SeqUtils, Align
from Bio.Align import substitution_matrices

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
    def fitness_get_seq(self):
        """
        From a fitness `OrderedDict`, extracts the amino acid sequence
        """
        seq_list = [ fitness_values['wildtype']['aa'] for idx, fitness_values in self.fitness_input.items() ]
        
        return "".join(seq_list)
    
    def complex_get_seq(self):
        """
        From a `PDB.Structure.Structure`, extracts the amino acid sequence
        """
        seq_list = [ SeqUtils.seq1(res.get_resname()) for res in self.complex.get_residues() ]

        return "".join(seq_list)        


    def fitness_reset_keys(self, alignment_shift_dict):
        """
        Given a fitness `OrderedDict` and an alignment shift dictionary (`{old index: new index, ..}`),
        reset the keys of the fitness `OrderedDict`
        """
        
    def get_alignment(self, fitness_seq, complex_seq):
        """
        Aligns two AA sequences with BioPython's `PairwiseAligner` (https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec128).
        We do a local alignment with BLOSUM to take evolutionary divergence into account. 
        """

        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -10 # set these to make gaps happen less. With fitness data we know there shouldn't
        aligner.extend_gap_score = -0.5 # really be any gaps.
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62") # SOTA alignment matrix
        alignments = aligner.align(fitness_seq, complex_seq) # produces a generator

        if len(alignments) > 1:
            raise NotImplementedError(f"More than 1 alignments found, this is currently not implemented.")
        logging.info(f"Found alignment:\n{alignments[0]}")

        return alignments[0]

    def align_fitness(self):
        """
        [Placeholder] align a fitness OrderedDict to a complex object
        """
        logger.info("Aligning fitness sequence to complex..\n")
        alignment = self.get_alignment(self.fitness_get_seq(), self.complex_get_seq())
        

if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_TRUNCATED
    from choppa.data.toy_data.resources import TOY_FITNESS_DATA_SECTIONED, TOY_FITNESS_DATA_COMPLETE_NOCONF

    from choppa.IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED, 
                                    # confidence_colname="confidence"
                                    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()

    AlignFactory(fitness_dict, complex).align_fitness()
    # print(type(fitness_dict), type(complex))
    # AlignFactory(fitness_dict, complex)

    
    