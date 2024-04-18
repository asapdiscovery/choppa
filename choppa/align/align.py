
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

    def complex_get_seqidcs(self):
        """
        From a `PDB.Structure.Structure`, extracts the sequence indices as stored in the PDB file
        """
        return [ res.get_id()[1] for res in self.complex.get_residues() ]      

    def get_fitness_alignment_shift_dict(self, alignment):
        """
        Given an input complex sequence with residue indices (may not start at 0) and the fitness-complex
        alignment, creates a dictionary with indices that should be used for the fitness data of the form
        {fitness_idx : aligned_idx}     
        """
        
        alignment_shift_dict = {}
        for i, (fitness_res, pdb_res, pdb_resid) in enumerate(zip(alignment[0], alignment[1], self.complex_get_seqidcs()), start=1):
            # do some checks before adding to the alignment dict. we do these checks at multiple layers to be 100% sure we're not mismatching the two sequences.
            if fitness_res == "-" and fitness_res != pdb_res:
                # the fitness data does not contain this residue in the PDB and alignment has created a gap -> good
                alignment_shift_dict[i] = pdb_resid
            elif fitness_res == pdb_res:
                # the fitness data does contain this residue in the PDB and alignment has matched it -> good
                alignment_shift_dict[i] = pdb_resid
            else:
                # either the fitness residue is mismatched to another residue type in the PDB, or the fitness residue
                # is not in the PDB and the alignment has suggested a gap in the PDB -> both bad
                raise ValueError(f"Alignment has matched a fitness residue to either a gap in the PDB or to a different residue type:\nPDB index: {pdb_resid}\nPDB residue: {pdb_res}\nFitness index: {i}\nFitness residue: {fitness_res}\n.. consider adjusting alignment hyperparameters.")

        return alignment_shift_dict
    
    def fitness_reset_keys(self, alignment):
        """
        Given a fitness `OrderedDict` and an `alignment`,
        reset the keys of the fitness `OrderedDict`. 

        NB: also fills the fitness dict with indices that exist in the PDB but not in the fitness data, i.e.
        represented as 'empty' dict entries. This way the fitness HTML view will have 'empty' fitness data
        for those residues.
        """
        alignment_dict = self.get_fitness_alignment_shift_dict(alignment)
        reset_dict = {}
        for _, fitness_data in self.fitness_input.items():
            # we build a new dict where keys are the aligned index, then the aligned/unaligned indices (provenance),
            # then wildtype data and then per-mutant fitness data

            if not fitness_data['fitness_csv_index'] in alignment_dict.keys():
                logger.warn(f"Fitness data found to have a residue (index {fitness_data['fitness_csv_index']}) not in the PDB - skipping.")
                continue

            reset_dict[alignment_dict[fitness_data['fitness_csv_index']]] = \
            {'fitness_aligned_index' : alignment_dict[fitness_data['fitness_csv_index']], **fitness_data}
        
        return reset_dict
    
    def fill_aligned_fitness(self, aligned_fitness_dict):
        """
        For an aligned fitness dict, there may be gaps with respect to the complex PDB. Fills these with empty data
        for easier parsing during visualization.
        """
        filled_aligned_fitness_dict = {}
        
        for complex_idx, complex_res in zip(self.complex_get_seqidcs(), self.complex_get_seq()):
            if complex_res == "X":
                continue # this is a ligand, we can skip because we don't show fitness for this anyway
            
            if complex_idx in aligned_fitness_dict:
                # this fitness-complex data matches, can just copy the data across
                filled_aligned_fitness_dict[complex_idx] = aligned_fitness_dict[complex_idx]
            else:
                # no fitness data for this residue in the complex, need to make an empty one
                filled_aligned_fitness_dict[complex_idx] = {'wildtype':{complex_res}}

        return filled_aligned_fitness_dict, len(filled_aligned_fitness_dict)-len(aligned_fitness_dict)

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
            # can implement a fix once we hit a system that has this; not sure how to handle other than take the top alignment
            raise NotImplementedError(f"More than 1 alignments found, this is currently not implemented.")
        logging.info(f"Found alignment:\n{str(alignments[0]).replace('target', 'CSV   ').replace('query', 'PDB  ')}")

        return alignments[0]

    def align_fitness(self):
        """
        Align a fitness OrderedDict to a complex object
        """
        logger.info("Aligning fitness sequence to complex..\n")
        alignment = self.get_alignment(self.fitness_get_seq(), self.complex_get_seq())

        aligned_fitness = self.fitness_reset_keys(alignment)

        filled_aligned_fitness_dict, num_filled = self.fill_aligned_fitness(aligned_fitness)
        logger.info(f"After aligning fitness data to PDB complex, filled {num_filled} empty entries in the fitness sequence (total entries in sequence: {len(self.complex_get_seq())}).\n")
        
        return filled_aligned_fitness_dict

if __name__ == "__main__":
    from choppa.data.toy_data.resources import TOY_COMPLEX, TOY_FITNESS_DATA_COMPLETE, TOY_FITNESS_DATA_TRUNCATED
    from choppa.data.toy_data.resources import TOY_FITNESS_DATA_SECTIONED, TOY_FITNESS_DATA_COMPLETE_NOCONF

    from choppa.IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(TOY_FITNESS_DATA_SECTIONED, 
                                    # confidence_colname="confidence"
                                    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()



    
    