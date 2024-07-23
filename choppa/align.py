import logging, sys
from collections import OrderedDict
from Bio import PDB, SeqUtils, Align
from Bio.Align import substitution_matrices

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger()


class AlignFactory:
    """
    Base class for aligning Fitness data with PDB complexes within `choppa`.
    """

    def __init__(
        self,
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
        seq_list = [
            fitness_values["wildtype"]["aa"]
            for _, fitness_values in self.fitness_input.items()
        ]

        return "".join(seq_list)

    def fitness_get_seqidcs(self):
        """
        From a fitness `OrderedDict`, extracts the amino acid sequence indices as a list
        """
        # could just do this with range but keep in this form if we ever need to revert to csv indices
        seq_idcs = [i for i, (_, _) in enumerate(self.fitness_input.items(), start=1)]

        return seq_idcs

    def complex_get_seq(self):
        """
        From a `PDB.Structure.Structure`, extracts the amino acid sequence
        """
        seq_list = [
            SeqUtils.seq1(res.get_resname()) for res in self.complex.get_residues()
        ]

        return "".join(seq_list)

    def complex_get_seqidcs(self):
        """
        From a `PDB.Structure.Structure`, extracts the sequence indices as stored in the PDB file
        """
        return [res.get_id()[1] for res in self.complex.get_residues()]

    def get_fitness_alignment_shift_dict(self, alignment):
        """
        Given an input complex sequence with residue indices (may not start at 0) and the fitness-complex
        alignment, creates a dictionary with indices that should be used for the fitness data of the form
        {fitness_idx : aligned_idx}
        """
        alignment_shift_dict = {}
        to_remove = []
        incrementer = 0
        for fitness_res, fitness_resid, pdb_res, pdb_resid in zip(
            alignment[0],
            self.fitness_get_seqidcs(),
            alignment[1],
            self.complex_get_seqidcs(),
        ):
            # do some checks before adding to the alignment dict. we do these checks at multiple layers to be 100% sure we're not mismatching the two sequences.
            if fitness_res != pdb_res:
                    to_remove.append(fitness_resid)
                    # incrementer += 1 # we shift the alignment by 1.
                # if pdb_res.isalpha() and fitness_res == "-": 
                #     # the fitness data does not contain this residue in the PDB and alignment has created a gap -> good
                #     alignment_shift_dict[fitness_resid + incrementer] = pdb_resid
                #     print(f"{fitness_res}{fitness_resid}->{pdb_res}{pdb_resid}") ## DEBUG
                # elif pdb_res == "-" and fitness_res.isalpha():
                #     # the fitness residue was matched to no residue in the PDB - this is OK (assuming the alignment is good)
                #     incrementer += 1 # we shift the alignment by 1.
                #     alignment_shift_dict[fitness_resid + incrementer] = pdb_resid
                #     print(f"{fitness_res}{fitness_resid}->{pdb_res}{pdb_resid}") ## DEBUG
                # elif pdb_res.isalpha() and fitness_res.isalpha():
                #     # this is slightly worrying because the alignment has been forced to match the wrong fitness residue with a PDB 
                #     # residue (i.e. a residue type mismatch). As long as the alignment is good, this should be sustained.
                #     logger.warn(
                #     f"Unable to match fitness residue {fitness_res} ({fitness_resid}) to PDB residue {pdb_res} ({pdb_resid})"
                #     )
                #     to_remove.append(fitness_resid)
                #     # incrementer += 1 # we shift the alignment by 1.
            else:
                #fitness_res == pdb_res:
                # the fitness data does contain this residue in the PDB and alignment has matched it -> good
                alignment_shift_dict[fitness_resid + incrementer] = pdb_resid
                print(f"{fitness_res}{fitness_resid}->{pdb_res}{pdb_resid}") ## DEBUG
            
            # else:
            #     raise ValueError(f"Encountered irregular alignment case for fitness {fitness_res}{fitness_resid} and PDB {pdb_res}{pdb_resid}. Please add intended behavior to conditionals in code blocks above this error catcher.")
        
        return alignment_shift_dict, to_remove

    def fitness_reset_keys(self, alignment):
        """
        Given a fitness `OrderedDict` and an `alignment`,
        reset the keys of the fitness `OrderedDict`.

        NB: also fills the fitness dict with indices that exist in the PDB but not in the fitness data, i.e.
        represented as 'empty' dict entries. This way the fitness HTML view will have 'empty' fitness data
        for those residues.
        """
        alignment_dict, to_remove = self.get_fitness_alignment_shift_dict(alignment)
        reset_dict = {}
        for _, fitness_data in self.fitness_input.items():
            # we build a new dict where keys are the aligned index, then the aligned/unaligned indices (provenance),
            # then wildtype data and then per-mutant fitness data

            if not fitness_data["fitness_csv_index"] in alignment_dict.keys():
                logger.warn(
                    f"Fitness data found to have a residue (index {fitness_data['fitness_csv_index']}) not in the PDB - skipping."
                )
                continue

            reset_dict[alignment_dict[fitness_data["fitness_csv_index"]]] = {
                "fitness_aligned_index": alignment_dict[
                    fitness_data["fitness_csv_index"]
                ],
                **fitness_data,
            }

        return reset_dict, to_remove

    def fill_aligned_fitness(self, aligned_fitness_dict, to_remove):
        """
        For an aligned fitness dict, there may be gaps with respect to the complex PDB. Fills these with empty data
        for easier parsing during visualization.
        """
        filled_aligned_fitness_dict = {}
        
        for complex_idx, complex_res in zip(
            self.complex_get_seqidcs(), self.complex_get_seq()
        ):
            if complex_res == "X":
                continue  # this is a ligand or water, we can skip because we don't show fitness for this anyway
            if complex_idx in to_remove:
                # wrong fitness data for this residue in the complex, need to make an empty one
                filled_aligned_fitness_dict[complex_idx] = {"wildtype": {complex_res}}
                continue
            


            #### now that we have reset both protein PDB and fitness sequence, this fucks up the below check. 
            # how to fix this? how can we assign the correct fitness datapoint to the correct protein index?


            if complex_idx in aligned_fitness_dict:
                # check that the fitness wildtype equals the protein PDB residue type
                if not self.complex_get_seq()[complex_idx] == aligned_fitness_dict[complex_idx]['wildtype']['aa']:
                    # hard stop - this is a critical alignment mismatch
                    raise ValueError(f"Alignment mismatch between wildtype and PDB!\n\nFitness: {aligned_fitness_dict[
                    complex_idx
                ]}\n\nProtein: {complex_idx}{complex_res}")
                # this fitness-complex data matches, can just copy the data across
                filled_aligned_fitness_dict[complex_idx] = aligned_fitness_dict[
                    complex_idx
                ]
            else:
                # no fitness data for this residue in the complex, need to make an empty one
                filled_aligned_fitness_dict[complex_idx] = {"wildtype": {complex_res}}

        return filled_aligned_fitness_dict, len(filled_aligned_fitness_dict) - len(
            aligned_fitness_dict
        )

    def get_alignment(self, fitness_seq, complex_seq):
        """
        Aligns two AA sequences with BioPython's `PairwiseAligner` (https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec128).
        We do a local alignment with BLOSUM to take evolutionary divergence into account.
        """

        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = (
            -10
        )  # set these to make gaps happen less. With fitness data we know there shouldn't
        aligner.extend_gap_score = -0.5  # really be any gaps.
        aligner.substitution_matrix = substitution_matrices.load(
            "BLOSUM62"
        )  # SOTA alignment matrix
        alignments = aligner.align(fitness_seq, complex_seq)  # produces a generator

        # if len(alignments) > 1:
        #     # can implement a fix once we hit a system that has this; not sure how to handle other than take the top alignment
        #     raise NotImplementedError(f"More than 1 alignments found, this is currently not implemented.")
        logging.info(
            f"Found alignment:\n{str(alignments[0]).replace('target', 'CSV   ').replace('query', 'PDB  ')}"
        )

        return alignments[0]

    def align_fitness(self):
        """
        Align a fitness OrderedDict to a complex object
        """
        logger.info("Aligning fitness sequence to complex..\n")
        alignment = self.get_alignment(self.fitness_get_seq(), self.complex_get_seq())
        
        aligned_fitness, to_remove = self.fitness_reset_keys(alignment)

        filled_aligned_fitness_dict, num_filled = self.fill_aligned_fitness(
            aligned_fitness, to_remove
        )
        logger.info(
            f"After aligning fitness data to PDB complex, filled {num_filled} empty entries in the fitness sequence (total entries in sequence: {len(self.complex_get_seq())}).\n"
        )

        return filled_aligned_fitness_dict


if __name__ == "__main__":
    from data.toy_data.resources import (
        TOY_COMPLEX,
        TOY_FITNESS_DATA_COMPLETE,
        TOY_FITNESS_DATA_TRUNCATED,
    )
    from data.toy_data.resources import (
        TOY_FITNESS_DATA_SECTIONED,
        TOY_FITNESS_DATA_COMPLETE_NOCONF,
    )

    from IO.input import FitnessFactory, ComplexFactory

    fitness_dict = FitnessFactory(
        TOY_FITNESS_DATA_SECTIONED,
        # confidence_colname="confidence"
    ).get_fitness_basedict()
    complex = ComplexFactory(TOY_COMPLEX).load_pdb()

    filled_aligned_fitness_dict = AlignFactory(fitness_dict, complex).align_fitness()
