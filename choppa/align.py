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
        {fitness_idx : aligned_idx}.

        This is complicated for multiple reasons, so we're iterating over each individual index. For example
        in the following alignment:
        CSV              50 DMYIERAGDITWEKDAEVTGNSPRLDVALDESGDFSLVEEDGPPMREIILKVVLMAICGM
                          0 |||||||||||||||||||||||||||||||||||||||---------------------
        PDB               0 DMYIERAGDITWEKDAEVTGNSPRLDVALDESGDFSLVE---------------------

        CSV             110 NPIAIPFAAGAWYVYVKTGKRSGALWDVPAPKEVKKGETTDGVYRVMTRRLLGSTQVGVG
                         60 ------------------------------------||||||||||||||||||||||||
        PDB              39 ------------------------------------GETTDGVYRVMTRRLLGSTQVGVG

        we need to 1) keep track of the starting indices of fitness ('CSV') and crystal structure ('PDB') (50 and 0, resp.)
        and 2) we need to be able to skip the gap in alignment.
        """

        alignment_shift_dict = {}
        to_remove = []

        # first we grab the starting indices for fitness ('CSV') and complex ('PDB') by slicing
        # the alignment object view. Hacky but robust, no method implemented in BioPython for this.
        start_idx_fitness = int(alignment.format().splitlines()[0].split()[1])
        start_idx_complex = int(alignment.format().splitlines()[2].split()[1])
        # print(f"Start fitness: {start_idx_fitness}, start complex: {start_idx_complex}") # DEBUG

        # now we will loop over the alignment. We need both the fitness and complex residues
        # and the original indices.

        # this might break if the PDB is longer than the fitness data?
        for fitness_res, complex_res in zip(alignment[0], alignment[1]):
            # print(start_idx_fitness, fitness_res, complex_res, start_idx_complex) # DEBUG
            if fitness_res == complex_res:
                # good match. Can add this fitness data to the dict. Can bump both.
                alignment_shift_dict[start_idx_fitness] = start_idx_complex
                start_idx_fitness += 1
                start_idx_complex += 1
            else:
                # bad match.
                to_remove.append(start_idx_complex)
                if complex_res == "-":
                    # there is a gap in the fitness data.
                    # skip over this fitness datapoint only
                    start_idx_fitness += 1
                else:
                    # the alignment matched a fitness residue to the wrong residue type, can happen in e.g. point mutations
                    # in this case we also need to skip over the protein residue
                    start_idx_complex += 1
                    start_idx_fitness += 1

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
                # logger.warn( # disabled for now, can spam a lot
                #     f"Fitness data found to have a residue (index {fitness_data['fitness_csv_index']}) not in the PDB - skipping."
                # )
                continue

            reset_dict[alignment_dict[fitness_data["fitness_csv_index"]]] = {
                "fitness_aligned_index": alignment_dict[
                    fitness_data["fitness_csv_index"]
                ],
                **fitness_data,
            }

        """
        are we just setting the wrong indexing somewhere? 
        looks like we might be taking the fitness residue 
        instead of complex?
        """
        print(to_remove)
        print(reset_dict.keys())
        for resi_to_remove in set(to_remove):
            # remove complex indices that we have no fitness data for. We'll fill these with empty fitness data
            # later on.
            if resi_to_remove in reset_dict:
                reset_dict.pop(resi_to_remove)

        return reset_dict

    def fill_aligned_fitness(self, aligned_fitness_dict):
        """
        For an aligned fitness dict, there may be gaps with respect to the complex PDB. Fills these with empty data
        for easier parsing during visualization.
        """
        filled_aligned_fitness_dict = {}
        print(aligned_fitness_dict.keys())
        for complex_idx, complex_res in zip(
            self.complex_get_seqidcs(), self.complex_get_seq()
        ):
            if complex_res == "X":
                continue  # this is a ligand or water, we can skip because we don't show fitness for this anyway

            if complex_idx not in aligned_fitness_dict:
                # no fitness data for this residue in the complex, need to make an empty one
                filled_aligned_fitness_dict[complex_idx] = {"wildtype": {complex_res}}
            elif complex_idx in aligned_fitness_dict:
                print(
                    f"{complex_idx}:{self.complex_get_seq()[complex_idx]} should be {aligned_fitness_dict[complex_idx]['wildtype']['aa']}:{aligned_fitness_dict[complex_idx]['fitness_csv_index']}"
                )  # DEBUG
                # check that the fitness wildtype equals the protein PDB residue type
                if (
                    not self.complex_get_seq()[complex_idx]
                    == aligned_fitness_dict[complex_idx]["wildtype"]["aa"]
                ):
                    # hard stop - this is a critical alignment mismatch
                    raise ValueError(
                        f"Alignment mismatch between wildtype and PDB!\n\nFitness: {aligned_fitness_dict[complex_idx]}\n\nProtein: {complex_idx}{complex_res}"
                    )
                # this fitness-complex data matches, can just copy the data across
                filled_aligned_fitness_dict[complex_idx] = aligned_fitness_dict[
                    complex_idx
                ]
            else:
                # no fitness data for this residue in the complex, need to make an empty one
                filled_aligned_fitness_dict[complex_idx] = {"wildtype": {complex_res}}

        for i, j in filled_aligned_fitness_dict.items():
            print()
            print(i, j)
            break
        # so alignment indices are correct now, but for some reason the wrong logoplots are showing up?
        # do the wildtypes in the filled_aligned_fitness_dict correspond? why
        # are the surface colors wrong??

        # for some reason the indices in filled_aligned_fitness_dict are fucked, fitness_csv_index should start at 50

        return filled_aligned_fitness_dict, len(filled_aligned_fitness_dict) - len(
            aligned_fitness_dict
        )

    def get_alignment(self, fitness_seq, complex_seq):
        """
        Aligns two AA sequences with BioPython's `PairwiseAligner` (https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec128).
        We do a local alignment with BLOSUM to take evolutionary divergence into account.
        """

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.open_gap_score = (
            -20
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

        aligned_fitness = self.fitness_reset_keys(alignment)

        filled_aligned_fitness_dict, num_filled = self.fill_aligned_fitness(
            aligned_fitness
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
