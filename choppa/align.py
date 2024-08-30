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

    def alignment_idx_to_original_idx(self, alignment):
        """
        The BioPython alignment shifts indices freely based on query/target overlap. This function
        return a dict that maps the new indices in the alignment back to the old/original index sequence.
        """
        aligned_to_idx = {}
        counter = 0
        for i, res in enumerate(alignment):
            if res == "-":
                aligned_to_idx[i] = None
            else:
                aligned_to_idx[i] = counter
                counter += 1

        return aligned_to_idx

    def get_fitness_alignment_shift_dict(self, alignment):
        """
        Given an input complex sequence with residue indices (may not start at 0) and the fitness-complex
        alignment, creates a dictionary with indices that should be used for the fitness data of the form
        {fitness_idx : aligned_idx}
        """
        fitness_idx_dict = self.alignment_idx_to_original_idx(alignment[0])
        complex_idx_dict = self.alignment_idx_to_original_idx(alignment[1])

        alignment_shift_dict = {}
        for i, (ali_fitness_res, ali_complex_res) in enumerate(
            zip(
                alignment[0],
                alignment[1],
            )
        ):
            if fitness_idx_dict[i] and complex_idx_dict[i]:
                
                # this means that there is an alignment on this index.
                # get the original fitness/complex residue types on this index so we can double-check
                ori_fitness_res = self.fitness_get_seq()[fitness_idx_dict[i]]
                ori_complex_res = self.complex_get_seq()[complex_idx_dict[i]]
                if ali_fitness_res == ali_complex_res:
                    # there is a residue match
                    if (  # these should always all be the same - otherwise the final fitness view will have mismatched fitness per residue
                        not ali_fitness_res
                        == ali_complex_res
                        == ori_fitness_res
                        == ori_complex_res
                    ):
                        raise ValueError(
                            "Alignment failed; unable to match aligned sequence indices back to original sequence indices - check your alignment."
                        )
                    # okay, so for the aligned sequences we need to map the ORIGINAL fitness sequence indices
                    # to the ORIGINAL complex sequence indices using the alignment. Because we've grabbed the
                    # ORIGINAL indices using self.alignment_idx_to_original_idx() we know these are correct.
                    alignment_shift_dict[
                        self.fitness_get_seqidcs()[fitness_idx_dict[i]]
                    ] = self.complex_get_seqidcs()[complex_idx_dict[i]]
                else:
                    # either a mismatch (point mutation) or ligand. Can just not add to alignment dict, will
                    # be ignored downstream and result in a black residue (no fitness data) in case there is
                    # a complex residue.
                    pass
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
        print(alignment_dict)
        raise Exception
        reset_dict = {}
        for _, fitness_data in self.fitness_input.items():
            # we build a new dict where keys are the aligned index, then the aligned/unaligned indices (provenance),
            # then wildtype data and then per-mutant fitness data
            if not fitness_data["fitness_csv_index"] in alignment_dict.keys():
                # logger.warn( # bit too chatty - disable for now.
                #     f"Fitness data found to have a residue (index {fitness_data['fitness_csv_index']}) not in the PDB - skipping."
                # )
                continue

            reset_dict[alignment_dict[fitness_data["fitness_csv_index"]]] = {
                "fitness_aligned_index": alignment_dict[
                    fitness_data["fitness_csv_index"]
                ],
                **fitness_data,
            }

        return reset_dict

    def fill_aligned_fitness(self, aligned_fitness_dict):
        """
        For an aligned fitness dict, there may be gaps with respect to the complex PDB. Fills these with empty data
        for easier parsing during visualization.
        """
        filled_aligned_fitness_dict = {}

        for complex_idx, complex_res in zip(
            self.complex_get_seqidcs(), self.complex_get_seq()
        ):
            if complex_res == "X":
                continue  # this is a ligand, we can skip because we don't show fitness for this anyway

            if complex_idx in aligned_fitness_dict:
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
        # dump the seqs
        # print("FITNESS_SEQ")
        # for i, f in enumerate(fitness_seq):
        #     print(i, f)
        # print("COMPLEX_SEQ")
        # for i, f in enumerate(complex_seq):
        #     print(i, f)
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"  # this means that both alignments (fitness/PDB) start at index 0. Easier for downstream.
        aligner.open_gap_score = (
            -20
        )  # set these to make gaps happen less. With fitness data we know there shouldn't really be any gaps.
        aligner.extend_gap_score = 2
        aligner.substitution_matrix = substitution_matrices.load(
            "PAM70"
        )  # found this algorithm by comparing all available in BioPython
        # produces a generator but we can just pick the top one
        alignment = aligner.align(fitness_seq, complex_seq)[0]

        # if len(alignments) > 1:
        #     # can implement a fix once we hit a system that has this; not sure how to handle other than take the top alignment
        #     raise NotImplementedError(f"More than 1 alignments found, this is currently not implemented.")
        logging.info(
            f"Found alignment:\n{str(alignment).replace('target', 'CSV   ').replace('query', 'PDB  ')}"
        )

        return alignment

    def align_fitness(self):
        """
        Align a fitness OrderedDict to a complex object
        """
        logger.info("Aligning fitness sequence to complex..\n")
        alignment = self.get_alignment(self.fitness_get_seq(), self.complex_get_seq())

        aligned_fitness = self.reset_keys2(alignment)

        # filled_aligned_fitness_dict, num_filled = self.fill_aligned_fitness(
        #     aligned_fitness
        # )
        # logger.info(
        #     f"After aligning fitness data to PDB complex, filled {num_filled} empty entries in the fitness sequence (total entries in sequence: {len(self.complex_get_seq())}).\n"
        # )

        return aligned_fitness


    def reset_keys2(self, alignment):
        # print(type(alignment))
        import numpy as np
        fitness_indices = np.asarray(alignment.indices[0])
        # print(fitness_indices)
        # print(np.asarray(alignment.inverse_indices[0]))

        # print(len(fitness_indices))
        complex_indices = np.asarray(alignment.indices[1])
        # print(complex_indices)
        print(np.asarray(alignment.inverse_indices[1]))
        inv_i = np.asarray(alignment.inverse_indices[1])
        translation_dict = {}
        for i, (fitness_idx, complex_idx) in enumerate(zip(fitness_indices, complex_indices)):
            if fitness_idx != -1 and complex_idx != -1 and alignment[1][i] != "X":
                translation_dict[fitness_idx] = complex_idx
        
        inv_translation_dict = {v: k for k, v in translation_dict.items()}
        
        print(inv_translation_dict)
        print(translation_dict)
        data = self.fitness_input.copy()

        newdata = {}
        for k, v in data.items():
            if k in translation_dict:
                print(v["wildtype"]["aa"])
                newdata[translation_dict[k]] = v
                print(translation_dict[k], v)
                         
        


        return newdata
            

            





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
